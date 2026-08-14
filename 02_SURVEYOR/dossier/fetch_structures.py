"""DOSSIER — pre-fetch + render canonical AlphaFold structure images.

For visual purposes only in the dossier header (not used for docking/
structural analysis -- that's master_surveyor's job, see
docs/MASTER_SURVEYOR_plan.md's m2_structures.py). AlphaFold DB's API
returns coordinate files (PDB/CIF) and a PAE plot, but no pre-rendered
cartoon image, so this renders one itself: fetch the canonical AlphaFold
model for the hit's uniprot_acc, load it in 3Dmol.js inside a headless
Chrome tab (real WebGL via the swiftshader software rasterizer -- plain
--headless does NOT expose WebGL without the angle/swiftshader flags),
screenshot, crop the whitespace margin, cache as one PNG per accession.

Cached by uniprot_acc, not by hit -- many hits share a gene, and the
canonical structure only depends on the gene, not the cell type or which
alt transcript triggered the hit. render.py's header_section() just
checks this cache and silently omits the image if the accession isn't
there (fetch failed, no AlphaFold entry, or this script hasn't run yet) --
a dossier never depends on this cache existing to render correctly.

Usage:
    python fetch_structures.py
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import requests
from PIL import Image

from dossier.config import AFDB_API, AFDB_TIMEOUT, ASSETS_DIR, CHROME_BIN, MASTER_SURVEYOR_GROUPS, STRUCTURE_CACHE_DIR
from dossier.data import load_hits

_CHROME_FLAGS = [
    "--headless", "--disable-gpu", "--no-sandbox",
    "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader",
    "--window-size=620,620", "--virtual-time-budget=4000",
]


def _log(msg: str) -> None:
    print(f"  [structures] {msg}", flush=True)


def _fetch_pdb_text(uniprot_acc: str) -> str | None:
    try:
        resp = requests.get(AFDB_API.format(acc=uniprot_acc), timeout=AFDB_TIMEOUT)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        entries = resp.json()
        if not entries:
            return None
        pdb_resp = requests.get(entries[0]["pdbUrl"], timeout=AFDB_TIMEOUT * 2)
        pdb_resp.raise_for_status()
        return pdb_resp.text
    except requests.RequestException as e:
        _log(f"{uniprot_acc}: fetch failed ({e})")
        return None


def _crop_whitespace(src: Path, dst: Path, pad: int = 14) -> None:
    img = Image.open(src).convert("RGB")
    arr = np.array(img)
    non_white = ~np.all(arr > 250, axis=2)
    if not non_white.any():
        img.save(dst)
        return
    ys, xs = non_white.any(axis=1), non_white.any(axis=0)
    y0, y1 = int(ys.argmax()), int(len(ys) - ys[::-1].argmax())
    x0, x1 = int(xs.argmax()), int(len(xs) - xs[::-1].argmax())
    box = (max(x0 - pad, 0), max(y0 - pad, 0), min(x1 + pad, img.width), min(y1 + pad, img.height))
    img.crop(box).save(dst)


def _render_png(pdb_text: str, out_path: Path) -> bool:
    js_lib = (ASSETS_DIR / "js" / "3Dmol-min.js").read_text()
    html = f"""<!doctype html><html><head><meta charset="utf-8"/>
<style>html,body{{margin:0;padding:0;background:#ffffff;}}#viewer{{width:600px;height:600px;}}</style>
</head><body>
<div id="viewer"></div>
<script>{js_lib}</script>
<script>
  var viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{backgroundColor: 'white'}});
  viewer.addModel({json.dumps(pdb_text)}, 'pdb');
  viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
  viewer.zoomTo();
  viewer.render();
</script>
</body></html>"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".html", delete=False) as f:
        f.write(html)
        html_path = Path(f.name)
    raw_path = out_path.with_name(out_path.stem + ".raw.png")
    try:
        subprocess.run(
            [CHROME_BIN, *_CHROME_FLAGS, f"--screenshot={raw_path}", f"file://{html_path}"],
            check=True, capture_output=True, timeout=60,
        )
        if not raw_path.exists():
            return False
        _crop_whitespace(raw_path, out_path)
        return True
    except (subprocess.SubprocessError, OSError) as e:
        _log(f"render failed: {e}")
        return False
    finally:
        html_path.unlink(missing_ok=True)
        raw_path.unlink(missing_ok=True)


def main() -> None:
    STRUCTURE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    df = load_hits()
    # Scoped to master_surveyor's own hit set (docs/MASTER_SURVEYOR_plan.md),
    # same boundary generate_all.py uses -- re-run with a broader df slice
    # later if a dossier outside this scope wants a cached image too; the
    # cache key is uniprot_acc, so nothing here is scope-specific by design.
    scoped = df[df["master_group"].isin(MASTER_SURVEYOR_GROUPS)]
    accessions = sorted(scoped["uniprot_acc"].dropna().unique())
    _log(f"{len(accessions):,} unique UniProt accessions in master_surveyor's {len(scoped):,}-row scope")

    rendered, cached, failed = 0, 0, []
    for i, acc in enumerate(accessions, 1):
        out_path = STRUCTURE_CACHE_DIR / f"{acc}.png"
        if out_path.exists():
            cached += 1
            continue
        pdb_text = _fetch_pdb_text(acc)
        if pdb_text is None:
            failed.append(acc)
        elif _render_png(pdb_text, out_path):
            rendered += 1
        else:
            failed.append(acc)
        if i % 10 == 0 or i == len(accessions):
            _log(f"{i:,}/{len(accessions):,} processed "
                 f"({rendered} rendered, {cached} already cached, {len(failed)} failed)")

    _log(f"done: {rendered} rendered, {cached} already cached, {len(failed)} failed")
    if failed:
        _log(f"failed accessions: {failed}")


if __name__ == "__main__":
    main()

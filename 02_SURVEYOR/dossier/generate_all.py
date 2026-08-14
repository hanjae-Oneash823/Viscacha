"""DOSSIER — batch-generate every hit that enters master_surveyor.

master_surveyor's own scope (docs/MASTER_SURVEYOR_plan.md) is
trial_failure_candidate + drug_repurposing_candidate only -- novel_target_candidate
is explicitly excluded there, so it's excluded here too.

True unique hit key is (gene_name, cell_type, hit_ENST_ID), NOT just
(gene_name, cell_type): 16 pairs carry two distinct hits (e.g. a
trial_failure hit and a separate drug_repurposing hit on different
transcripts of the same gene/cell_type). Output filenames disambiguate only
those 16 -- every other file keeps the plain <gene>_<cell_type>.html name.

Usage:
    python generate_all.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier.config import MASTER_SURVEYOR_GROUPS, OUT_DIR
from dossier.data import dossier_filename, load_hits
from dossier.render import render


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    hits = load_hits()
    hits = hits[hits["master_group"].isin(MASTER_SURVEYOR_GROUPS)]
    keys = hits[["gene_name", "cell_type", "hit_ENST_ID", "master_group"]].drop_duplicates()
    keys = keys.sort_values(["gene_name", "cell_type"])

    n = len(keys)
    print(f"{n:,} hits to render "
          f"({(keys['master_group'] == 'trial_failure_candidate').sum():,} trial_failure + "
          f"{(keys['master_group'] == 'drug_repurposing_candidate').sum():,} drug_repurposing)")

    ok, failed = 0, []
    t0 = time.time()
    for i, (_, r) in enumerate(keys.iterrows(), 1):
        gene, ct, hit_enst = r["gene_name"], r["cell_type"], r["hit_ENST_ID"]
        out_path = OUT_DIR / dossier_filename(gene, ct, hit_enst)
        try:
            out_path.write_text(render(gene, ct, hit_enst))
            ok += 1
        except Exception as e:
            failed.append((gene, ct, hit_enst, repr(e)))
        if i % 20 == 0 or i == n:
            elapsed = time.time() - t0
            print(f"  {i:,}/{n:,} done ({elapsed:.0f}s, {len(failed)} failed so far)")

    print(f"\n{ok:,}/{n:,} dossiers written to {OUT_DIR}")
    if failed:
        print(f"{len(failed):,} failed:")
        for gene, ct, hit_enst, err in failed:
            print(f"  {gene} / {ct} / {hit_enst}: {err}")


if __name__ == "__main__":
    main()

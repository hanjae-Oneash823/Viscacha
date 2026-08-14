"""M3 — per-hit BIOVIA-Discovery-Studio-ready export folder.

Bundles M1's ligand SDF, M2's canonical + alt structures, and M2b's
structure-confidence read into one folder per hit, plus a manifest.json
carrying everything a user would want to know before opening the files in
Discovery Studio (including the cutoffs in effect at export time, for
reproducibility -- matches this codebase's existing audit-log philosophy).
"""

from __future__ import annotations

import json
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor.config import DOCKING_DIR


def _json_safe(obj):
    """Recursively replace NaN/Infinity with None. json.dumps happily emits
    the literal tokens NaN/Infinity for these (Python-specific leniency),
    but that is not valid JSON -- a standards-compliant parser (or another
    language's tooling) reading manifest.json would fail on it. hit_row
    fields like alt_usage_delta are legitimately NaN for hits where that
    field doesn't apply (e.g. drug_repurposing rows), so this can't just be
    avoided upstream -- it has to be sanitized at serialization time.
    """
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_safe(v) for v in obj]
    return obj


def _log(msg: str) -> None:
    print(f"  [m3] {msg}", file=sys.stderr, flush=True)


def hit_folder_name(gene: str, cell_type: str, transcript: str) -> str:
    safe_transcript = transcript.replace("/", "_")
    return f"{gene}_{cell_type}_{safe_transcript}"


def export_hit(
    hit_row: dict,
    canonical_structure_path: Path,
    canonical_structure_source: str,
    alt_top_model_path: Path,
    structure_confidence: dict,
    ligand_statuses: list[dict],
    ligands_sdf_path: Path,
    thresholds_used: dict,
) -> Path:
    """Returns the created folder's path.

    hit_row is a plain dict (a hits_deep.csv row, already resolved to one
    hit via m0_select.representative_row) -- the manifest fields below are
    exactly the set docs/MASTER_SURVEYOR_plan.md specifies, plus M2b's
    structure_confidence which didn't exist when that plan was written.
    """
    gene = hit_row["gene_name"]
    cell_type = hit_row["cell_type"]
    transcript = hit_row.get("hit_transcript_name") or hit_row["hit_ENST_ID"]

    folder = DOCKING_DIR / hit_folder_name(gene, cell_type, transcript)
    folder.mkdir(parents=True, exist_ok=True)

    shutil.copy(canonical_structure_path, folder / "protein_canonical.pdb")
    shutil.copy(alt_top_model_path, folder / "protein_alt.pdb")
    if ligands_sdf_path.exists():
        shutil.copy(ligands_sdf_path, folder / "ligands.sdf")

    manifest = {
        "gene_name": gene,
        "cell_type": cell_type,
        "transcript": transcript,
        "hit_ENST_ID": hit_row.get("hit_ENST_ID"),
        "master_group": hit_row.get("master_group"),
        "protein_change_type": hit_row.get("protein_change_type"),
        "affected_domain": hit_row.get("affected_domain"),
        "delta_usage": hit_row.get("alt_usage_delta"),
        "chi_padj": hit_row.get("chi_padj"),
        "uniprot_acc": hit_row.get("uniprot_acc"),
        "chembl_target_id": hit_row.get("chembl_target_id"),
        "canonical_structure_source": canonical_structure_source,
        "structure_confidence": structure_confidence,
        "ligands": ligand_statuses,
        "thresholds_used": thresholds_used,
    }
    (folder / "manifest.json").write_text(json.dumps(_json_safe(manifest), indent=2, default=str))

    n_ligands_ok = sum(1 for s in ligand_statuses if s.get("ok"))
    _log(f"exported {folder.name} "
         f"(structure_confidence={structure_confidence.get('verdict', 'unknown')}, "
         f"{n_ligands_ok}/{len(ligand_statuses)} ligands)")
    return folder


def zip_folder(folder: Path, out_zip: Path | None = None) -> Path:
    """Zips folder for download -- Discovery Studio runs on the user's own
    desktop, not this server, so exported files need to leave the box.
    """
    out_zip = out_zip or folder.with_suffix(".zip")
    archive_base = str(out_zip.with_suffix(""))
    shutil.make_archive(archive_base, "zip", root_dir=folder.parent, base_dir=folder.name)
    return Path(archive_base + ".zip")

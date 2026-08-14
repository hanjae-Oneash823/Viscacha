"""MASTER_SURVEYOR — batch CLI entry point.

Runs m0 (filter+collapse) -> m1 (ligands) -> m2 (structures) -> m2b (QC) ->
m3 (export) over every hit in scope, using default thresholds. Mirrors
junior_surveyor/run_junior_surveyor.py's orchestrator shape. This is the
headless/batch path; dossier_server drives the same m0-m3 functions
per-cart-item interactively instead of over the whole shortlist at once.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor import m0_select, m1_ligands, m2_structures, m2b_structure_qc, m3_export
from master_surveyor.config import DOCKING_DIR, SHORTLIST_CSV


def _log(msg: str) -> None:
    print(f"[run_master_surveyor] {msg}", file=sys.stderr, flush=True)


def run_one(hit_row: dict, thresholds_used: dict) -> Path | None:
    gene, cell_type = hit_row["gene_name"], hit_row["cell_type"]
    transcript = hit_row.get("hit_transcript_name") or hit_row["hit_ENST_ID"]
    _log(f"{gene} / {cell_type} / {transcript}")

    try:
        canonical = m2_structures.get_canonical_structure(
            hit_row.get("uniprot_acc", ""), hit_row["canonical_protein_seq"],
        )
        alt = m2_structures.get_alt_structure(hit_row["alt_protein_seq"])
        confidence = m2b_structure_qc.assess(
            alt, canonical.path,
            hit_row["canonical_protein_seq"], hit_row["alt_protein_seq"],
            int(hit_row.get("changed_aa_start") or 0), int(hit_row.get("changed_aa_end") or 0),
        )

        drug_names = [d for d in (hit_row.get("drug_names") or "").split("|") if d]
        ligands_sdf = DOCKING_DIR / m3_export.hit_folder_name(gene, cell_type, transcript) / "ligands.sdf"
        ligand_statuses = m1_ligands.build_ligands_sdf(drug_names, ligands_sdf) if drug_names else []

        return m3_export.export_hit(
            hit_row, canonical.path, canonical.source, alt.seed_models[0],
            confidence, ligand_statuses, ligands_sdf, thresholds_used,
        )
    except Exception as exc:  # noqa: BLE001 -- one hit's failure must not abort the batch
        _log(f"  FAILED: {exc}")
        return None


def main() -> None:
    shortlist = m0_select.run()
    SHORTLIST_CSV.parent.mkdir(parents=True, exist_ok=True)
    shortlist.to_csv(SHORTLIST_CSV, index=False)
    _log(f"{len(shortlist):,} hits in shortlist.csv")

    exported, failed = 0, 0
    for _, row in shortlist.iterrows():
        result = run_one(row.to_dict(), thresholds_used={})
        if result is not None:
            exported += 1
        else:
            failed += 1

    _log(f"done: {exported:,} exported, {failed:,} failed")


if __name__ == "__main__":
    main()

"""
Surveyor (Layer 2) — orchestrator (plan 4.1).

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python -m layer2.run_layer2
     (from /home/welcome3/Viscacha_pipeline)

Steps:
  1. Validate + load inputs; compute the Braak-condition correlation.
  2. Gene deduplication into gene / (gene, cell_type) work lists.
  3. Gene-level modules (M01 -> M02 -> M03, then M04, M05, M06) once per gene.
  4. Cell-type modules (M07) then M08 (dossier + report) per (gene, cell_type).
  5. Global outputs: run summary + audit log.

M04-M06 are run sequentially (correctness over speed); each external module is
individually non-blocking per the plan's failure policy, so one database being
down degrades a section rather than failing the run.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from layer2.config import CONFIG
from layer2.inputs import load_inputs, InputValidationError
from layer2.utils.api_client import APIClient
from layer2.utils.audit import AuditLog
from layer2 import m01_sequence as m01_mod
from layer2 import m02_domain as m02_mod
from layer2 import m03_splice as m03_mod
from layer2 import m04_drug as m04_mod
from layer2 import m05_disease as m05_mod
from layer2 import m06_pathway as m06_mod
from layer2 import m07_expression as m07_mod
from layer2 import m_vis as mvis_mod
from layer2 import m08_report as m08_mod
from layer2.m02_domain import M02Result


def _unavailable_m02(gene_id: str, uniprot, reason: str) -> M02Result:
    return M02Result(gene_id=gene_id, uniprot=uniprot, reviewed=False,
                     canonical_length=0, available=False, warnings=[reason])


def main():
    CONFIG.ensure_output_dirs()
    audit = AuditLog(CONFIG.audit_log_path, verbose=True)
    client = APIClient(audit)

    # --- Step 1-2: inputs ---
    try:
        inp = load_inputs()
    except InputValidationError as e:
        audit.error("inputs", str(e))
        audit.write()
        print(f"\nINPUT VALIDATION FAILED: {e}", file=sys.stderr)
        sys.exit(1)

    if inp.braak_condition_correlation is not None:
        audit.value("braak_condition_correlation", round(inp.braak_condition_correlation, 4))
    for w in inp.warnings:
        audit.warning("inputs", w)
    audit.info(f"{inp.n_genes} genes, {inp.n_candidates} candidates")

    # --- Step 3: gene-level modules ---
    gene_results: dict[str, tuple | None] = {}
    for gene in inp.gene_work_list:
        audit.info(f"--- gene {gene.gene_id} ---")
        try:
            m01 = m01_mod.run(gene, client, audit, inp.enst_map)
        except m01_mod.M01Error as e:
            audit.error("M01", f"{gene.gene_id}: {e}")
            gene_results[gene.gene_id] = None
            continue

        try:
            m02 = m02_mod.run(m01, client, audit)
        except m02_mod.M02Error as e:
            audit.error("M02", f"{gene.gene_id}: {e}")
            m02 = _unavailable_m02(gene.gene_id, m01.canonical_uniprot, str(e))

        m03 = m03_mod.run(m01, m02)
        m04 = m04_mod.run(gene, client, audit)
        m05 = m05_mod.run(gene, client, audit)
        m06 = m06_mod.run(gene, m01.canonical_uniprot, client, audit)
        gene_results[gene.gene_id] = (m01, m02, m03, m04, m05, m06)

    # --- Step 4: cell-type modules + report ---
    summary_rows = []
    for cand in inp.cell_type_work_list:
        gr = gene_results.get(cand.gene_id)
        if gr is None:
            summary_rows.append({"gene": cand.gene_id, "cell_type": cand.cell_type,
                                 "type": "FAILED", "tier": None, "verdict": "failed",
                                 "reason": "M01 sequence retrieval failed"})
            continue
        m01, m02, m03, m04, m05, m06 = gr
        m07 = m07_mod.run(cand, m01, client, audit)
        try:
            viz = mvis_mod.run(m07)
        except Exception as e:                       # viz must never break a dossier
            audit.warning("M_VIS", f"{cand.gene_id}/{cand.cell_type}: {e}")
            viz = {}
        ctx = m08_mod.CandidateContext(
            candidate=cand, m01=m01, m02=m02, m03=m03, m04=m04, m05=m05,
            m06=m06, m07=m07, braak_correlation=inp.braak_condition_correlation,
            input_checksums=inp.input_file_checksums,
            database_versions=audit.database_versions, viz_svgs=viz)
        summary_rows.append(m08_mod.run(ctx))

    # --- Step 5: global outputs ---
    _write_run_summary(summary_rows)
    audit.write()

    print("\n" + "=" * 64)
    print("Surveyor Layer 2 — COMPLETE")
    print(f"  dossiers: {CONFIG.candidates_dir}")
    print(f"  reports:  {CONFIG.reports_dir}")
    print(f"  summary:  {CONFIG.run_summary_path}")
    print("=" * 64)
    for r in summary_rows:
        print(f"  {r['gene']:9} {r['cell_type']:18} Type {r['type']:4} "
              f"tier {str(r['tier'] or '-'):2} -> {r['verdict'].upper():12} {r['reason']}")


def _write_run_summary(rows: list[dict]):
    lines = ["# Surveyor Layer 2 — Run Summary", "",
             "| Gene | Cell type | Type | Tier | Verdict | Rationale |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['gene']} | {r['cell_type']} | {r['type']} | "
                     f"{r['tier'] or '-'} | {r['verdict'].upper()} | {r['reason']} |")
    CONFIG.run_summary_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

"""
Viscacha Layer 0 — orchestrator
Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python -m layer0.run_layer0
(from /home/welcome3/Viscacha_pipeline)
"""

import sys
from pathlib import Path

# Ensure the project root is on sys.path when run directly
sys.path.insert(0, str(Path(__file__).parent.parent))

from layer0.config import OUT_META, OUT_ADATA, OUT_QC
from layer0.utils.qc_log import QCLogger
import layer0.step01_metadata         as step01
import layer0.step02_barcode_merge    as step02
import layer0.step03_barcode_filter   as step03
import layer0.step04_prevalence_filter as step04
import layer0.step05_covariate_audit  as step05
import layer0.step06_pseudobulk       as step06
import layer0.plots                   as plots


def main():
    qc_log = QCLogger()

    print("=" * 60)
    print("Viscacha Layer 0 — starting")
    print("=" * 60)

    # 0.1 — Metadata harmonization
    print("\n[Layer 0] Step 0.1: Metadata harmonization")
    unified_meta = step01.run(qc_log)
    OUT_META.mkdir(parents=True, exist_ok=True)
    unified_meta.to_csv(OUT_META / "unified_metadata.csv")
    print(f"  Saved unified_metadata.csv ({len(unified_meta)} donors)")
    plots.plot_step01(unified_meta)

    # 0.2 — Barcode merge
    print("\n[Layer 0] Step 0.2: Barcode merge and obs enrichment")
    adata_tx = step02.run(unified_meta, qc_log)

    # 0.3 — Barcode confidence filter
    print("\n[Layer 0] Step 0.3: Barcode confidence filter")
    adata_filtered = step03.run(adata_tx, qc_log)
    OUT_ADATA.mkdir(parents=True, exist_ok=True)
    adata_filtered.write_h5ad(OUT_ADATA / "adata_tx_step03.h5ad")
    print(f"  Saved adata_tx_step03.h5ad ({adata_filtered.n_obs} cells)")
    plots.plot_step03(adata_tx, adata_filtered)

    # Free unfiltered copy
    del adata_tx

    # 0.4 — Prevalence filter (diagnostic only — no longer applied to the pseudobulk
    # output; gene/transcript-count filtering moved entirely to satuRn's "use all
    # transcripts" model, so this is kept just to report what a prevalence cut would
    # have dropped)
    print("\n[Layer 0] Step 0.4: Transcript prevalence filter (diagnostic)")
    prevalence_masks = step04.run(adata_filtered, qc_log)
    plots.plot_step04(prevalence_masks, n_total=adata_filtered.n_vars)

    # 0.5 — Covariate completeness audit
    print("\n[Layer 0] Step 0.5: Covariate completeness audit")
    step05.run(adata_filtered, qc_log)

    # 0.6 — Pseudo-bulk aggregation (all transcripts, no prevalence subsetting)
    print("\n[Layer 0] Step 0.6: Pseudo-bulk aggregation")
    step06.run(adata_filtered, qc_log)
    plots.plot_step06(adata_filtered)

    # Write QC log
    print("\n[Layer 0] Writing QC log...")
    qc_log.write(OUT_QC)

    print("\n" + "=" * 60)
    print("Viscacha Layer 0 — COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()

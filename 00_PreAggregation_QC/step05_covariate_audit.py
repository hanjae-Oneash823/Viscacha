"""
Step 0.5 — Covariate completeness audit
Verifies all required DRIMSeq covariates are non-null for all donors.
No AnnData modification — audit only.
"""

import pandas as pd
import anndata as ad

from .utils.qc_log import QCLogger

REQUIRED_COVARIATES = ['age', 'sex', 'braak_stage', 'pct_counts_mt']


def run(adata_tx: ad.AnnData, qc_log: QCLogger) -> None:
    # Build donor-level covariate table
    # pct_counts_mt is per-barcode; we just check it's non-null across barcodes
    donor_cols = [c for c in REQUIRED_COVARIATES if c != 'pct_counts_mt']
    donor_meta = (
        adata_tx.obs[['donor'] + donor_cols]
        .drop_duplicates(subset='donor')
        .set_index('donor')
    )

    all_ok = True
    for cov in donor_cols:
        missing_donors = donor_meta[donor_meta[cov].isna()].index.tolist()
        if missing_donors:
            qc_log.log_covariate_missing(missing_donors, cov)
            all_ok = False
        else:
            print(f"  [OK] '{cov}' complete for all {len(donor_meta)} donors")

    # pct_counts_mt: check barcode-level
    n_mt_null = adata_tx.obs['pct_counts_mt'].isna().sum()
    if n_mt_null > 0:
        qc_log.flag("step_05", f"{n_mt_null} barcodes have null pct_counts_mt")
        all_ok = False
    else:
        n_donors = adata_tx.obs['donor'].nunique()
        print(f"  [OK] 'pct_counts_mt' non-null across all barcodes ({n_donors} donors)")

    qc_log.log_rin_absent()

    if all_ok:
        print("[Step 0.5] All covariates complete. No issues.")
    else:
        print("[Step 0.5] Covariate audit COMPLETE — see QC log for flags.")

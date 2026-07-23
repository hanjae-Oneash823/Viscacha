"""
Step 0.2 — Barcode merge and obs enrichment
Transfers doublet_score and pct_counts_mt from adata_sr into adata_tx,
then joins unified metadata onto adata_tx.obs via donor column.
"""

import anndata as ad
import pandas as pd

from layer0.config import SR_PATH, TX_PATH, TRANSFER_COLS
from layer0.utils.qc_log import QCLogger


def run(unified_meta: pd.DataFrame, qc_log: QCLogger) -> ad.AnnData:
    # Load adata_sr in backed mode — only access .obs, never .X
    print("[Step 0.2] Loading adata_sr (backed='r') for obs transfer...")
    adata_sr = ad.read_h5ad(SR_PATH, backed='r')
    sr_transfer = adata_sr.obs[TRANSFER_COLS].copy()
    adata_sr.file.close()
    print(f"  adata_sr obs extracted ({len(sr_transfer)} barcodes). File handle closed.")

    # Load adata_tx fully (74 MB)
    print("[Step 0.2] Loading adata_tx...")
    adata_tx = ad.read_h5ad(TX_PATH)
    print(f"  adata_tx loaded: {adata_tx.shape}")

    # Transfer doublet_score + pct_counts_mt via barcode index join
    n_before = adata_tx.n_obs
    adata_tx.obs = adata_tx.obs.join(sr_transfer, how='left')

    n_unmatched = adata_tx.obs[TRANSFER_COLS[0]].isna().sum()
    qc_log.log_unmatched_barcodes(int(n_unmatched), n_before)
    if n_unmatched > 0:
        qc_log.flag("step_02", f"{n_unmatched} barcodes had no match in adata_sr — these will be removed in Step 0.3 (null cell_type or doublet score)")

    # Join unified metadata onto obs via donor column using map()
    # (DataFrame.join(on=col) is unreliable with AnnData's obs; use map per column instead)
    for col in unified_meta.columns:
        adata_tx.obs[col] = adata_tx.obs['donor'].map(unified_meta[col].to_dict())

    # Report donors in AnnData not found in metadata
    adata_donors = set(adata_tx.obs['donor'].unique())
    meta_donors  = set(unified_meta.index)
    missing_in_meta = adata_donors - meta_donors
    if missing_in_meta:
        qc_log.flag("step_02", f"Donors in AnnData but not in metadata (will have NaN covariates): {sorted(missing_in_meta)}")

    extra_in_meta = meta_donors - adata_donors
    if extra_in_meta:
        print(f"  Metadata donors not in AnnData (expected — excluded samples): {sorted(extra_in_meta)}")

    print(f"[Step 0.2] obs enrichment complete. adata_tx obs columns: {list(adata_tx.obs.columns)}")
    return adata_tx

"""
Step 0.3 — Barcode confidence filter
Removes barcodes with doublet_score >= 0.3 or null cell_type / doublet_score.
Logs dropout per cell_type × donor and flags combinations losing >20%.
"""

import anndata as ad

from layer0.config import DOUBLET_THRESHOLD, CELL_TYPES, DROPOUT_FLAG_PCT
from layer0.utils.qc_log import QCLogger


def run(adata_tx: ad.AnnData, qc_log: QCLogger) -> ad.AnnData:
    obs = adata_tx.obs

    keep_mask = (
        obs['doublet_score'].notna() &
        obs['cell_type'].notna() &
        (obs['doublet_score'] < DOUBLET_THRESHOLD)
    )

    n_before_total = adata_tx.n_obs
    n_after_total  = keep_mask.sum()
    print(f"[Step 0.3] Total barcodes: {n_before_total} → {n_after_total} "
          f"({n_before_total - n_after_total} removed)")

    # Log per cell_type × donor
    for ct in obs['cell_type'].dropna().unique():
        for donor in obs['donor'].unique():
            ct_donor_mask   = (obs['cell_type'] == ct) & (obs['donor'] == donor)
            n_before = int(ct_donor_mask.sum())
            if n_before == 0:
                continue
            n_after = int((ct_donor_mask & keep_mask).sum())
            qc_log.log_barcode_drop(
                ct, donor, n_before, n_after,
                reason='doublet_score>=0.3 or null cell_type/doublet_score',
                dropout_flag_pct=DROPOUT_FLAG_PCT,
            )

    adata_filtered = adata_tx[keep_mask].copy()
    return adata_filtered

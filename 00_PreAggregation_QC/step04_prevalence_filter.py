"""
Step 0.4 — Transcript prevalence filter
Per cell type: keep transcripts detected (count > 0) in >= 40% of donors.
Returns a dict {cell_type: [transcript_names]} — does NOT modify the AnnData.
"""

import json
import numpy as np
import anndata as ad

from .config import MIN_PREVALENCE, OUT_META
from .utils.qc_log import QCLogger
from .utils.sparse_utils import donor_detection_vector


def run(adata_tx: ad.AnnData, qc_log: QCLogger) -> dict:
    prevalence_masks = {}
    n_total = adata_tx.n_vars

    for ct in adata_tx.obs['cell_type'].dropna().unique():
        sub  = adata_tx[adata_tx.obs['cell_type'] == ct]
        donors = sub.obs['donor'].unique()
        n_donors = len(donors)

        if n_donors == 0:
            continue

        # Accumulate detection counts across donors
        detection = np.zeros(n_total, dtype=np.int32)
        for donor in donors:
            d_mask = sub.obs['donor'] == donor
            rows   = np.where(d_mask)[0]
            chunk  = sub.X[rows, :]
            detection += donor_detection_vector(chunk).astype(np.int32)

        prevalence = detection / n_donors
        keep_idx   = np.where(prevalence >= MIN_PREVALENCE)[0]
        keep_names = adata_tx.var_names[keep_idx].tolist()

        prevalence_masks[ct] = keep_names
        qc_log.log_prevalence_filter(ct, n_total, len(keep_names))
        print(f"  [{ct}] {n_total} → {len(keep_names)} transcripts "
              f"({len(keep_names)/n_total*100:.1f}% kept, {n_donors} donors)")

    # Save mask to disk so it can be inspected or reloaded without re-running
    out_path = OUT_META / "prevalence_transcript_masks.json"
    with open(out_path, 'w') as f:
        json.dump(prevalence_masks, f)
    print(f"[Step 0.4] Prevalence masks saved to {out_path}")

    return prevalence_masks

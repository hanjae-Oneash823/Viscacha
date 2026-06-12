"""
Step 0.6 — Pseudo-bulk aggregation
For each cell type: sum transcript counts per donor, compute median_pct_mt.
Outputs two CSV pairs per cell type:
  counts_{stem}.csv        — donors × transcripts (AD + Control only, primary DRIMSeq input)
  metadata_{stem}.csv      — donor covariates
  counts_{stem}_sensitivity.csv   — Active control donors only
  metadata_{stem}_sensitivity.csv
"""

import numpy as np
import pandas as pd
import anndata as ad

from layer0.config import (
    CELL_TYPES, MIN_CELLS_PB, COND_ACTIVE, OUT_PB,
)
from layer0.utils.qc_log import QCLogger
from layer0.utils.sparse_utils import sparse_sum_rows


def _ct_stem(cell_type: str) -> str:
    return cell_type.replace(' ', '_')


def _aggregate_cell_type(
    adata: ad.AnnData,
    cell_type: str,
    prevalence_mask: list,
    qc_log: QCLogger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (count_df, meta_df) for all donors in this cell type.
    count_df: donors × filtered_transcripts (raw int counts)
    meta_df:  donors × covariates
    """
    sub = adata[adata.obs['cell_type'] == cell_type]

    # Build column index for prevalence-filtered transcripts (applied once)
    keep_set = set(prevalence_mask)
    keep_idx = np.array([i for i, v in enumerate(adata.var_names) if v in keep_set],
                        dtype=np.int64)
    keep_names = adata.var_names[keep_idx]

    pb_counts = {}
    pb_meta   = {}

    for donor in sub.obs['donor'].unique():
        d_mask  = sub.obs['donor'] == donor
        n_cells = int(d_mask.sum())

        if n_cells < MIN_CELLS_PB:
            qc_log.log_min_cell_exclusion(cell_type, donor, n_cells)
            continue

        rows  = np.where(d_mask)[0]
        chunk = sub.X[rows, :][:, keep_idx]   # slice rows then cols (adata_tx loaded fully)
        pb_counts[donor] = sparse_sum_rows(chunk)

        obs_slice = sub.obs[d_mask]
        pb_meta[donor] = {
            'condition':      obs_slice['condition'].iloc[0],
            'age':            obs_slice['age'].iloc[0],
            'sex':            obs_slice['sex'].iloc[0],
            'braak_stage':    obs_slice['braak_stage'].iloc[0],
            'median_pct_mt':  float(obs_slice['pct_counts_mt'].median()),
            'n_cells':        n_cells,
        }

    if not pb_counts:
        return pd.DataFrame(), pd.DataFrame()

    count_df = pd.DataFrame(pb_counts, index=keep_names).T   # donors × transcripts
    meta_df  = pd.DataFrame(pb_meta).T

    # Ensure counts are integer (sum of int32 sparse may be float64)
    count_df = count_df.astype(np.int64)

    return count_df, meta_df


def run(adata_tx: ad.AnnData, prevalence_masks: dict, qc_log: QCLogger) -> None:
    OUT_PB.mkdir(parents=True, exist_ok=True)

    for ct in CELL_TYPES:
        if ct not in prevalence_masks:
            qc_log.flag("step_06", f"No prevalence mask for cell type '{ct}' — skipping")
            continue

        mask = prevalence_masks[ct]
        count_df, meta_df = _aggregate_cell_type(adata_tx, ct, mask, qc_log)

        if count_df.empty:
            qc_log.flag("step_06", f"No donors met min_cells threshold for '{ct}'")
            continue

        stem = _ct_stem(ct)

        # Split primary (AD + Control) vs sensitivity (Active control)
        primary_donors    = meta_df[meta_df['condition'] != COND_ACTIVE].index
        sensitivity_donors = meta_df[meta_df['condition'] == COND_ACTIVE].index

        # Primary
        if len(primary_donors) > 0:
            count_df.loc[primary_donors].to_csv(OUT_PB / f"counts_{stem}.csv")
            meta_df.loc[primary_donors].to_csv(OUT_PB / f"metadata_{stem}.csv")

        # Sensitivity
        if len(sensitivity_donors) > 0:
            count_df.loc[sensitivity_donors].to_csv(OUT_PB / f"counts_{stem}_sensitivity.csv")
            meta_df.loc[sensitivity_donors].to_csv(OUT_PB / f"metadata_{stem}_sensitivity.csv")

        n_primary     = len(primary_donors)
        n_sensitivity = len(sensitivity_donors)
        n_tx = count_df.shape[1]
        print(f"  [{ct}] {n_primary} primary + {n_sensitivity} sensitivity donors, "
              f"{n_tx} transcripts → saved as '{stem}'")

    print(f"[Step 0.6] Pseudo-bulk outputs written to {OUT_PB}")

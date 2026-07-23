#!/usr/bin/env python3
"""
Build gene-level pseudobulk (donors x genes) from the SHORT-READ h5ad
(adata_sr.h5ad), mirroring layer0/step06_pseudobulk.py's donor-aggregation
logic exactly (same MIN_CELLS_PB, same primary/sensitivity donor split),
so its output is a drop-in gene-level counterpart to the long-read
transcript pseudobulk used by step15_gene_de.R.

adata_sr.h5ad already has donor/condition/cell_type/pct_counts_mt natively
(these are the source that gets transferred INTO the long-read object in
step02_barcode_merge.py) — it is just missing age/sex/braak_stage, which
are joined from outputs/layer0/metadata/unified_metadata.csv.

Output: outputs/layer0/pseudobulk_sr/{counts,metadata}_{cell_type}.csv
        (+ _sensitivity variants), same schema as outputs/layer0/pseudobulk/.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python build_sr_gene_pseudobulk.py
"""

import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

from layer0.config import SR_PATH, OUT_META, CELL_TYPES, MIN_CELLS_PB, COND_ACTIVE
from layer0.utils.sparse_utils import sparse_sum_rows

OUT_PB_SR = Path("/home/welcome3/Viscacha_pipeline/outputs/layer0/pseudobulk_sr")


def _ct_stem(cell_type: str) -> str:
    return cell_type.replace(" ", "_")


def _aggregate_cell_type(adata: ad.AnnData, cell_type: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    sub = adata[adata.obs["cell_type"] == cell_type]
    keep_names = adata.var_names

    pb_counts = {}
    pb_meta = {}

    for donor in sub.obs["donor"].unique():
        d_mask = (sub.obs["donor"] == donor).values
        n_cells = int(d_mask.sum())
        if n_cells < MIN_CELLS_PB:
            continue

        rows = np.where(d_mask)[0]
        chunk = sub.X[rows, :]
        pb_counts[donor] = sparse_sum_rows(chunk)

        obs_slice = sub.obs[d_mask]
        pb_meta[donor] = {
            "condition": obs_slice["condition"].iloc[0],
            "age": obs_slice["age"].iloc[0],
            "sex": obs_slice["sex"].iloc[0],
            "braak_stage": obs_slice["braak_stage"].iloc[0],
            "median_pct_mt": float(obs_slice["pct_counts_mt"].median()),
            "n_cells": n_cells,
        }

    if not pb_counts:
        return pd.DataFrame(), pd.DataFrame()

    count_df = pd.DataFrame(pb_counts, index=keep_names).T
    meta_df = pd.DataFrame(pb_meta).T
    count_df = count_df.astype(np.int64)
    return count_df, meta_df


def main() -> None:
    OUT_PB_SR.mkdir(parents=True, exist_ok=True)

    print("Loading unified_metadata.csv ...")
    unified_meta = pd.read_csv(OUT_META / "unified_metadata.csv", index_col="donor_id")

    print(f"Loading adata_sr (backed='r') from {SR_PATH} ...")
    adata = ad.read_h5ad(SR_PATH, backed="r")
    print(f"  {adata.shape[0]} cells x {adata.shape[1]} genes")

    # adata_sr already has donor/condition/cell_type/pct_counts_mt.
    # Only age/sex/braak_stage need to be joined in, same as step02_barcode_merge.py.
    for col in ["age", "sex", "braak_stage"]:
        adata.obs[col] = adata.obs["donor"].map(unified_meta[col].to_dict())

    for ct in CELL_TYPES:
        count_df, meta_df = _aggregate_cell_type(adata, ct)
        if count_df.empty:
            print(f"  [{ct}] no donors met MIN_CELLS_PB — skipping")
            continue

        stem = _ct_stem(ct)
        primary_donors = meta_df[meta_df["condition"] != COND_ACTIVE].index
        sensitivity_donors = meta_df[meta_df["condition"] == COND_ACTIVE].index

        if len(primary_donors) > 0:
            count_df.loc[primary_donors].to_csv(OUT_PB_SR / f"counts_{stem}.csv")
            meta_df.loc[primary_donors].to_csv(OUT_PB_SR / f"metadata_{stem}.csv")
        if len(sensitivity_donors) > 0:
            count_df.loc[sensitivity_donors].to_csv(OUT_PB_SR / f"counts_{stem}_sensitivity.csv")
            meta_df.loc[sensitivity_donors].to_csv(OUT_PB_SR / f"metadata_{stem}_sensitivity.csv")

        print(f"  [{ct}] {len(primary_donors)} primary + {len(sensitivity_donors)} sensitivity donors, "
              f"{count_df.shape[1]} genes -> saved as '{stem}'")

    print(f"\nDone. SR gene pseudobulk written to {OUT_PB_SR}")


if __name__ == "__main__":
    main()

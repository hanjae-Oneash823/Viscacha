#!/usr/bin/env python3
"""
Step 15 (Wilcoxon variant): Gene-level differential expression via
scanpy.tl.rank_genes_groups(method="wilcoxon"), run on the SAME donor-level
pseudobulk matrices as step15_gene_de_sr.R / step15_gene_de_gx.R (DESeq2) --
a nonparametric complement/sanity-check, not a per-cell test. A per-cell
Wilcoxon on raw single cells would be pseudoreplicated (cells from the same
donor are not independent draws); donors remain the unit of replication here.

Mirrors step15_gene_de_sr.R / step15_gene_de_gx.R's per-cell-type loop,
MIN_GENE_TOTAL=10 gene filter, and NA-covariate donor dropping. Wilcoxon has
no covariate-adjustment mechanism (unlike DESeq2's GLM), so age/sex/
median_pct_mt are NOT adjusted for here -- a trade-off, not an oversight.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python 01_ViscachaDTU_Analysis/step15_gene_de_wilcoxon.py
     (from /home/welcome3/Viscacha_pipeline)
Output: outputs/01_ViscachaDTU_Analysis/gene_level_de_results_wilcoxon_SR.csv
        outputs/01_ViscachaDTU_Analysis/gene_level_de_results_wilcoxon_GX.csv
  columns: gene_id, cell_type, log2FC, pval, padj, score
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path

REPO_ROOT = Path("/home/welcome3/Viscacha_pipeline")
OUT_DIR = REPO_ROOT / "outputs/01_ViscachaDTU_Analysis"
PB_ROOT = REPO_ROOT / "outputs/00_PreAggregation_QC"

# Same cell types as step15_gene_de_sr.R / _gx.R (Lymphocyte dropped: too few donors)
CELL_TYPES = [
    "Excitatory_neuron",
    "Inhibitory_neuron",
    "Oligodendrocyte",
    "OPC",
    "Astrocyte",
    "Microglia",
    "Vascular_cell",
]
MIN_GENE_TOTAL = 10
COND_LEVELS = ["Control", "AD"]
MODEL_VARS = ["condition", "age", "sex", "median_pct_mt"]


def run_one(in_dir: Path, cell_type: str) -> pd.DataFrame | None:
    f_counts = in_dir / f"counts_{cell_type}.csv"
    f_meta = in_dir / f"metadata_{cell_type}.csv"
    if not f_counts.exists() or not f_meta.exists():
        print(f"  [{cell_type}] pseudobulk not found -- skipping")
        return None

    counts = pd.read_csv(f_counts, index_col=0)  # donors x genes
    meta = pd.read_csv(f_meta, index_col=0)  # donors x covariates

    # Light pre-filter on genes (all donors), same as R's rowSums(gene_mat) >= MIN_GENE_TOTAL
    gene_totals = counts.sum(axis=0)
    counts = counts.loc[:, gene_totals >= MIN_GENE_TOTAL]

    # Align metadata; drop donors with NA in any model covariate
    meta = meta.loc[counts.index]
    keep = meta[MODEL_VARS].notna().all(axis=1)
    n_dropped = int((~keep).sum())
    if n_dropped:
        print(f"  [{cell_type}] dropping {n_dropped} donor(s) with NA covariates")
    counts = counts.loc[keep]
    meta = meta.loc[keep].copy()

    meta["condition"] = pd.Categorical(meta["condition"], categories=COND_LEVELS)
    if meta["condition"].dropna().nunique() < 2 or counts.shape[0] < 4:
        print(f"  [{cell_type}] insufficient samples/conditions -- skipping")
        return None

    adata = ad.AnnData(
        X=counts.values.astype(np.float64),
        obs=meta,
        var=pd.DataFrame(index=counts.columns),
    )
    adata.var_names_make_unique()

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby="condition",
        groups=["AD"],
        reference="Control",
        method="wilcoxon",
        corr_method="benjamini-hochberg",
    )

    res = sc.get.rank_genes_groups_df(adata, group="AD")
    res = res.rename(
        columns={
            "names": "gene_id",
            "logfoldchanges": "log2FC",
            "pvals": "pval",
            "pvals_adj": "padj",
            "scores": "score",
        }
    )
    res["cell_type"] = cell_type
    n_sig = int((res["padj"] < 0.05).sum())
    print(f"  [{cell_type}] {len(res)} genes tested, {n_sig} DE at padj<0.05")
    return res[["gene_id", "cell_type", "log2FC", "pval", "padj", "score"]]


def run_source(in_dir: Path, label: str, out_name: str) -> None:
    print("=" * 60)
    print(f"Step 15 (Wilcoxon) -- Gene-level DE on {label} pseudobulk")
    print("=" * 60)
    all_de = []
    for ct in CELL_TYPES:
        print(f"\n{ct}")
        de = run_one(in_dir, ct)
        if de is not None:
            all_de.append(de)
    if not all_de:
        raise RuntimeError(f"No gene-level DE results produced for {label}.")
    combined = pd.concat(all_de, ignore_index=True)
    out_path = OUT_DIR / out_name
    combined.to_csv(out_path, index=False)
    print(f"\nSaved {len(combined)} gene x cell_type rows to {out_path}")


if __name__ == "__main__":
    run_source(PB_ROOT / "pseudobulk_sr", "SR (short-read)", "gene_level_de_results_wilcoxon_SR.csv")
    run_source(PB_ROOT / "pseudobulk_gx", "GX (long-read gene-level)", "gene_level_de_results_wilcoxon_GX.csv")

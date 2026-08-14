#!/usr/bin/env python3
"""
Method-comparison (not platform-comparison): for each of SR and GX gene-level
pseudobulk separately, do DESeq2 (parametric GLM, existing step15_gene_de_sr.R
/ _gx.R) and scanpy rank_genes_groups Wilcoxon (step15_gene_de_wilcoxon.py) --
run on the SAME donor x gene counts -- agree on which genes look different
between AD and Control?

padj can't be compared directly: Wilcoxon's discreteness floor (~1.7e-4 to
2.2e-4 at n=8-13 donors/group, see step15_gene_de_wilcoxon.py run log) means
it never survives BH correction genome-wide, regardless of DESeq2's calls.
So this compares RAW p-value rank and log2FC sign/magnitude instead -- do the
two methods point at the same genes even though only one of them can call
significance at this sample size.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python 01_ViscachaDTU_Analysis/compare_deg_wilcoxon_vs_deseq2.py
     (from /home/welcome3/Viscacha_pipeline)
Output: outputs/01_ViscachaDTU_Analysis/wilcoxon_vs_deseq2_{SR,GX}_merged.csv
        printed per-cell-type summary (Spearman rho on log2FC and -log10 pval,
        log2FC sign concordance, top-50-by-pval Jaccard overlap)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

OUT_DIR = Path("/home/welcome3/Viscacha_pipeline/outputs/01_ViscachaDTU_Analysis")
TOP_N = 50


def compare_source(label: str, deseq2_file: str, wilcoxon_file: str, out_name: str) -> None:
    deseq2 = pd.read_csv(OUT_DIR / deseq2_file)
    wilcoxon = pd.read_csv(OUT_DIR / wilcoxon_file)

    merged = deseq2.merge(
        wilcoxon, on=["gene_id", "cell_type"], suffixes=("_deseq2", "_wilcoxon")
    )
    merged.to_csv(OUT_DIR / out_name, index=False)

    print("=" * 70)
    print(f"{label}: DESeq2 vs Wilcoxon (same donor x gene counts, n={len(merged)} rows)")
    print("=" * 70)

    rows = []
    for ct, sub in merged.groupby("cell_type"):
        sub = sub.dropna(subset=["pval_deseq2", "pval_wilcoxon", "log2FC_deseq2", "log2FC_wilcoxon"])
        if len(sub) < 10:
            continue

        rho_lfc, _ = spearmanr(sub["log2FC_deseq2"], sub["log2FC_wilcoxon"])
        rho_p, _ = spearmanr(-np.log10(sub["pval_deseq2"]), -np.log10(sub["pval_wilcoxon"]))
        sign_agree = np.mean(np.sign(sub["log2FC_deseq2"]) == np.sign(sub["log2FC_wilcoxon"]))

        top_deseq2 = set(sub.nsmallest(TOP_N, "pval_deseq2")["gene_id"])
        top_wilcoxon = set(sub.nsmallest(TOP_N, "pval_wilcoxon")["gene_id"])
        jaccard = len(top_deseq2 & top_wilcoxon) / len(top_deseq2 | top_wilcoxon)

        rows.append(
            {
                "cell_type": ct,
                "n_genes": len(sub),
                "spearman_rho_log2FC": round(rho_lfc, 3),
                "spearman_rho_neglog10p": round(rho_p, 3),
                "log2FC_sign_agree_frac": round(sign_agree, 3),
                f"top{TOP_N}_pval_jaccard": round(jaccard, 3),
            }
        )

    summary = pd.DataFrame(rows)
    print(summary.to_string(index=False))
    print()


if __name__ == "__main__":
    compare_source(
        "SR (short-read)",
        "gene_level_de_results_SR.csv",
        "gene_level_de_results_wilcoxon_SR.csv",
        "wilcoxon_vs_deseq2_SR_merged.csv",
    )
    compare_source(
        "GX (long-read gene-level)",
        "gene_level_de_results_GX.csv",
        "gene_level_de_results_wilcoxon_GX.csv",
        "wilcoxon_vs_deseq2_GX_merged.csv",
    )

"""
Diagnostic: per-gene pseudobulk concordance between adata_sr (short-read,
gene-level) and the two long-read views -- adata_transcript_loose (transcript-
level, collapsed to gene via ENSG_ID) and adata_gene_loose (gene-level).

Unlike explore_sr_tx_depth.py / explore_sr_gx_depth.py (which compare total
counts PER CELL), this compares total counts PER GENE, summed across cells
(pseudobulk). Genes are matched across files via Ensembl ID (ENSG_ID), with
gene-symbol fallback for SR var_names that aren't already an ENSG ID.

TX and GX barcodes are both proper subsets of SR's barcodes (verified: 0
orphans). To keep each pairwise comparison apples-to-apples, the SR side of
each comparison is summed only over the cells that are also present in the
other file -- e.g. the "SR vs TX" pseudobulk uses SR counts from
SR-cells-that-are-also-in-TX, not all of SR. Both matched sums are
accumulated in a single chunked pass over the 37 GB SR matrix (backed mode,
row-chunked, boolean-masked per chunk) so SR.X is only read once.

Standalone exploratory script, not part of the run_layer0 pipeline.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python -m 00_PreAggregation_QC.explore_gene_pseudobulk_concordance
     (from /home/welcome3/Viscacha_pipeline)
Output: outputs/00_PreAggregation_QC/plots/gene_pseudobulk/*.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

from .config import OUT_DIR, SR_PATH, TX_PATH, GX_PATH

OUT_PLOTS = OUT_DIR / "plots" / "gene_pseudobulk"
SR_CHUNK_SIZE = 20000


def _savefig(fig: plt.Figure, name: str) -> Path:
    OUT_PLOTS.mkdir(parents=True, exist_ok=True)
    path = OUT_PLOTS / name
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] saved {path.name}")
    return path


# ---------------------------------------------------------------------------
# Data loading / pseudobulk construction
# ---------------------------------------------------------------------------
def load_tx_gx():
    print("Loading adata_transcript_loose (full, for X sum)...")
    tx = ad.read_h5ad(TX_PATH)
    print("Loading adata_gene_loose (full, for X sum)...")
    gx = ad.read_h5ad(GX_PATH)
    return tx, gx


def tx_gene_pseudobulk(tx: ad.AnnData) -> pd.Series:
    """Collapse transcript-level pseudobulk to gene-level via ENSG_ID."""
    transcript_total = np.asarray(tx.X.sum(axis=0)).ravel()
    s = pd.Series(transcript_total, index=tx.var["ENSG_ID"].values)
    return s.groupby(level=0).sum()


def gx_gene_pseudobulk(gx: ad.AnnData) -> pd.Series:
    gene_total = np.asarray(gx.X.sum(axis=0)).ravel()
    s = pd.Series(gene_total, index=gx.var["ENSG_ID"].values)
    return s.groupby(level=0).sum()


def build_sr_to_ensg_map(sr_var_names: pd.Index, tx: ad.AnnData, gx: ad.AnnData) -> dict:
    """SR var_names are a mix of gene symbols and raw ENSG IDs (used when a
    symbol is duplicated/missing). Map each to a canonical ENSG_ID using the
    TX and GX panels, which both carry an explicit ENSG_ID column."""
    symbol_to_ensg = {}
    symbol_to_ensg.update(dict(zip(tx.var["gene_name"], tx.var["ENSG_ID"])))
    symbol_to_ensg.update(dict(zip(gx.var_names, gx.var["ENSG_ID"])))
    known_ensg = set(tx.var["ENSG_ID"]) | set(gx.var["ENSG_ID"])

    mapping = {}
    for v in sr_var_names:
        if v.startswith("ENSG") and v in known_ensg:
            mapping[v] = v
        elif v in symbol_to_ensg:
            mapping[v] = symbol_to_ensg[v]
    return mapping


def sr_pseudobulk_matched(tx_barcodes: set, gx_barcodes: set) -> tuple[pd.Series, pd.Series]:
    """Single chunked pass over the 37 GB SR matrix (backed mode). Accumulates
    two per-gene sums: SR counts from cells also present in TX, and SR counts
    from cells also present in GX."""
    print("Loading adata_sr (backed, chunked X sum)...")
    sr = ad.read_h5ad(SR_PATH, backed="r")
    n_obs, n_genes = sr.shape
    sum_tx = np.zeros(n_genes)
    sum_gx = np.zeros(n_genes)
    obs_names = sr.obs_names

    for start in range(0, n_obs, SR_CHUNK_SIZE):
        end = min(start + SR_CHUNK_SIZE, n_obs)
        bc_chunk = obs_names[start:end]
        mask_tx = bc_chunk.isin(tx_barcodes)
        mask_gx = bc_chunk.isin(gx_barcodes)
        if not (mask_tx.any() or mask_gx.any()):
            continue
        chunk_X = sr.X[start:end]
        if mask_tx.any():
            sum_tx += np.asarray(chunk_X[mask_tx].sum(axis=0)).ravel()
        if mask_gx.any():
            sum_gx += np.asarray(chunk_X[mask_gx].sum(axis=0)).ravel()
        print(f"  ...processed {end:,}/{n_obs:,} SR cells", end="\r")
    print()
    sr.file.close()

    sr_var_names = sr.var_names
    return (pd.Series(sum_tx, index=sr_var_names), pd.Series(sum_gx, index=sr_var_names))


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_scatter(merged: pd.DataFrame, other_label: str, suffix: str) -> None:
    sub = merged[(merged["sr"] > 0) & (merged["other"] > 0)]
    r_raw, _ = pearsonr(sub["sr"], sub["other"])
    r_log, _ = pearsonr(np.log1p(sub["sr"]), np.log1p(sub["other"]))
    rho, _ = spearmanr(sub["sr"], sub["other"])

    fig, ax = plt.subplots(figsize=(7, 6))
    hb = ax.hexbin(sub["sr"], sub["other"], gridsize=50, bins="log", cmap="turbo",
                    xscale="log", yscale="log", mincnt=1)
    fig.colorbar(hb, ax=ax, label="genes (log count)")
    ax.set_xlabel("SR per-gene pseudobulk total counts (log scale)")
    ax.set_ylabel(f"{other_label} per-gene pseudobulk total counts (log scale)")
    ax.set_title(f"Per-gene pseudobulk: SR vs {other_label}", weight="semibold", pad=12)
    ax.text(0.03, 0.97,
            f"Pearson r (raw) = {r_raw:.3f}\nPearson r (log-log) = {r_log:.3f}\nSpearman ρ = {rho:.3f}\nn = {len(sub):,} genes",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#cccccc"))
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, f"01_sr_vs_{suffix}_gene_pseudobulk_scatter.png")


def plot_marginal_histograms(merged_tx: pd.DataFrame, merged_gx: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    panels = [
        (axes[0, 0], merged_tx["sr"], "SR total counts (matched to TX cells)", "#3A86FF"),
        (axes[0, 1], merged_tx["other"], "TX-collapsed-to-gene total counts", "#FF006E"),
        (axes[1, 0], merged_gx["sr"], "SR total counts (matched to GX cells)", "#3A86FF"),
        (axes[1, 1], merged_gx["other"], "GX total counts", "#FF006E"),
    ]
    for ax, vals, label, color in panels:
        vals = vals[vals > 0]
        bins = np.logspace(np.log10(vals.min()), np.log10(vals.max()), 50)
        ax.hist(vals, bins=bins, color=color, edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.axvline(vals.median(), color="black", linestyle="--", linewidth=1,
                   label=f"median = {vals.median():.0f}")
        ax.set_xlabel(f"{label} (log scale)")
        ax.set_ylabel("# genes")
        ax.set_title(label, weight="semibold", fontsize=10)
        ax.legend(frameon=False, fontsize=8)
        sns.despine(ax=ax)
    fig.suptitle("Per-gene pseudobulk depth distributions", weight="semibold", y=1.0)
    fig.tight_layout()
    _savefig(fig, "02_marginal_gene_pseudobulk_histograms.png")


def plot_ratio_distributions(merged_tx: pd.DataFrame, merged_gx: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, merged, label, color in [
        (axes[0], merged_tx, "TX", "#8338EC"),
        (axes[1], merged_gx, "GX", "#FB5607"),
    ]:
        sub = merged[(merged["sr"] > 0) & (merged["other"] > 0)]
        ratio = sub["sr"] / sub["other"]
        bins = np.logspace(np.log10(ratio.min()), np.log10(ratio.max()), 60)
        ax.hist(ratio, bins=bins, color=color, edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.axvline(ratio.median(), color="black", linestyle="--", linewidth=1,
                   label=f"median = {ratio.median():.1f}x")
        ax.set_xlabel(f"SR / {label} pseudobulk ratio per gene (log scale)")
        ax.set_ylabel("# genes")
        ax.set_title(f"SR vs {label}: fold-difference per gene", weight="semibold", fontsize=10)
        ax.legend(frameon=False)
        sns.despine(ax=ax)
    fig.suptitle("Per-gene pseudobulk ratio distributions", weight="semibold", y=1.02)
    fig.tight_layout()
    _savefig(fig, "03_gene_pseudobulk_ratio_distributions.png")


def plot_efficiency_vs_expression(merged_tx: pd.DataFrame, merged_gx: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, merged, label, color in [
        (axes[0], merged_tx, "TX", "#06D6A0"),
        (axes[1], merged_gx, "GX", "#FFB703"),
    ]:
        df = merged[(merged["sr"] > 0) & (merged["other"] > 0)].copy()
        df["efficiency"] = df["other"] / df["sr"]
        df["sr_bin"] = pd.qcut(df["sr"], q=20, duplicates="drop")
        binned = df.groupby("sr_bin", observed=True).agg(
            sr_mid=("sr", "median"),
            eff_median=("efficiency", "median"),
            eff_q25=("efficiency", lambda x: x.quantile(0.25)),
            eff_q75=("efficiency", lambda x: x.quantile(0.75)),
        ).reset_index(drop=True)

        ax.plot(binned["sr_mid"], binned["eff_median"], color=color, linewidth=2.5, marker="o", markersize=4)
        ax.fill_between(binned["sr_mid"], binned["eff_q25"], binned["eff_q75"], color=color, alpha=0.25, label="IQR")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("SR per-gene pseudobulk (log scale, binned)")
        ax.set_ylabel(f"{label}/SR capture efficiency (log scale)")
        ax.set_title(f"Does {label} capture depend on gene expression level?", weight="semibold", fontsize=10)
        ax.legend(frameon=False)
        sns.despine(ax=ax)
    fig.suptitle("Capture efficiency vs. gene expression level", weight="semibold", y=1.02)
    fig.tight_layout()
    _savefig(fig, "04_capture_efficiency_vs_gene_expression.png")


def plot_panel_missingness(sr_matched: pd.Series, in_panel: pd.Series, label: str, suffix: str) -> None:
    """Are genes absent from the TX/GX panel simply the lowly-expressed ones in SR?"""
    df = pd.DataFrame({"sr_total": sr_matched, "in_panel": in_panel})
    df = df[df["sr_total"] > 0]

    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.logspace(np.log10(df["sr_total"].min()), np.log10(df["sr_total"].max()), 50)
    for status, color in [(True, "#06D6A0"), (False, "#FB5607")]:
        sub = df[df["in_panel"] == status]["sr_total"]
        lbl = f"In {label} panel" if status else f"Not in {label} panel"
        ax.hist(sub, bins=bins, alpha=0.55, color=color, label=f"{lbl} (n={len(sub):,})",
                density=True, edgecolor="white", linewidth=0.2)
    ax.set_xscale("log")
    ax.set_xlabel("SR per-gene pseudobulk total counts (log scale)")
    ax.set_ylabel("density")
    ax.set_title(f"Are genes missing from {label} simply the lowly-expressed ones in SR?",
                 weight="semibold", pad=12, fontsize=11)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, f"05_gene_missingness_by_sr_expression_{suffix}.png")


def main():
    tx, gx = load_tx_gx()
    tx_bc = set(tx.obs_names)
    gx_bc = set(gx.obs_names)

    tx_gene_pb = tx_gene_pseudobulk(tx)
    gx_gene_pb = gx_gene_pseudobulk(gx)
    tx_ensg_set = set(tx.var["ENSG_ID"])
    gx_ensg_set = set(gx.var["ENSG_ID"])

    sr = ad.read_h5ad(SR_PATH, backed="r")
    sr_var_names = sr.var_names
    sr_to_ensg = build_sr_to_ensg_map(sr_var_names, tx, gx)
    sr.file.close()

    sr_matched_tx, sr_matched_gx = sr_pseudobulk_matched(tx_bc, gx_bc)

    in_tx_panel = pd.Series({v: (sr_to_ensg.get(v) in tx_ensg_set) for v in sr_var_names})
    in_gx_panel = pd.Series({v: (sr_to_ensg.get(v) in gx_ensg_set) for v in sr_var_names})

    print(f"\nn SR genes: {len(sr_var_names):,}")
    print(f"n SR genes mapped to an ENSG_ID: {len(sr_to_ensg):,}")
    print(f"n SR genes present in TX panel: {in_tx_panel.sum():,}")
    print(f"n SR genes present in GX panel: {in_gx_panel.sum():,}\n")

    # Collapse SR matched-pseudobulk from var_name -> ENSG_ID (sum any collisions)
    sr_matched_tx_ensg = sr_matched_tx.rename(index=sr_to_ensg).dropna()
    sr_matched_tx_ensg = sr_matched_tx_ensg[sr_matched_tx_ensg.index.notna()]
    sr_matched_tx_ensg = pd.Series(sr_matched_tx_ensg.values, index=sr_matched_tx_ensg.index).groupby(level=0).sum()

    sr_matched_gx_ensg = sr_matched_gx.rename(index=sr_to_ensg).dropna()
    sr_matched_gx_ensg = pd.Series(sr_matched_gx_ensg.values, index=sr_matched_gx_ensg.index).groupby(level=0).sum()

    merged_tx = pd.concat(
        [sr_matched_tx_ensg.rename("sr"), tx_gene_pb.rename("other")], axis=1, join="inner"
    )
    merged_gx = pd.concat(
        [sr_matched_gx_ensg.rename("sr"), gx_gene_pb.rename("other")], axis=1, join="inner"
    )
    print(f"n genes matched SR<->TX (via ENSG_ID): {len(merged_tx):,}")
    print(f"n genes matched SR<->GX (via ENSG_ID): {len(merged_gx):,}\n")

    plot_scatter(merged_tx, "TX (transcript-collapsed)", "tx")
    plot_scatter(merged_gx, "GX (gene-level)", "gx")
    plot_marginal_histograms(merged_tx, merged_gx)
    plot_ratio_distributions(merged_tx, merged_gx)
    plot_efficiency_vs_expression(merged_tx, merged_gx)
    plot_panel_missingness(sr_matched_tx, in_tx_panel, "TX", "tx")
    plot_panel_missingness(sr_matched_gx, in_gx_panel, "GX", "gx")

    print(f"\nAll plots saved to: {OUT_PLOTS}")


if __name__ == "__main__":
    main()

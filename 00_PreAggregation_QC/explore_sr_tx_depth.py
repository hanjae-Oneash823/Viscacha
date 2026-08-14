"""
Diagnostic: count-depth relationship between adata_sr (short-read, gene-level)
and adata_transcript_loose_filtering_for_bulk_analysis (long-read, transcript-
level). Both raw inputs already carry a precomputed `total_counts` in .obs, so
this never touches the (37 GB) adata_sr.X matrix -- only .obs is read, in
backed mode.

Standalone exploratory script, not part of the run_layer0 pipeline and does
not write anything the pipeline consumes downstream.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python -m 00_PreAggregation_QC.explore_sr_tx_depth
     (from /home/welcome3/Viscacha_pipeline)
Output: outputs/00_PreAggregation_QC/plots/sr_tx_depth/*.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import anndata as ad
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

from .config import OUT_DIR, SR_PATH, TX_PATH

OUT_PLOTS_DEPTH = OUT_DIR / "plots" / "sr_tx_depth"


def _savefig(fig: plt.Figure, name: str) -> Path:
    OUT_PLOTS_DEPTH.mkdir(parents=True, exist_ok=True)
    path = OUT_PLOTS_DEPTH / name
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] saved {path.name}")
    return path


def load_obs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Backed-mode obs-only load of both raw AnnData files. Never touches .X."""
    print("Loading adata_sr (backed, obs only)...")
    sr = ad.read_h5ad(SR_PATH, backed="r")
    sr_obs = sr.obs[["donor", "condition", "cell_type", "library_kit", "total_counts",
                      "n_genes_by_counts", "pct_counts_mt"]].copy()
    sr.file.close()

    print("Loading adata_transcript_loose (backed, obs only)...")
    tx = ad.read_h5ad(TX_PATH, backed="r")
    tx_obs = tx.obs[["donor", "condition", "cell_type", "library_kit", "total_counts",
                      "n_genes_by_counts"]].copy()
    tx.file.close()

    return sr_obs, tx_obs


def build_merged(sr_obs: pd.DataFrame, tx_obs: pd.DataFrame) -> pd.DataFrame:
    merged = sr_obs.join(tx_obs, lsuffix="_sr", rsuffix="_tx", how="inner")
    merged["ratio_sr_to_tx"] = (
        merged["total_counts_sr"] / merged["total_counts_tx"].replace(0, np.nan)
    )
    return merged


# ---------------------------------------------------------------------------
# 1. SR vs TX depth — hexbin scatter, log-log, with correlation annotated
# ---------------------------------------------------------------------------
def plot_depth_scatter(merged: pd.DataFrame) -> None:
    r_raw, _   = pearsonr(merged["total_counts_sr"], merged["total_counts_tx"])
    r_log, _   = pearsonr(np.log1p(merged["total_counts_sr"]),
                           np.log1p(merged["total_counts_tx"]))
    rho, _     = spearmanr(merged["total_counts_sr"], merged["total_counts_tx"])

    fig, ax = plt.subplots(figsize=(7, 6))
    hb = ax.hexbin(merged["total_counts_sr"], merged["total_counts_tx"],
                    gridsize=60, bins="log", cmap="turbo",
                    xscale="log", yscale="log", mincnt=1)
    fig.colorbar(hb, ax=ax, label="cells (log count)")
    ax.set_xlabel("SR total_counts (log scale)")
    ax.set_ylabel("TX total_counts (log scale)")
    ax.set_title("Per-cell depth: short-read vs long-read transcript", weight="semibold", pad=12)
    ax.text(0.03, 0.97,
            f"Pearson r (raw) = {r_raw:.3f}\nPearson r (log-log) = {r_log:.3f}\nSpearman ρ = {rho:.3f}\nn = {len(merged):,} cells",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#cccccc"))
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "01_sr_vs_tx_depth_scatter.png")


# ---------------------------------------------------------------------------
# 2. Marginal histograms of SR depth and TX depth
# ---------------------------------------------------------------------------
def plot_marginal_histograms(merged: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, col, label, color in [
        (axes[0], "total_counts_sr", "SR total_counts", "#3A86FF"),
        (axes[1], "total_counts_tx", "TX total_counts", "#FF006E"),
    ]:
        vals = merged[col]
        bins = np.logspace(np.log10(max(vals.min(), 1)), np.log10(vals.max()), 50)
        ax.hist(vals, bins=bins, color=color, edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.axvline(vals.median(), color="black", linestyle="--", linewidth=1,
                   label=f"median = {vals.median():.0f}")
        ax.set_xlabel(f"{label} (log scale)")
        ax.set_ylabel("# cells")
        ax.set_title(label, weight="semibold")
        ax.legend(frameon=False, fontsize=8)
        sns.despine(ax=ax)
    fig.suptitle("Depth distributions", weight="semibold", y=1.02)
    fig.tight_layout()
    _savefig(fig, "02_marginal_depth_histograms.png")


# ---------------------------------------------------------------------------
# 3. SR/TX ratio distribution
# ---------------------------------------------------------------------------
def plot_ratio_distribution(merged: pd.DataFrame) -> None:
    ratio = merged["ratio_sr_to_tx"].dropna()
    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.logspace(np.log10(max(ratio.min(), 0.1)), np.log10(ratio.max()), 60)
    ax.hist(ratio, bins=bins, color="#8338EC", edgecolor="white", linewidth=0.3)
    ax.set_xscale("log")
    ax.axvline(ratio.median(), color="black", linestyle="--", linewidth=1,
               label=f"median = {ratio.median():.1f}x")
    ax.set_xlabel("SR / TX depth ratio per cell (log scale)")
    ax.set_ylabel("# cells")
    ax.set_title("How many fold deeper is SR than TX, per cell?", weight="semibold", pad=12)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "03_sr_tx_ratio_distribution.png")


# ---------------------------------------------------------------------------
# 4. Capture efficiency (TX/SR) vs SR depth, binned
# ---------------------------------------------------------------------------
def plot_efficiency_vs_depth(merged: pd.DataFrame) -> None:
    df = merged.copy()
    df["efficiency"] = df["total_counts_tx"] / df["total_counts_sr"]
    df["sr_bin"] = pd.qcut(df["total_counts_sr"], q=25, duplicates="drop")
    binned = df.groupby("sr_bin", observed=True).agg(
        sr_mid=("total_counts_sr", "median"),
        eff_median=("efficiency", "median"),
        eff_q25=("efficiency", lambda x: x.quantile(0.25)),
        eff_q75=("efficiency", lambda x: x.quantile(0.75)),
    ).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.plot(binned["sr_mid"], binned["eff_median"], color="#06D6A0", linewidth=2.5, marker="o", markersize=4)
    ax.fill_between(binned["sr_mid"], binned["eff_q25"], binned["eff_q75"],
                     color="#06D6A0", alpha=0.25, label="IQR")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("SR total_counts (log scale, binned)")
    ax.set_ylabel("TX/SR capture efficiency (log scale)")
    ax.set_title("Does long-read capture efficiency change with depth?", weight="semibold", pad=12)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "04_capture_efficiency_vs_depth.png")


# ---------------------------------------------------------------------------
# 5. Per-cell-type scatter facets
# ---------------------------------------------------------------------------
def plot_by_cell_type(merged: pd.DataFrame) -> None:
    cell_types = sorted(merged["cell_type_tx"].dropna().unique())
    n = len(cell_types)
    n_col = min(4, n)
    n_row = int(np.ceil(n / n_col))
    fig, axes = plt.subplots(n_row, n_col, figsize=(n_col * 3.3, n_row * 3.3), squeeze=False)

    for i, ct in enumerate(cell_types):
        ax = axes[i // n_col][i % n_col]
        sub = merged[merged["cell_type_tx"] == ct]
        r_log, _ = pearsonr(np.log1p(sub["total_counts_sr"]), np.log1p(sub["total_counts_tx"]))
        gridsize = max(10, min(40, int(np.sqrt(len(sub)))))
        ax.hexbin(sub["total_counts_sr"], sub["total_counts_tx"],
                  gridsize=gridsize, bins="log", cmap="rainbow",
                  xscale="log", yscale="log", mincnt=1)
        ax.set_title(f"{ct}\nn={len(sub):,}, r(log)={r_log:.2f}", fontsize=9)
        sns.despine(ax=ax)

    for j in range(n, n_row * n_col):
        axes[j // n_col][j % n_col].axis("off")

    fig.suptitle("SR vs TX depth, by cell type", weight="semibold", y=1.02)
    fig.tight_layout()
    _savefig(fig, "05_depth_scatter_by_celltype.png")


# ---------------------------------------------------------------------------
# 6. Per-donor scatter facets
# ---------------------------------------------------------------------------
def plot_by_donor(merged: pd.DataFrame) -> None:
    donors = sorted(merged["donor_tx"].dropna().unique())
    n = len(donors)
    n_col = 5
    n_row = int(np.ceil(n / n_col))
    fig, axes = plt.subplots(n_row, n_col, figsize=(n_col * 2.6, n_row * 2.6), squeeze=False)

    for i, donor in enumerate(donors):
        ax = axes[i // n_col][i % n_col]
        sub = merged[merged["donor_tx"] == donor]
        if len(sub) >= 3:
            r_log, _ = pearsonr(np.log1p(sub["total_counts_sr"]), np.log1p(sub["total_counts_tx"]))
            title = f"{donor}\nn={len(sub):,}, r(log)={r_log:.2f}"
        else:
            title = f"{donor}\nn={len(sub):,}"
        gridsize = max(8, min(30, int(np.sqrt(len(sub)))))
        ax.hexbin(sub["total_counts_sr"], sub["total_counts_tx"],
                  gridsize=gridsize, bins="log", cmap="hot",
                  xscale="log", yscale="log", mincnt=1)
        ax.set_title(title, fontsize=7.5)
        ax.tick_params(labelsize=6)
        sns.despine(ax=ax)

    for j in range(n, n_row * n_col):
        axes[j // n_col][j % n_col].axis("off")

    fig.suptitle("SR vs TX depth, by donor", weight="semibold", y=1.01)
    fig.tight_layout()
    _savefig(fig, "06_depth_scatter_by_donor.png")


# ---------------------------------------------------------------------------
# 7. SR depth: cells with vs without a TX match
# ---------------------------------------------------------------------------
def plot_sr_depth_by_match_status(sr_obs: pd.DataFrame, tx_obs: pd.DataFrame) -> None:
    has_match = sr_obs.index.isin(tx_obs.index)
    df = sr_obs.copy()
    df["has_tx_match"] = np.where(has_match, "Matched in TX", "No TX match")

    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.logspace(np.log10(max(df["total_counts"].min(), 1)),
                        np.log10(df["total_counts"].max()), 50)
    for label, color in [("Matched in TX", "#06D6A0"), ("No TX match", "#FB5607")]:
        sub = df[df["has_tx_match"] == label]["total_counts"]
        ax.hist(sub, bins=bins, alpha=0.55, color=color, label=f"{label} (n={len(sub):,})",
                density=True, edgecolor="white", linewidth=0.2)
    ax.set_xscale("log")
    ax.set_xlabel("SR total_counts (log scale)")
    ax.set_ylabel("density")
    ax.set_title("Are the SR cells missing from TX simply the lowest-depth ones?",
                 weight="semibold", pad=12, fontsize=11)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "07_sr_depth_by_tx_match_status.png")


def main():
    sr_obs, tx_obs = load_obs()
    merged = build_merged(sr_obs, tx_obs)

    print(f"\nn SR cells: {len(sr_obs):,}")
    print(f"n TX cells: {len(tx_obs):,}")
    print(f"n matched (inner join): {len(merged):,}")
    print(f"frac TX barcodes found in SR: {len(merged) / len(tx_obs):.3f}")
    print(f"frac SR barcodes found in TX: {len(merged) / len(sr_obs):.3f}\n")

    plot_depth_scatter(merged)
    plot_marginal_histograms(merged)
    plot_ratio_distribution(merged)
    plot_efficiency_vs_depth(merged)
    plot_by_cell_type(merged)
    plot_by_donor(merged)
    plot_sr_depth_by_match_status(sr_obs, tx_obs)

    print(f"\nAll plots saved to: {OUT_PLOTS_DEPTH}")


if __name__ == "__main__":
    main()

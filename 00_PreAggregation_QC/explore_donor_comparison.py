"""
Diagnostic: compare donors to each other on the SR/TX depth metrics from
explore_sr_tx_depth.py, and check whether donor-to-donor variability lines
up with library_kit (a real batch confound -- PO donors = 3'v4, SMC donors
= 3'v3) or with condition (AD / Control / Active control). If it tracks
library_kit, that's a normalization/batch issue. If it tracks condition,
that's a confound that would undermine 01_ViscachaDTU_Analysis's AD-vs-Control comparison.

Standalone exploratory script. Reuses load_obs/build_merged from
explore_sr_tx_depth.py so the underlying data and definitions match exactly.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python -m 00_PreAggregation_QC.explore_donor_comparison
     (from /home/welcome3/Viscacha_pipeline)
Output: outputs/00_PreAggregation_QC/plots/donor_comparison/*.png
        outputs/00_PreAggregation_QC/plots/donor_comparison/donor_summary.csv
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from itertools import combinations
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, mannwhitneyu, kruskal, ks_2samp
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess

from .config import OUT_DIR
from .explore_sr_tx_depth import load_obs, build_merged

OUT_PLOTS_DONOR = OUT_DIR / "plots" / "donor_comparison"

COND_ORDER  = ["Control", "Active control", "AD"]
COND_COLORS = {"AD": "#EF233C", "Control": "#06D6A0", "Active control": "#FFBE0B"}
KIT_COLORS  = {"3'v3": "#3A86FF", "3'v4": "#FB5607"}

METRICS = [
    ("median_sr_depth",   "Median SR depth",          True),
    ("median_tx_depth",   "Median TX depth",          True),
    ("median_efficiency", "Median TX/SR efficiency",  True),
    ("r_log_sr_tx",       "SR-TX log-log correlation", False),
]


def _savefig(fig: plt.Figure, name: str) -> Path:
    OUT_PLOTS_DONOR.mkdir(parents=True, exist_ok=True)
    path = OUT_PLOTS_DONOR / name
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] saved {path.name}")
    return path


def build_donor_summary(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor, sub in merged.groupby("donor_tx"):
        eff = sub["total_counts_tx"] / sub["total_counts_sr"]
        if len(sub) >= 3:
            r_log, _ = pearsonr(np.log1p(sub["total_counts_sr"]), np.log1p(sub["total_counts_tx"]))
        else:
            r_log = np.nan
        rows.append({
            "donor": donor,
            "condition": sub["condition_sr"].iloc[0],
            "library_kit": sub["library_kit_sr"].iloc[0],
            "n_cells_matched": len(sub),
            "median_sr_depth": sub["total_counts_sr"].median(),
            "median_tx_depth": sub["total_counts_tx"].median(),
            "median_efficiency": eff.median(),
            "r_log_sr_tx": r_log,
        })
    return pd.DataFrame(rows).sort_values("donor").reset_index(drop=True)


def _strip_box(ax, df, col, group_col, order, palette, log=False):
    sns.boxplot(data=df, x=group_col, y=col, order=order, hue=group_col, palette=palette,
                ax=ax, width=0.5, fliersize=0, boxprops={"alpha": 0.45}, legend=False)
    sns.stripplot(data=df, x=group_col, y=col, order=order, hue=group_col, palette=palette,
                  ax=ax, size=6, alpha=0.9, edgecolor="white", linewidth=0.6, jitter=0.15,
                  legend=False)
    if log:
        ax.set_yscale("log")
    sns.despine(ax=ax)


# ---------------------------------------------------------------------------
# A. Technical metrics by library_kit (batch check)
# ---------------------------------------------------------------------------
def plot_metrics_by_kit(summary: pd.DataFrame) -> None:
    kits = sorted(summary["library_kit"].dropna().unique())
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    for ax, (col, label, log) in zip(axes, METRICS):
        df = summary.dropna(subset=[col])
        _strip_box(ax, df, col, "library_kit", kits, KIT_COLORS, log=log)
        groups = [df.loc[df["library_kit"] == k, col] for k in kits]
        if all(len(g) > 0 for g in groups) and len(kits) == 2:
            u, p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
            ax.set_title(f"{label}\nMann-Whitney p = {p:.3f}", fontsize=9.5)
        else:
            ax.set_title(label, fontsize=9.5)
        ax.set_xlabel(None)
    fig.suptitle("Technical metrics by library kit (PO=3'v4, SMC=3'v3)", weight="semibold", y=1.04)
    fig.tight_layout()
    _savefig(fig, "A_metrics_by_library_kit.png")


# ---------------------------------------------------------------------------
# B. Technical metrics by condition (biology confound check)
# ---------------------------------------------------------------------------
def plot_metrics_by_condition(summary: pd.DataFrame) -> None:
    conds = [c for c in COND_ORDER if c in summary["condition"].unique()]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    for ax, (col, label, log) in zip(axes, METRICS):
        df = summary.dropna(subset=[col])
        _strip_box(ax, df, col, "condition", conds, COND_COLORS, log=log)
        groups = [df.loc[df["condition"] == c, col] for c in conds]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) >= 2:
            h, p = kruskal(*groups)
            ax.set_title(f"{label}\nKruskal-Wallis p = {p:.3f}", fontsize=9.5)
        else:
            ax.set_title(label, fontsize=9.5)
        ax.set_xlabel(None)
        plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    fig.suptitle("Technical metrics by condition", weight="semibold", y=1.04)
    fig.tight_layout()
    _savefig(fig, "B_metrics_by_condition.png")


# ---------------------------------------------------------------------------
# C. Donor cell yield (n_cells_matched) by kit and by condition
# ---------------------------------------------------------------------------
def plot_ncells_comparisons(summary: pd.DataFrame) -> None:
    kits  = sorted(summary["library_kit"].dropna().unique())
    conds = [c for c in COND_ORDER if c in summary["condition"].unique()]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))

    _strip_box(axes[0], summary, "n_cells_matched", "library_kit", kits, KIT_COLORS, log=True)
    groups = [summary.loc[summary["library_kit"] == k, "n_cells_matched"] for k in kits]
    if len(kits) == 2:
        u, p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
        axes[0].set_title(f"By library kit\nMann-Whitney p = {p:.3f}", fontsize=9.5)
    axes[0].set_xlabel(None)

    _strip_box(axes[1], summary, "n_cells_matched", "condition", conds, COND_COLORS, log=True)
    groups = [summary.loc[summary["condition"] == c, "n_cells_matched"] for c in conds]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        h, p = kruskal(*groups)
        axes[1].set_title(f"By condition\nKruskal-Wallis p = {p:.3f}", fontsize=9.5)
    axes[1].set_xlabel(None)
    plt.setp(axes[1].get_xticklabels(), rotation=20, ha="right")

    fig.suptitle("Donor cell yield (matched SR+TX cells)", weight="semibold", y=1.03)
    fig.tight_layout()
    _savefig(fig, "C_ncells_matched_comparisons.png")


# ---------------------------------------------------------------------------
# D. Donor profile heatmap (z-scored metrics, row-annotated by kit/condition)
# ---------------------------------------------------------------------------
def plot_donor_profile_heatmap(summary: pd.DataFrame) -> None:
    cols = ["median_sr_depth", "median_tx_depth", "median_efficiency",
            "r_log_sr_tx", "n_cells_matched"]
    df = summary.set_index("donor")[cols]
    z = (df - df.mean()) / df.std()
    z = z.loc[summary.sort_values(["library_kit", "condition"])["donor"].values]

    fig, (ax_anno, ax_heat) = plt.subplots(
        1, 2, figsize=(9, 8), gridspec_kw={"width_ratios": [0.5, 6]}
    )

    kit_codes  = summary.set_index("donor").loc[z.index, "library_kit"].map(KIT_COLORS)
    cond_codes = summary.set_index("donor").loc[z.index, "condition"].map(COND_COLORS)
    anno = pd.DataFrame({"kit": kit_codes, "condition": cond_codes})
    for i, col in enumerate(anno.columns):
        ax_anno.bar(i, len(z), bottom=0, width=1, color="none")
    for row_i, donor in enumerate(z.index):
        ax_anno.add_patch(plt.Rectangle((0, row_i), 1, 1, color=anno.loc[donor, "kit"]))
        ax_anno.add_patch(plt.Rectangle((1, row_i), 1, 1, color=anno.loc[donor, "condition"]))
    ax_anno.set_xlim(0, 2)
    ax_anno.set_ylim(0, len(z))
    ax_anno.set_xticks([0.5, 1.5])
    ax_anno.set_xticklabels(["kit", "condition"], rotation=45, ha="right", fontsize=8)
    ax_anno.set_yticks(np.arange(len(z)) + 0.5)
    ax_anno.set_yticklabels(z.index, fontsize=7)
    ax_anno.invert_yaxis()
    for spine in ax_anno.spines.values():
        spine.set_visible(False)

    sns.heatmap(z, ax=ax_heat, cmap="Spectral_r", center=0, cbar_kws={"label": "z-score"},
                yticklabels=False, xticklabels=[c.replace("_", " ") for c in cols])
    ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=30, ha="right")
    ax_heat.set_ylabel(None)

    fig.suptitle("Per-donor technical profile (rows sorted by kit, then condition)",
                 weight="semibold", y=1.0)
    fig.tight_layout()
    _savefig(fig, "D_donor_profile_heatmap.png")


# ---------------------------------------------------------------------------
# E. Donor-level SR vs TX depth scatter, colored by kit, sized by n_cells
# ---------------------------------------------------------------------------
def plot_donor_level_scatter(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    for kit, color in KIT_COLORS.items():
        sub = summary[summary["library_kit"] == kit]
        ax.scatter(sub["median_sr_depth"], sub["median_tx_depth"],
                   s=sub["n_cells_matched"] / 15, color=color, alpha=0.75,
                   edgecolor="white", linewidth=0.8, label=kit)
    for _, row in summary.iterrows():
        ax.annotate(row["donor"], (row["median_sr_depth"], row["median_tx_depth"]),
                    fontsize=6, alpha=0.7, xytext=(3, 3), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Donor median SR depth (log scale)")
    ax.set_ylabel("Donor median TX depth (log scale)")
    ax.set_title("Donor-level depth summary\n(point size = # matched cells)", weight="semibold", pad=12)
    ax.legend(title="library_kit", frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "E_donor_level_depth_scatter.png")


# ---------------------------------------------------------------------------
# F/G/H/I. Distributional comparison: pooled KS tests + ECDF overlays.
# No per-donor summarization -- compares the actual cell-level distributions.
# ---------------------------------------------------------------------------
def compute_ks_tests(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in ["total_counts_sr", "total_counts_tx", "efficiency"]:
        kits = sorted(merged["library_kit_sr"].dropna().unique())
        a = merged.loc[merged["library_kit_sr"] == kits[0], col]
        b = merged.loc[merged["library_kit_sr"] == kits[1], col]
        stat, p = ks_2samp(a, b)
        rows.append({"metric": col, "comparison": f"{kits[0]} vs {kits[1]}",
                     "group_by": "library_kit", "n1": len(a), "n2": len(b),
                     "ks_stat": stat, "p_value": p})

        conds = [c for c in COND_ORDER if c in merged["condition_sr"].unique()]
        for c1, c2 in combinations(conds, 2):
            a = merged.loc[merged["condition_sr"] == c1, col]
            b = merged.loc[merged["condition_sr"] == c2, col]
            stat, p = ks_2samp(a, b)
            rows.append({"metric": col, "comparison": f"{c1} vs {c2}",
                         "group_by": "condition", "n1": len(a), "n2": len(b),
                         "ks_stat": stat, "p_value": p})
    return pd.DataFrame(rows)


def _ecdf_overlay(merged: pd.DataFrame, value_col: str, group_col: str,
                   group_colors: dict, group_order: list, xlabel: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 6.5))

    # Thin per-donor ECDFs -- the actual distributions, no summarization
    for donor, sub in merged.groupby("donor_tx"):
        grp = sub[group_col].iloc[0]
        vals = np.sort(sub[value_col].values)
        y = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, y, color=group_colors.get(grp, "gray"), alpha=0.35, linewidth=1.0)

    # Bold pooled-group ECDFs on top
    for grp in group_order:
        vals = np.sort(merged.loc[merged[group_col] == grp, value_col].values)
        if len(vals) == 0:
            continue
        y = np.arange(1, len(vals) + 1) / len(vals)
        ax.plot(vals, y, color=group_colors[grp], linewidth=3.2,
                label=f"{grp} (pooled, n={len(vals):,})")

    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("ECDF")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    sns.despine(ax=ax)
    fig.tight_layout()
    return fig


def plot_ecdf_by_kit(merged: pd.DataFrame, ks: pd.DataFrame) -> None:
    kits = sorted(merged["library_kit_sr"].dropna().unique())
    for col, label, fname in [
        ("total_counts_sr", "SR total_counts (log scale)", "F_ecdf_sr_depth_by_kit.png"),
        ("total_counts_tx", "TX total_counts (log scale)", "G_ecdf_tx_depth_by_kit.png"),
        ("efficiency", "TX/SR efficiency per cell (log scale)", "N_ecdf_efficiency_by_kit.png"),
    ]:
        fig = _ecdf_overlay(merged, col, "library_kit_sr", KIT_COLORS, kits, label)
        row = ks[(ks["metric"] == col) & (ks["group_by"] == "library_kit")].iloc[0]
        fig.axes[0].set_title(
            f"Per-donor ECDFs ({col}), by library kit\n"
            f"pooled KS: D={row['ks_stat']:.3f}, p={row['p_value']:.2e}",
            weight="semibold", pad=12)
        _savefig(fig, fname)


def plot_ecdf_by_condition(merged: pd.DataFrame, ks: pd.DataFrame) -> None:
    conds = [c for c in COND_ORDER if c in merged["condition_sr"].unique()]
    for col, label, fname in [
        ("total_counts_sr", "SR total_counts (log scale)", "H_ecdf_sr_depth_by_condition.png"),
        ("total_counts_tx", "TX total_counts (log scale)", "I_ecdf_tx_depth_by_condition.png"),
        ("efficiency", "TX/SR efficiency per cell (log scale)", "O_ecdf_efficiency_by_condition.png"),
    ]:
        fig = _ecdf_overlay(merged, col, "condition_sr", COND_COLORS, conds, label)
        sub_ks = ks[(ks["metric"] == col) & (ks["group_by"] == "condition")]
        ks_lines = "\n".join(
            f"{r.comparison}: D={r.ks_stat:.3f}, p={r.p_value:.2e}" for r in sub_ks.itertuples()
        )
        fig.axes[0].set_title(
            f"Per-donor ECDFs ({col}), by condition\n{ks_lines}",
            weight="semibold", pad=12, fontsize=10)
        _savefig(fig, fname)


# ---------------------------------------------------------------------------
# J. Per-cell SR->TX trend (LOWESS), by kit. Thin lines = each donor's own
# fitted trend through their actual per-cell points (no summarization);
# bold lines = pooled trend per kit group. Line width/alpha ~ donor n_cells.
# ---------------------------------------------------------------------------
def plot_lowess_by_kit(merged: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 7))

    n_by_donor = merged.groupby("donor_tx").size()
    max_n = n_by_donor.max()

    for donor, sub in merged.groupby("donor_tx"):
        if len(sub) < 10:
            continue
        kit = sub["library_kit_sr"].iloc[0]
        x = np.log10(sub["total_counts_sr"].values)
        y = np.log10(sub["total_counts_tx"].values)
        fit = sm_lowess(y, x, frac=0.4, delta=0.01 * (x.max() - x.min()))
        n = len(sub)
        lw    = 0.7 + 1.6 * (n / max_n)
        alpha = 0.25 + 0.45 * (n / max_n)
        ax.plot(10 ** fit[:, 0], 10 ** fit[:, 1],
                color=KIT_COLORS.get(kit, "gray"), linewidth=lw, alpha=alpha)

    for kit, color in KIT_COLORS.items():
        sub = merged[merged["library_kit_sr"] == kit]
        x = np.log10(sub["total_counts_sr"].values)
        y = np.log10(sub["total_counts_tx"].values)
        fit = sm_lowess(y, x, frac=0.2, delta=0.01 * (x.max() - x.min()))
        ax.plot(10 ** fit[:, 0], 10 ** fit[:, 1], color=color, linewidth=3.5,
                label=f"{kit} (pooled, n={len(sub):,})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("SR total_counts (log scale)")
    ax.set_ylabel("TX total_counts (log scale)")
    ax.set_title(
        "Per-donor SR → TX trend (LOWESS), by library kit\n"
        "thin = individual donor (width/alpha ∝ n cells)  |  bold = pooled per kit",
        weight="semibold", fontsize=10.5, pad=12)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "J_lowess_trend_by_kit.png")


# ---------------------------------------------------------------------------
# K. 2D density contours of (SR, TX) depth, by kit -- joint distribution
# comparison without binning or summarizing to a single value. Filled +
# marginal density panels (joint-grid style) for readability.
# ---------------------------------------------------------------------------
def plot_kde_contours_by_kit(merged: pd.DataFrame) -> None:
    df = merged.copy()
    df["log_sr"] = np.log10(df["total_counts_sr"])
    df["log_tx"] = np.log10(df["total_counts_tx"])

    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(4, 4, hspace=0.08, wspace=0.08)
    ax_main  = fig.add_subplot(gs[1:4, 0:3])
    ax_top   = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    for kit, color in KIT_COLORS.items():
        sub = df[df["library_kit_sr"] == kit]
        sns.kdeplot(data=sub, x="log_sr", y="log_tx", ax=ax_main, fill=True,
                    color=color, alpha=0.45, levels=4, thresh=0.05)
        sns.kdeplot(data=sub, x="log_sr", y="log_tx", ax=ax_main, fill=False,
                    color=color, levels=4, thresh=0.05, linewidths=1.2, alpha=0.9)
        sns.kdeplot(data=sub, x="log_sr", ax=ax_top, color=color, fill=True, alpha=0.35, linewidth=1.5)
        sns.kdeplot(data=sub, y="log_tx", ax=ax_right, color=color, fill=True, alpha=0.35, linewidth=1.5)

    handles = [Line2D([0], [0], color=c, lw=2.5, label=k) for k, c in KIT_COLORS.items()]
    ax_main.legend(handles=handles, frameon=False, title="library_kit", loc="upper left")
    ax_main.set_xlabel("log10(SR total_counts)")
    ax_main.set_ylabel("log10(TX total_counts)")
    sns.despine(ax=ax_main)

    for ax in (ax_top, ax_right):
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.suptitle("2D density: SR vs TX depth, by library kit\nfilled = density, marginals = 1D depth distributions",
                 weight="semibold", y=1.0)
    _savefig(fig, "K_kde_contours_by_kit.png")


# ---------------------------------------------------------------------------
# L. Density DIFFERENCE between kits -- directly shows where one kit's
# density exceeds the other's, instead of requiring a visual diff of overlays.
# ---------------------------------------------------------------------------
def plot_kde_difference_by_kit(merged: pd.DataFrame) -> None:
    from scipy.stats import gaussian_kde

    df = merged.copy()
    df["log_sr"] = np.log10(df["total_counts_sr"])
    df["log_tx"] = np.log10(df["total_counts_tx"])

    xmin, xmax = df["log_sr"].min(), df["log_sr"].max()
    ymin, ymax = df["log_tx"].min(), df["log_tx"].max()
    xx, yy = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
    grid = np.vstack([xx.ravel(), yy.ravel()])

    kits = list(KIT_COLORS.keys())  # ["3'v3", "3'v4"]
    densities = {}
    for kit in kits:
        sub = df[df["library_kit_sr"] == kit]
        kde = gaussian_kde(np.vstack([sub["log_sr"], sub["log_tx"]]))
        dens = kde(grid).reshape(xx.shape)
        densities[kit] = dens / dens.sum()  # normalize so both kits have equal total mass

    diff = densities[kits[1]] - densities[kits[0]]
    vmax = np.abs(diff).max()

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    cf = ax.contourf(xx, yy, diff, levels=20, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    fig.colorbar(cf, ax=ax, label=f"density({kits[1]}) − density({kits[0]})")
    ax.set_xlabel("log10(SR total_counts)")
    ax.set_ylabel("log10(TX total_counts)")
    ax.set_title(f"Where does {kits[1]} (red) exceed {kits[0]} (blue) in density?",
                 weight="semibold", fontsize=11, pad=12)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "L_kde_difference_by_kit.png")


# ---------------------------------------------------------------------------
# M. Per-donor contours -- one thin single-level outline per donor (colored
# by kit, width/alpha ~ n cells), overlaid on faint pooled-kit fills for
# context. Mirrors the LOWESS-by-kit thin-line convention.
# ---------------------------------------------------------------------------
def plot_kde_contours_by_donor(merged: pd.DataFrame) -> None:
    df = merged.copy()
    df["log_sr"] = np.log10(df["total_counts_sr"])
    df["log_tx"] = np.log10(df["total_counts_tx"])

    fig, ax = plt.subplots(figsize=(8.5, 7.5))

    for kit, color in KIT_COLORS.items():
        sub = df[df["library_kit_sr"] == kit]
        sns.kdeplot(data=sub, x="log_sr", y="log_tx", ax=ax, fill=True,
                    color=color, alpha=0.12, levels=3, thresh=0.1)

    n_by_donor = df.groupby("donor_tx").size()
    max_n = n_by_donor.max()
    for donor, sub in df.groupby("donor_tx"):
        if len(sub) < 10:
            continue
        kit = sub["library_kit_sr"].iloc[0]
        n = len(sub)
        lw    = 0.8 + 1.4 * (n / max_n)
        alpha = 0.35 + 0.5 * (n / max_n)
        sns.kdeplot(data=sub, x="log_sr", y="log_tx", ax=ax, fill=False,
                    color=KIT_COLORS.get(kit, "gray"), levels=2, thresh=0.4,
                    linewidths=lw, alpha=alpha)

    handles = [Line2D([0], [0], color=c, lw=2.5, label=k) for k, c in KIT_COLORS.items()]
    ax.legend(handles=handles, frameon=False, title="library_kit", loc="upper left")
    ax.set_xlabel("log10(SR total_counts)")
    ax.set_ylabel("log10(TX total_counts)")
    ax.set_title(
        "Per-donor density contours (outer-mass outline), colored by kit\n"
        "faint fill = pooled per-kit density  |  line width/alpha ∝ donor n cells",
        weight="semibold", fontsize=10.5, pad=12)
    sns.despine(ax=ax)
    fig.tight_layout()
    _savefig(fig, "M_kde_contours_by_donor.png")


# ---------------------------------------------------------------------------
# N2. Ridgeline (joyplot): one density curve per donor, stacked, ordered by
# kit then median efficiency. No specialized package available (joypy not
# installed) -- built manually with gaussian_kde.
# ---------------------------------------------------------------------------
def plot_ridgeline_efficiency_by_donor(merged: pd.DataFrame) -> None:
    from scipy.stats import gaussian_kde

    df = merged.copy()
    df["log_efficiency"] = np.log10(df["efficiency"])

    order_df = (df.groupby("donor_tx")
                  .agg(kit=("library_kit_sr", "first"), med=("log_efficiency", "median"))
                  .reset_index()
                  .sort_values(["kit", "med"]))
    donors_ordered = order_df["donor_tx"].tolist()
    n = len(donors_ordered)

    xmin, xmax = df["log_efficiency"].min(), df["log_efficiency"].max()
    xs = np.linspace(xmin - 0.3, xmax + 0.3, 400)  # KDE grid, in log10 units (bandwidth needs log space)
    xs_plot = 10 ** xs  # actual efficiency values, for a true log-scale x-axis

    fig, ax = plt.subplots(figsize=(8.5, 0.42 * n + 1.5))
    step = 1.0
    yticks, yticklabels = [], []
    for i, donor in enumerate(donors_ordered):
        sub = df[df["donor_tx"] == donor]
        kit = sub["library_kit_sr"].iloc[0]
        color = KIT_COLORS.get(kit, "gray")
        kde = gaussian_kde(sub["log_efficiency"].values)
        density = kde(xs)
        density = density / density.max() * 1.3  # slight overlap into the row above
        y_base = i * step
        ax.fill_between(xs_plot, y_base, y_base + density, color=color, alpha=0.8,
                         linewidth=0.8, edgecolor="white", zorder=n - i)
        yticks.append(y_base)
        yticklabels.append(donor)

    ax.set_xscale("log")
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=7)
    ax.set_xlabel("TX/SR efficiency per cell (log scale)")
    ax.set_xlim(10 ** (xmin - 0.3), 10 ** (xmax + 0.3))
    ax.set_ylim(-0.2, n * step + 1.0)
    ax.xaxis.grid(True, which="major", color="#cccccc", linewidth=0.8, zorder=0)
    ax.xaxis.grid(True, which="minor", color="#e8e8e8", linewidth=0.5, zorder=0)
    handles = [Line2D([0], [0], color=c, lw=4, label=k) for k, c in KIT_COLORS.items()]
    ax.legend(handles=handles, frameon=False, loc="upper right", title="library_kit")
    ax.set_title("Per-donor efficiency distributions (ridgeline)\nordered by kit, then by median within kit",
                 weight="semibold", pad=12)
    sns.despine(ax=ax, left=True)
    fig.tight_layout()
    _savefig(fig, "P_ridgeline_efficiency_by_donor.png")


# ---------------------------------------------------------------------------
# O2. Raincloud plots (violin + boxplot + jittered subsample) of per-cell
# efficiency, by kit and by condition. Points are subsampled per group
# (max 400) purely for jitter-plot legibility -- the violin/box use all cells.
# ---------------------------------------------------------------------------
def _raincloud(ax, df, group_col, order, palette, rng):
    sns.violinplot(data=df, x=group_col, y="log_efficiency", order=order,
                    hue=group_col, palette=palette, inner=None, cut=0,
                    ax=ax, legend=False, alpha=0.45, linewidth=1.0)
    sns.boxplot(data=df, x=group_col, y="log_efficiency", order=order,
                hue=group_col, palette=palette, ax=ax, width=0.12, legend=False,
                boxprops={"zorder": 3, "alpha": 0.9}, fliersize=0,
                whiskerprops={"zorder": 3}, capprops={"zorder": 3},
                medianprops={"zorder": 3, "color": "black"})
    sample = df.groupby(group_col, group_keys=False).apply(
        lambda g: g.sample(min(400, len(g)), random_state=42), include_groups=False
    )
    sample[group_col] = df.loc[sample.index, group_col]
    sns.stripplot(data=sample, x=group_col, y="log_efficiency", order=order,
                  hue=group_col, palette=palette, ax=ax, size=2.5, alpha=0.35,
                  jitter=0.25, legend=False)
    ax.set_xlabel(None)
    ax.set_ylabel("log10(TX/SR efficiency per cell)")
    sns.despine(ax=ax)


def plot_raincloud_efficiency(merged: pd.DataFrame) -> None:
    df = merged.copy()
    df["log_efficiency"] = np.log10(df["efficiency"])
    rng = np.random.default_rng(42)

    kits = sorted(df["library_kit_sr"].dropna().unique())
    fig, ax = plt.subplots(figsize=(7, 6))
    _raincloud(ax, df, "library_kit_sr", kits, KIT_COLORS, rng)
    ax.set_title("Efficiency by library kit\n(violin + box + 400 sampled points per group)",
                 weight="semibold", pad=12)
    fig.tight_layout()
    _savefig(fig, "Q_raincloud_efficiency_by_kit.png")

    conds = [c for c in COND_ORDER if c in df["condition_sr"].unique()]
    fig, ax = plt.subplots(figsize=(8, 6))
    _raincloud(ax, df, "condition_sr", conds, COND_COLORS, rng)
    ax.set_title("Efficiency by condition\n(violin + box + 400 sampled points per group)",
                 weight="semibold", pad=12)
    plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
    fig.tight_layout()
    _savefig(fig, "R_raincloud_efficiency_by_condition.png")


def main():
    sr_obs, tx_obs = load_obs()
    merged = build_merged(sr_obs, tx_obs)
    merged["efficiency"] = merged["total_counts_tx"] / merged["total_counts_sr"]
    summary = build_donor_summary(merged)

    OUT_PLOTS_DONOR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_PLOTS_DONOR / "donor_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"\nDonor summary ({len(summary)} donors) written to: {csv_path}\n")
    print(summary.to_string(index=False))
    print()

    plot_metrics_by_kit(summary)
    plot_metrics_by_condition(summary)
    plot_ncells_comparisons(summary)
    plot_donor_profile_heatmap(summary)
    plot_donor_level_scatter(summary)

    ks = compute_ks_tests(merged)
    ks_path = OUT_PLOTS_DONOR / "ks_tests.csv"
    ks.to_csv(ks_path, index=False)
    print(f"\nKS test results written to: {ks_path}\n")
    print(ks.to_string(index=False))

    plot_ecdf_by_kit(merged, ks)
    plot_ecdf_by_condition(merged, ks)

    plot_lowess_by_kit(merged)
    plot_kde_contours_by_kit(merged)
    plot_kde_difference_by_kit(merged)
    plot_kde_contours_by_donor(merged)

    plot_ridgeline_efficiency_by_donor(merged)
    plot_raincloud_efficiency(merged)

    print(f"\nAll plots saved to: {OUT_PLOTS_DONOR}")


if __name__ == "__main__":
    main()

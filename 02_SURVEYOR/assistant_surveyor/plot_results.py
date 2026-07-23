"""
ASSISTANT_SURVEYOR visualization suite — 18 diagnostic plots.

Run from the repo root:
    python assistant_surveyor/plot_results.py
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
import seaborn as sns

# ── paths ──────────────────────────────────────────────────────────────────
HITS_CSV = Path("outputs/assistant_surveyor/hits_enriched.csv")
OUT_DIR  = Path("outputs/assistant_surveyor/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── neon colormap (same as expression_vs_psi plots) ───────────────────────
NEON_CMAP = LinearSegmentedColormap.from_list(
    "neon", ["#0a0020", "#5522cc", "#cc1188", "#cc5500", "#bbbb00"]
)

def _neon(t: float) -> str:
    """Sample a hex color from the neon cmap at position t ∈ [0, 1]."""
    return mcolors.to_hex(NEON_CMAP(t))

# ── consistent colour palettes derived from the neon cmap ─────────────────
PROXY_COLORS = {
    "C":   "#cc1155",   # crimson-rose      ~330° — CDS + domain
    "D":   "#0099bb",   # teal-cyan         ~195° — CDS, no domain
    "NMD": "#cc8800",   # amber             ~ 40° — NMD isoform
    "N":   "#999999",   # neutral grey              — no protein change
}
GATE_COLORS  = {True: "#17a2a8", False: "#cccccc"}   # junior_pass: pass (turquoise) / drop
BIO_COLORS   = {
    "PC_CDS":    "#146b3a",   # dark green   — protein-coding, CDS-altering
    "PC_UTR":    "#3fa66d",   # medium green — protein-coding, UTR-only
    "PC_CDS_ND": "#8fd9b6",   # light green  — protein-coding, CDS not defined
    "novel":     _neon(0.75), # neon orange  — unchanged
    "RI":        "#4a1f7a",   # dark purple  — retained intron
    "NMD":       "#7c4aa8",   # medium purple — nonsense-mediated decay
    "TEC":       "#a97fc9",   # light purple — to be experimentally confirmed
    "other":     "#d4bce6",   # lightest purple — dropped, minor GENCODE biotypes
}
BIO_ORDER = ["PC_CDS", "PC_UTR", "PC_CDS_ND", "novel", "TEC", "RI", "NMD", "other"]
OT_COLORS = {"supported": _neon(0.50), "emerging": _neon(0.75), "novel": "#cccccc"}
GROUP_COLORS = {
    "trial_failure_candidate": "#2a78d6",   # blue — dominant/canonical, down in AD
    "new_target_candidate":    "#e34948",   # red  — minor/alternate, up in AD
}

ANNOT_COLOR  = "#003399"
THRESH_COLOR = "red"
ZERO_COLOR   = "black"

CELL_ORDER = [
    "Excitatory_neuron", "Inhibitory_neuron", "Astrocyte",
    "Oligodendrocyte", "OPC", "Microglia", "Vascular_cell", "Lymphocyte",
]

BG   = "#ffffff"
FG   = "#1a1a1a"
GRID = "#e0e0e0"

plt.rcParams.update({
    "font.family":          "sans-serif",
    "font.size":            10,
    "text.color":           FG,
    "axes.facecolor":       BG,
    "figure.facecolor":     BG,
    "axes.edgecolor":       "#cccccc",
    "axes.labelcolor":      FG,
    "xtick.color":          FG,
    "ytick.color":          FG,
    "axes.spines.top":      False,
    "axes.spines.right":    False,
    "axes.spines.left":     True,
    "axes.spines.bottom":   True,
    "axes.grid":            False,
    "legend.facecolor":     BG,
    "legend.edgecolor":     "#cccccc",
    "legend.labelcolor":    FG,
    "figure.dpi":           150,
})


def save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  saved → {path.name}")


def load() -> pd.DataFrame:
    df = pd.read_csv(HITS_CSV)
    df["abs_delta"]      = df["delta_usage"].abs()
    df["neg_log10_pval"] = -np.log10(df["permutation_pval"].clip(lower=1e-4))
    return df


def push_apart_y(ys: list[float], min_gap: float = 0.25,
                 ylim: tuple[float, float] = (0.0, 4.5)) -> list[float]:
    """Iteratively push label y-positions apart so annotations don't overlap."""
    ys = list(ys)
    for _ in range(300):
        moved = False
        order = sorted(range(len(ys)), key=lambda i: ys[i])
        for k in range(len(order) - 1):
            i, j = order[k], order[k + 1]
            gap = ys[j] - ys[i]
            if gap < min_gap:
                push   = (min_gap - gap) / 2
                ys[i] -= push
                ys[j] += push
                moved  = True
        if not moved:
            break
    margin = 0.05
    return [max(ylim[0] + margin, min(ylim[1] - margin, y)) for y in ys]


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 1 — OVERVIEW
# ═══════════════════════════════════════════════════════════════════════════

def plot_junior_gate_overview(df: pd.DataFrame) -> None:
    """01 — Stacked bar: hits per junior-layer gate outcome coloured by proxy_type."""
    proxy_order = ["C", "D", "NMD", "N"]
    gate_order  = [True, False]
    gate_labels = {True: "Pass", False: "Drop"}

    counts = (
        df.groupby(["junior_pass", "proxy_type"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=gate_order, columns=proxy_order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(7, 4))
    bottoms = np.zeros(len(gate_order))
    for pt in proxy_order:
        vals = counts[pt].values
        bars = ax.bar(
            [gate_labels[g] for g in gate_order], vals,
            bottom=bottoms, color=PROXY_COLORS[pt], label=pt,
            edgecolor=BG, linewidth=0.5,
        )
        for bar, val, bot in zip(bars, vals, bottoms):
            if val > 20:
                ax.text(
                    bar.get_x() + bar.get_width() / 2, bot + val / 2,
                    str(val), ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold",
                )
        bottoms += vals

    for i, tot in enumerate(counts.sum(axis=1)):
        ax.text(i, tot + 10, str(tot), ha="center", va="bottom",
                fontsize=10, fontweight="bold")

    ax.set_ylabel("Hits (transcript × cell-type)")
    ax.set_title("Junior-layer gate outcome by proxy type", fontweight="bold")
    ax.legend(title="proxy_type", bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "01_junior_gate_overview.png")


def plot_biotype_by_celltype(df: pd.DataFrame) -> None:
    """02 — Stacked bar: biotype_class per cell type."""
    bio_order = BIO_ORDER
    cells     = [c for c in CELL_ORDER if c in df["cell_type"].unique()]

    counts = (
        df.groupby(["cell_type", "biotype_class"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=cells, columns=bio_order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    bottoms = np.zeros(len(cells))
    for bc in bio_order:
        vals = counts[bc].values
        ax.bar(cells, vals, bottom=bottoms, color=BIO_COLORS[bc],
               label=bc, edgecolor=BG, linewidth=0.5)
        bottoms += vals

    ax.set_xticklabels(cells, rotation=35, ha="right")
    ax.set_ylabel("Hits")
    ax.set_title("Biotype class distribution per cell type", fontweight="bold")
    ax.legend(title="biotype_class", bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "02_biotype_by_celltype.png")


def plot_proxy_by_celltype(df: pd.DataFrame) -> None:
    """03 — Stacked bar: proxy_type per cell type."""
    proxy_order = ["C", "D", "NMD", "N"]
    cells       = [c for c in CELL_ORDER if c in df["cell_type"].unique()]

    counts = (
        df.groupby(["cell_type", "proxy_type"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=cells, columns=proxy_order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    bottoms = np.zeros(len(cells))
    for pt in proxy_order:
        vals = counts[pt].values
        bars = ax.bar(cells, vals, bottom=bottoms, color=PROXY_COLORS[pt],
                      label=pt, edgecolor=BG, linewidth=0.5)
        for bar, val, bot in zip(bars, vals, bottoms):
            if val > 15:
                ax.text(
                    bar.get_x() + bar.get_width() / 2, bot + val / 2,
                    str(val), ha="center", va="center", fontsize=8, color="white",
                )
        bottoms += vals

    ax.set_xticklabels(cells, rotation=35, ha="right")
    ax.set_ylabel("Hits")
    ax.set_title("Proxy type distribution per cell type", fontweight="bold")
    ax.legend(title="proxy_type", bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "03_proxy_by_celltype.png")


def plot_usage_direction(df: pd.DataFrame) -> None:
    """04 — Diverging bar: AD-enriched up, CT-enriched down, per cell type."""
    cells     = [c for c in CELL_ORDER if c in df["cell_type"].unique()]
    ad_counts = (df[df["usage_direction"] == "AD_enriched"]
                 .groupby("cell_type").size().reindex(cells, fill_value=0))
    ct_counts = (df[df["usage_direction"] == "CT_enriched"]
                 .groupby("cell_type").size().reindex(cells, fill_value=0))

    AD_COLOR = _neon(0.50)   # muted magenta
    CT_COLOR = _neon(0.22)   # muted violet

    x    = np.arange(len(cells))
    ymax = max(ad_counts.max(), ct_counts.max()) * 1.18

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.bar(x,  ad_counts.values,  color=AD_COLOR, edgecolor=BG, label="AD-enriched",  width=0.6)
    ax.bar(x, -ct_counts.values,  color=CT_COLOR, edgecolor=BG, label="CT-enriched",  width=0.6)

    ax.axhline(0, color=ZERO_COLOR, lw=0.8, alpha=0.5)

    # count labels
    for i, (ad, ct) in enumerate(zip(ad_counts.values, ct_counts.values)):
        if ad > 0:
            ax.text(i, ad + ymax * 0.02, str(ad), ha="center", va="bottom",
                    fontsize=8, color=FG)
        if ct > 0:
            ax.text(i, -(ct + ymax * 0.02), str(ct), ha="center", va="top",
                    fontsize=8, color=FG)

    # direction labels in chart area
    ax.text(len(cells) - 0.5, ymax * 0.88, "AD-enriched ↑",
            ha="right", color=AD_COLOR, fontsize=10, fontweight="bold")
    ax.text(len(cells) - 0.5, -ymax * 0.88, "↓ CT-enriched",
            ha="right", color=CT_COLOR, fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(cells, rotation=35, ha="right")
    ax.set_ylim(-ymax, ymax)

    # show absolute values on y axis
    yticks = [t for t in ax.get_yticks() if abs(t) <= ymax]
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(abs(int(t))) if t == int(t) else "" for t in yticks])
    ax.set_ylabel("Hits")

    ax.set_title("Isoform usage direction per cell type", fontweight="bold")
    save(fig, "04_usage_direction_by_celltype.png")


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 2 — EFFECT SIZE & SIGNIFICANCE
# ═══════════════════════════════════════════════════════════════════════════

def plot_volcano(df: pd.DataFrame) -> None:
    """05 — Volcano styled after expression_vs_psi: colour = proxy_type, black edge = junior pass."""
    ymax = df["neg_log10_pval"].max() * 1.08

    fig, ax = plt.subplots(figsize=(9, 7))

    # Layer 1 — N hits: grey background, small, no edge
    sub_n = df[df["proxy_type"] == "N"]
    ax.scatter(
        sub_n["delta_usage"], sub_n["neg_log10_pval"],
        c="#cccccc", s=6, alpha=0.35, linewidths=0,
        rasterized=True, zorder=1,
    )

    # Layer 2 — NMD / D / C (not junior-pass): coloured, medium, no edge
    for pt in ["NMD", "D", "C"]:
        sub = df[(df["proxy_type"] == pt) & (~df["junior_pass"])]
        ax.scatter(
            sub["delta_usage"], sub["neg_log10_pval"],
            c=PROXY_COLORS[pt], s=14, alpha=0.75, linewidths=0,
            rasterized=True, zorder=2,
        )

    # Layer 3 — junior pass: coloured, larger, black edge
    sub_t1 = df[df["junior_pass"]]
    ax.scatter(
        sub_t1["delta_usage"], sub_t1["neg_log10_pval"],
        c=[PROXY_COLORS[pt] for pt in sub_t1["proxy_type"]],
        s=38, alpha=1.0, linewidths=0.7, edgecolors="black",
        rasterized=True, zorder=3,
    )

    # Reference lines
    ax.axhline(0,                   color=ZERO_COLOR,   lw=0.6, ls="--", alpha=0.35)
    ax.axvline(0,                   color=ZERO_COLOR,   lw=0.6, ls="--", alpha=0.35)
    ax.axvline( 0.25,               color=THRESH_COLOR, lw=0.8, ls="--", alpha=0.7)
    ax.axvline(-0.25,               color=THRESH_COLOR, lw=0.8, ls="--", alpha=0.7)
    ax.axhline(-np.log10(0.05),     color=ZERO_COLOR,   lw=0.6, ls="--", alpha=0.35)

    # Annotate top-6 junior-pass C hits, push labels apart vertically
    top = (
        df[df["junior_pass"] & (df["proxy_type"] == "C")]
        .drop_duplicates("gene_name")
        .nlargest(6, "abs_delta")
    )
    if len(top):
        pts_x    = top["delta_usage"].tolist()
        pts_y    = top["neg_log10_pval"].tolist()
        names    = top["gene_name"].tolist()
        label_ys = push_apart_y(pts_y, min_gap=0.30, ylim=(0.0, ymax))

        for px, py, ly, name in zip(pts_x, pts_y, label_ys, names):
            dx = 0.10 if px >= 0 else -0.10
            ax.annotate(
                name,
                xy=(px, py), xytext=(px + dx, ly),
                fontsize=7.5, color=ANNOT_COLOR,
                arrowprops=dict(arrowstyle="-", color=ANNOT_COLOR, lw=0.7),
                ha="left" if dx > 0 else "right",
                va="center", zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=ANNOT_COLOR, linewidth=0.7, alpha=0.9),
            )

    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Δ usage (AD − Control)", fontsize=9)
    ax.set_ylabel("−log₁₀ (permutation p-value)", fontsize=9)
    ax.set_title(
        "Δ usage vs significance  |  colour = proxy type  |  black edge = junior pass",
        fontsize=10, pad=6,
    )
    ax.tick_params(labelsize=8)

    proxy_patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
                     for pt in ["C", "D", "NMD", "N"]]
    pass_patch    = mpatches.Patch(facecolor="white", edgecolor="black",
                                   label="Junior pass", linewidth=0.8)
    ax.legend(handles=proxy_patches + [pass_patch],
              title="proxy_type", bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)

    save(fig, "05_volcano.png")


def plot_psi_scatter(df: pd.DataFrame) -> None:
    """06 — PSI scatter: AD vs Control coloured by junior-layer gate outcome."""
    fig, ax = plt.subplots(figsize=(7, 7))
    for gate, label in [(False, "Drop"), (True, "Pass")]:
        sub = df[df["junior_pass"] == gate]
        ax.scatter(
            sub["Control"], sub["AD"],
            c=GATE_COLORS[gate], s=10, alpha=0.5, linewidths=0,
            label=label, rasterized=True,
        )
    ax.plot([0, 1], [0, 1], "--", color=ZERO_COLOR, lw=0.6, alpha=0.4)
    ax.set_xlabel("PSI — Control")
    ax.set_ylabel("PSI — AD")
    ax.set_title("PSI scatter: AD vs Control", fontweight="bold")
    ax.legend(title="Junior gate", markerscale=2)
    ax.set_aspect("equal")
    save(fig, "06_psi_scatter.png")


def plot_delta_violin(df: pd.DataFrame) -> None:
    """07 — Violin: |delta_usage| by proxy_type."""
    proxy_order = ["C", "D", "NMD", "N"]
    data   = [df[df["proxy_type"] == pt]["abs_delta"].values for pt in proxy_order]
    colors = [PROXY_COLORS[pt] for pt in proxy_order]

    fig, ax = plt.subplots(figsize=(7, 5))
    parts = ax.violinplot(data, positions=range(4), showmedians=True, showextrema=False)
    for pc, col in zip(parts["bodies"], colors):
        pc.set_facecolor(col)
        pc.set_alpha(0.75)
    parts["cmedians"].set_colors("white")
    parts["cmedians"].set_linewidth(2)

    ax.set_xticks(range(4))
    ax.set_xticklabels(proxy_order)
    ax.set_ylabel("|Δ usage|")
    ax.set_title("|Δ usage| distribution by proxy type", fontweight="bold")
    save(fig, "07_delta_usage_violin.png")


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 3 — OPENTARGETS & EVIDENCE
# ═══════════════════════════════════════════════════════════════════════════

def plot_ot_vs_effect(df: pd.DataFrame) -> None:
    """08 — Scatter: |delta_usage| vs OT score, size = ct_count."""
    proxy_order = ["N", "NMD", "D", "C"]
    fig, ax = plt.subplots(figsize=(8, 6))

    for pt in proxy_order:
        sub = df[df["proxy_type"] == pt]
        ax.scatter(
            sub["abs_delta"], sub["ad_ot_score"],
            c=PROXY_COLORS[pt], s=sub["ct_count"] * 15 + 5,
            alpha=0.55, linewidths=0, label=pt, rasterized=True,
        )

    for n_ct in [1, 2, 3]:
        ax.scatter([], [], c="grey", s=n_ct * 15 + 5,
                   label=f"ct_count={n_ct}", alpha=0.7)

    ax.set_xlabel("|Δ usage|")
    ax.set_ylabel("AD OpenTargets score")
    ax.set_title("Effect size vs OpenTargets AD association", fontweight="bold")
    ax.legend(title="proxy_type / ct", bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "08_ot_vs_effect.png")


def plot_ot_specificity(df: pd.DataFrame) -> None:
    """09 — OT score vs AD specificity. Top-4 OT (navy) and top-4 specificity (teal) annotated."""
    SPEC_CAP   = 35.0
    ANNOT_OT   = "#003399"   # navy  — top OT score
    ANNOT_SPEC = "#007755"   # teal  — top AD specificity

    df_ot = df[df["ad_ot_score"] > 0].copy()
    df_ot["spec_clipped"] = df_ot["ad_specificity"].clip(upper=SPEC_CAP)

    passing  = df[df["junior_pass"]].drop_duplicates("gene_name")
    top_ot   = passing.nlargest(4, "ad_ot_score").copy()
    top_spec = passing[passing["ad_ot_score"] > 0].nlargest(4, "ad_specificity").copy()
    top_ot  ["spec_clipped"] = top_ot  ["ad_specificity"].clip(upper=SPEC_CAP)
    top_spec["spec_clipped"] = top_spec["ad_specificity"].clip(upper=SPEC_CAP)

    fig, ax = plt.subplots(figsize=(9, 7))

    # Layer 1 — N proxy_type: grey background
    sub_n = df_ot[df_ot["proxy_type"] == "N"]
    ax.scatter(
        sub_n["ad_ot_score"], sub_n["spec_clipped"],
        c="#cccccc", s=6, alpha=0.35, linewidths=0,
        rasterized=True, zorder=1,
    )

    # Layer 2 — NMD / D / C (not junior-pass): coloured, no edge
    for pt in ["NMD", "D", "C"]:
        sub = df_ot[(df_ot["proxy_type"] == pt) & (~df_ot["junior_pass"])]
        ax.scatter(
            sub["ad_ot_score"], sub["spec_clipped"],
            c=PROXY_COLORS[pt], s=14, alpha=0.75, linewidths=0,
            rasterized=True, zorder=2,
        )

    # Layer 3 — junior pass: coloured, larger, black edge
    sub_t1 = df_ot[df_ot["junior_pass"]]
    ax.scatter(
        sub_t1["ad_ot_score"], sub_t1["spec_clipped"],
        c=[PROXY_COLORS[pt] for pt in sub_t1["proxy_type"]],
        s=38, alpha=1.0, linewidths=0.7, edgecolors="black",
        rasterized=True, zorder=3,
    )

    # Reference line: specificity = 1
    ax.axhline(1.0, color=ZERO_COLOR, lw=0.6, ls="--", alpha=0.35)

    def _annotate_group(rows: pd.DataFrame, color: str,
                        min_gap: float, dx_sign: int) -> None:
        pts_y    = rows["spec_clipped"].tolist()
        label_ys = push_apart_y(pts_y, min_gap=min_gap, ylim=(0.0, SPEC_CAP))
        for (_, r), ly in zip(rows.iterrows(), label_ys):
            px = r["ad_ot_score"]
            dx = dx_sign * 0.04
            ax.annotate(
                r["gene_name"],
                xy=(px, r["spec_clipped"]), xytext=(px + dx, ly),
                fontsize=7.5, color=color,
                arrowprops=dict(arrowstyle="-", color=color, lw=0.7),
                ha="right" if dx < 0 else "left",
                va="center", zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=color, linewidth=0.9, alpha=0.92),
            )

    # Top-4 OT score: labels go left (high x region)
    _annotate_group(top_ot,   ANNOT_OT,   min_gap=1.2, dx_sign=-1)
    # Top-4 specificity: labels go right (moderate x region)
    _annotate_group(top_spec, ANNOT_SPEC, min_gap=1.8, dx_sign=1)

    ax.set_xlim(-0.02, 0.60)
    ax.set_ylim(-0.8, SPEC_CAP * 1.05)
    ax.set_xlabel("AD OpenTargets score", fontsize=9)
    ax.set_ylabel("AD specificity  (AD / mean(PD, FTD, ALS))", fontsize=9)
    ax.set_title(
        "OT association vs disease specificity  |  colour = proxy type  |  black edge = junior pass",
        fontsize=10, pad=6,
    )
    ax.tick_params(labelsize=8)

    proxy_patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
                     for pt in ["C", "D", "NMD", "N"]]
    pass_patch    = mpatches.Patch(facecolor="white", edgecolor="black",
                                   label="Junior pass", linewidth=0.8)
    ot_patch      = mpatches.Patch(facecolor="white", edgecolor=ANNOT_OT,
                                   label="top-4 OT score", linewidth=0.9)
    spec_patch    = mpatches.Patch(facecolor="white", edgecolor=ANNOT_SPEC,
                                   label="top-4 AD specificity", linewidth=0.9)
    ax.legend(handles=proxy_patches + [pass_patch, ot_patch, spec_patch],
              bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)

    save(fig, "09_ot_specificity.png")


def plot_ot_label_by_celltype(df: pd.DataFrame) -> None:
    """10 — Stacked bar: OT label per cell type."""
    label_order = ["supported", "emerging", "novel"]
    cells       = [c for c in CELL_ORDER if c in df["cell_type"].unique()]

    counts = (
        df.groupby(["cell_type", "ad_ot_label"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=cells, columns=label_order, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    bottoms = np.zeros(len(cells))
    for lbl in label_order:
        vals = counts[lbl].values
        ax.bar(cells, vals, bottom=bottoms, color=OT_COLORS[lbl],
               label=lbl, edgecolor=BG, linewidth=0.5)
        bottoms += vals

    ax.set_xticklabels(cells, rotation=35, ha="right")
    ax.set_ylabel("Hits")
    ax.set_title("OpenTargets AD association label per cell type", fontweight="bold")
    ax.legend(title="OT label", bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "10_ot_label_by_celltype.png")


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 4 — MULTI-CELL-TYPE
# ═══════════════════════════════════════════════════════════════════════════

def plot_ct_count_dist(df: pd.DataFrame) -> None:
    """11 — Grouped bar: ct_count distribution by proxy_type."""
    proxy_order = ["C", "D", "NMD", "N"]
    ct_vals     = [1, 2, 3]

    counts = (
        df.groupby(["ct_count", "proxy_type"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=ct_vals, columns=proxy_order, fill_value=0)
    )

    x, w = np.arange(len(ct_vals)), 0.18
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, pt in enumerate(proxy_order):
        ax.bar(x + (i - 1.5) * w, counts[pt].values, width=w,
               color=PROXY_COLORS[pt], label=pt, edgecolor=BG)

    ax.set_xticks(x)
    ax.set_xticklabels(["1 cell type", "2 cell types", "3 cell types"])
    ax.set_ylabel("Transcripts")
    ax.set_title("Cell-type breadth distribution by proxy type", fontweight="bold")
    ax.legend(title="proxy_type")
    save(fig, "11_ct_count_distribution.png")


def plot_multict_heatmap(df: pd.DataFrame) -> None:
    """12 — Binary heatmap: multi-CT transcripts × cell types."""
    multi = df[df["multi_ct"]].copy()
    cells = [c for c in CELL_ORDER if c in multi["cell_type"].unique()]

    pivot = (
        multi.groupby(["transcript_name", "cell_type"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=cells, fill_value=0)
        .clip(upper=1)
    )
    pivot["_n"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("_n", ascending=False).drop(columns="_n")

    fig, ax = plt.subplots(figsize=(10, max(6, len(pivot) * 0.18)))
    sns.heatmap(
        pivot, ax=ax, cmap=NEON_CMAP, linewidths=0.3, linecolor=BG,
        cbar_kws={"label": "Hit present", "shrink": 0.4},
        vmin=0, vmax=1,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Transcript")
    ax.set_title(f"Multi-cell-type transcripts (n={len(pivot)}) × cell type", fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
    ax.tick_params(axis="y", labelsize=6)
    save(fig, "12_multict_heatmap.png")


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 5 — STRUCTURAL / DOMAIN
# ═══════════════════════════════════════════════════════════════════════════

def plot_cds_diff_hist(df: pd.DataFrame) -> None:
    """13 — Histogram: cds_diff_bp (full range + zoomed ≤1000 bp)."""
    cds = df[df["cds_diff_bp"] > 0].copy()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    bio_pairs = [("PC_CDS", BIO_COLORS["PC_CDS"]),
                 ("novel",  BIO_COLORS["novel"]),
                 ("PC_CDS_ND", BIO_COLORS["PC_CDS_ND"])]

    for ax, (title, subset) in zip(axes, [
        ("Full range", cds),
        ("≤1000 bp",  cds[cds["cds_diff_bp"] <= 1000]),
    ]):
        for bc, color in bio_pairs:
            sub = subset[subset["biotype_class"] == bc]
            if len(sub):
                ax.hist(sub["cds_diff_bp"], bins=40, color=color,
                        alpha=0.7, label=bc, edgecolor=BG)
        ax.set_xlabel("CDS symmetric diff (bp)")
        ax.set_ylabel("Hits")
        ax.set_title(title, fontweight="bold")
        ax.legend()

    fig.suptitle("CDS-altering hits: CDS symmetric difference", fontweight="bold", y=1.01)
    save(fig, "13_cds_diff_histogram.png")


def plot_domain_names_tier1(df: pd.DataFrame) -> None:
    """14 — Frequency bar: top domain names in junior-pass proxy-C hits."""
    t1c = df[
        df["junior_pass"] & (df["proxy_type"] == "C") &
        df["domain_names"].notna() & (df["domain_names"] != "")
    ]

    counter: Counter = Counter()
    for names in t1c["domain_names"]:
        for n in str(names).split(", "):
            n = n.strip()
            if n:
                counter[n] += 1

    if not counter:
        print("  [14] No domain names found in junior-pass C hits — skipping")
        return

    top = counter.most_common(20)
    names_, counts_ = zip(*top)

    fig, ax = plt.subplots(figsize=(8, 6))
    y = range(len(names_))
    ax.barh(list(y), list(counts_), color="#c0392b", edgecolor=BG)
    ax.set_yticks(list(y))
    ax.set_yticklabels(list(names_), fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Hits with this domain (junior pass, proxy C)")
    ax.set_title("Most common UniProt domains in junior-pass proxy-C hits", fontweight="bold")
    save(fig, "14_domain_names_tier1.png")


def plot_structural_feat_by_biotype(df: pd.DataFrame) -> None:
    """15 — Stacked bar: structural feature presence by biotype_class."""
    bio_order   = BIO_ORDER
    true_counts = (df.groupby("biotype_class")["has_structural_feat"]
                   .sum().reindex(bio_order, fill_value=0))
    false_counts = ((~df["has_structural_feat"])
                    .groupby(df["biotype_class"]).sum()
                    .reindex(bio_order, fill_value=0))

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(bio_order, true_counts.values,  color=_neon(0.50),
           label="Has structural feat", edgecolor=BG)
    ax.bar(bio_order, false_counts.values, color="#dddddd",
           label="No structural feat",  bottom=true_counts.values, edgecolor=BG)
    ax.set_ylabel("Hits")
    ax.set_title("UniProt structural feature presence by biotype class", fontweight="bold")
    ax.legend()
    save(fig, "15_structural_feat_by_biotype.png")


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 6 — JUNIOR-PASS DEEP-DIVE
# ═══════════════════════════════════════════════════════════════════════════

def plot_junior_pass_lollipop(df: pd.DataFrame) -> None:
    """16 — Lollipop: top 40 junior-pass genes by |Δ usage|."""
    t1 = (
        df[df["junior_pass"]]
        .sort_values("abs_delta", ascending=False)
        .drop_duplicates(subset="gene_name", keep="first")
        .head(40)
        .sort_values("abs_delta")
    )

    fig, ax = plt.subplots(figsize=(8, 12))
    y      = range(len(t1))
    colors = [PROXY_COLORS[pt] for pt in t1["proxy_type"]]

    ax.hlines(list(y), 0, t1["abs_delta"].values, color="#dddddd", lw=1)
    ax.scatter(t1["abs_delta"].values, list(y), c=colors, s=60, zorder=3)

    for i, (_, r) in enumerate(t1.iterrows()):
        ax.text(
            r["abs_delta"] + 0.006, i,
            f"OT={r['ad_ot_score']:.2f}",
            va="center", fontsize=7, color="#aaaaaa",
        )

    ax.set_yticks(list(y))
    ax.set_yticklabels(t1["gene_name"].values, fontsize=8)
    ax.set_xlabel("|Δ usage|")
    ax.set_xlim(0, 1.08)
    ax.set_title("Top 40 junior-pass genes by |Δ usage|", fontweight="bold")

    patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
               for pt in ["C", "D", "NMD", "N"]]
    ax.legend(handles=patches, title="proxy_type", loc="lower right")
    save(fig, "16_tier1_lollipop.png")


def plot_junior_pass_feature_heatmap(df: pd.DataFrame) -> None:
    """17 — Binary feature matrix: top 30 junior-pass genes × features."""
    t1 = (
        df[df["junior_pass"]]
        .sort_values("abs_delta", ascending=False)
        .drop_duplicates(subset="gene_name", keep="first")
        .head(30)
    )

    features = {
        "proxy_C":      (t1["proxy_type"] == "C").astype(int).values,
        "proxy_D":      (t1["proxy_type"] == "D").astype(int).values,
        "has_domain":   t1["has_structural_feat"].astype(int).values,
        "multi_ct":     t1["multi_ct"].astype(int).values,
        "ad_prior":     t1["ad_prior_flag"].astype(int).values,
        "OT_supported": (t1["ad_ot_label"] == "supported").astype(int).values,
        "OT_emerging":  (t1["ad_ot_label"] == "emerging").astype(int).values,
        "OT_genetic>0": (t1["ad_ot_genetic"] > 0).astype(int).values,
        "OT_lit>0":     (t1["ad_ot_lit"] > 0).astype(int).values,
    }

    mat = pd.DataFrame(features, index=t1["gene_name"].values)

    fig, ax = plt.subplots(figsize=(9, 9))
    sns.heatmap(
        mat, ax=ax, cmap=NEON_CMAP, linewidths=0.5, linecolor=BG,
        cbar_kws={"label": "Feature present", "shrink": 0.4},
        vmin=0, vmax=1,
    )
    ax.set_title("Junior-pass top-30 genes: feature matrix", fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
    save(fig, "17_tier1_feature_heatmap.png")


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 7 — THREE NEW EVIDENCE PLOTS
# ═══════════════════════════════════════════════════════════════════════════

def plot_log2fc_vs_delta(df: pd.DataFrame) -> None:
    """20 — log2FC vs delta_usage, expression_vs_psi style, colour = proxy_type."""
    YLIM = (-1.05, 1.05)
    sub  = df.dropna(subset=["log2FC"]).copy()

    fig, ax = plt.subplots(figsize=(9, 7))

    # Layer 1 — N: grey background
    sn = sub[sub["proxy_type"] == "N"]
    ax.scatter(sn["log2FC"], sn["delta_usage"],
               c="#cccccc", s=6, alpha=0.35, linewidths=0, rasterized=True, zorder=1)

    # Layer 2 — NMD / D / C (not junior-pass)
    for pt in ["NMD", "D", "C"]:
        s = sub[(sub["proxy_type"] == pt) & (~sub["junior_pass"])]
        ax.scatter(s["log2FC"], s["delta_usage"],
                   c=PROXY_COLORS[pt], s=14, alpha=0.75, linewidths=0,
                   rasterized=True, zorder=2)

    # Layer 3 — junior pass: black edge
    t1 = sub[sub["junior_pass"]]
    ax.scatter(t1["log2FC"], t1["delta_usage"],
               c=[PROXY_COLORS[pt] for pt in t1["proxy_type"]],
               s=38, alpha=1.0, linewidths=0.7, edgecolors="black",
               rasterized=True, zorder=3)

    ax.axhline(0,     color=ZERO_COLOR,   lw=0.6, ls="--", alpha=0.35)
    ax.axvline(0,     color=ZERO_COLOR,   lw=0.6, ls="--", alpha=0.35)
    ax.axhline( 0.25, color=THRESH_COLOR, lw=0.8, ls="--", alpha=0.7)
    ax.axhline(-0.25, color=THRESH_COLOR, lw=0.8, ls="--", alpha=0.7)

    # Annotate top-6 junior-pass C hits by |Δ usage|
    top = (sub[sub["junior_pass"] & (sub["proxy_type"] == "C")]
           .drop_duplicates("gene_name")
           .nlargest(6, "abs_delta"))
    if len(top):
        pts_y    = top["delta_usage"].tolist()
        label_ys = push_apart_y(pts_y, min_gap=0.13, ylim=YLIM)
        for (_, r), ly in zip(top.iterrows(), label_ys):
            px = r["log2FC"]
            dx = 0.5 if px >= 0 else -0.5
            ax.annotate(
                r["gene_name"],
                xy=(px, r["delta_usage"]), xytext=(px + dx, ly),
                fontsize=7.5, color=ANNOT_COLOR,
                arrowprops=dict(arrowstyle="-", color=ANNOT_COLOR, lw=0.7),
                ha="left" if dx > 0 else "right",
                va="center", zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=ANNOT_COLOR, linewidth=0.7, alpha=0.9),
            )

    ax.set_ylim(*YLIM)
    ax.set_xlabel("log₂FC (AD / Control)", fontsize=9)
    ax.set_ylabel("Δ usage (AD − Control)", fontsize=9)
    ax.set_title(
        "Expression change vs isoform usage change  |  colour = proxy type  |  black edge = junior pass",
        fontsize=10, pad=6,
    )
    ax.tick_params(labelsize=8)

    proxy_patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
                     for pt in ["C", "D", "NMD", "N"]]
    pass_patch = mpatches.Patch(facecolor="white", edgecolor="black",
                                 label="Junior pass", linewidth=0.8)
    ax.legend(handles=proxy_patches + [pass_patch],
              title="proxy_type", bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)
    save(fig, "20_log2fc_vs_delta_usage.png")


def plot_ot_genetic_vs_lit(df: pd.DataFrame) -> None:
    """21 — OT genetic association vs literature mining score."""
    sub = df[df["ad_ot_score"] > 0].copy()
    lim = max(sub["ad_ot_genetic"].max(), sub["ad_ot_lit"].max()) * 1.06

    fig, ax = plt.subplots(figsize=(8, 7))

    # Layer 1 — N
    sn = sub[sub["proxy_type"] == "N"]
    ax.scatter(sn["ad_ot_genetic"], sn["ad_ot_lit"],
               c="#cccccc", s=6, alpha=0.35, linewidths=0, rasterized=True, zorder=1)

    # Layer 2 — NMD / D / C (not junior-pass)
    for pt in ["NMD", "D", "C"]:
        s = sub[(sub["proxy_type"] == pt) & (~sub["junior_pass"])]
        ax.scatter(s["ad_ot_genetic"], s["ad_ot_lit"],
                   c=PROXY_COLORS[pt], s=14, alpha=0.75, linewidths=0,
                   rasterized=True, zorder=2)

    # Layer 3 — junior pass: black edge
    t1 = sub[sub["junior_pass"]]
    ax.scatter(t1["ad_ot_genetic"], t1["ad_ot_lit"],
               c=[PROXY_COLORS[pt] for pt in t1["proxy_type"]],
               s=38, alpha=1.0, linewidths=0.7, edgecolors="black",
               rasterized=True, zorder=3)

    # Diagonal: genetic = lit
    ax.plot([0, lim], [0, lim], "--", color=ZERO_COLOR, lw=0.6, alpha=0.4)
    ax.text(lim * 0.97, lim * 0.03, "genetic dominant →",
            ha="right", va="bottom", fontsize=8, color="#aaaaaa", style="italic")
    ax.text(lim * 0.03, lim * 0.97, "← literature dominant",
            ha="left",  va="top",    fontsize=8, color="#aaaaaa", style="italic")

    # Annotate top junior-pass genes with OT evidence
    top = (df[df["junior_pass"] & (df["ad_ot_score"] > 0.05)]
           .drop_duplicates("gene_name")
           .nlargest(6, "abs_delta"))
    seen: set[str] = set()
    for _, r in top.iterrows():
        if r["gene_name"] in seen:
            continue
        px, py = r["ad_ot_genetic"], r["ad_ot_lit"]
        dx = 0.03 if px < lim * 0.55 else -0.03
        ax.annotate(
            r["gene_name"],
            xy=(px, py), xytext=(px + dx, py),
            fontsize=7.5, color=ANNOT_COLOR,
            arrowprops=dict(arrowstyle="-", color=ANNOT_COLOR, lw=0.7),
            ha="left" if dx > 0 else "right",
            va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor=ANNOT_COLOR, linewidth=0.7, alpha=0.9),
        )
        seen.add(r["gene_name"])

    ax.set_xlim(-0.02, lim)
    ax.set_ylim(-0.02, lim)
    ax.set_xlabel("OT genetic association score  (GWAS + rare variants)", fontsize=9)
    ax.set_ylabel("OT literature mining score", fontsize=9)
    ax.set_title(
        "AD evidence type decomposition  |  colour = proxy type  |  black edge = junior pass",
        fontsize=10, pad=6,
    )
    ax.tick_params(labelsize=8)

    proxy_patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
                     for pt in ["C", "D", "NMD", "N"]]
    pass_patch = mpatches.Patch(facecolor="white", edgecolor="black",
                                 label="Junior pass", linewidth=0.8)
    ax.legend(handles=proxy_patches + [pass_patch],
              title="proxy_type", bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)
    save(fig, "21_ot_genetic_vs_lit.png")


def plot_cpm_scatter(df: pd.DataFrame) -> None:
    """22 — mean_CPM_AD vs mean_CPM_CT (log10+1), colour = proxy_type, junior-pass bordered."""
    sub = df.copy()
    sub["log_cpm_ct"] = np.log10(sub["mean_CPM_CT"] + 1)
    sub["log_cpm_ad"] = np.log10(sub["mean_CPM_AD"] + 1)
    lim = max(sub["log_cpm_ct"].max(), sub["log_cpm_ad"].max()) * 1.05

    fig, ax = plt.subplots(figsize=(8, 7))

    # Layer 1 — N
    sn = sub[sub["proxy_type"] == "N"]
    ax.scatter(sn["log_cpm_ct"], sn["log_cpm_ad"],
               c="#cccccc", s=6, alpha=0.30, linewidths=0, rasterized=True, zorder=1)

    # Layer 2 — NMD / D / C (not junior-pass)
    for pt in ["NMD", "D", "C"]:
        s = sub[(sub["proxy_type"] == pt) & (~sub["junior_pass"])]
        ax.scatter(s["log_cpm_ct"], s["log_cpm_ad"],
                   c=PROXY_COLORS[pt], s=14, alpha=0.75, linewidths=0,
                   rasterized=True, zorder=2)

    # Layer 3 — junior pass: black edge
    t1 = sub[sub["junior_pass"]]
    ax.scatter(t1["log_cpm_ct"], t1["log_cpm_ad"],
               c=[PROXY_COLORS[pt] for pt in t1["proxy_type"]],
               s=38, alpha=1.0, linewidths=0.7, edgecolors="black",
               rasterized=True, zorder=3)

    # Diagonal: equal expression
    ax.plot([0, lim], [0, lim], "--", color=ZERO_COLOR, lw=0.6, alpha=0.4)

    # Annotate top junior-pass C hits by |Δ usage|
    top = (sub[sub["junior_pass"] & (sub["proxy_type"] == "C")]
           .drop_duplicates("gene_name")
           .nlargest(6, "abs_delta"))
    if len(top):
        pts_y    = top["log_cpm_ad"].tolist()
        label_ys = push_apart_y(pts_y, min_gap=0.15, ylim=(0, lim))
        for (_, r), ly in zip(top.iterrows(), label_ys):
            px = r["log_cpm_ct"]
            dx = 0.15 if px < lim * 0.6 else -0.15
            ax.annotate(
                r["gene_name"],
                xy=(px, r["log_cpm_ad"]), xytext=(px + dx, ly),
                fontsize=7.5, color=ANNOT_COLOR,
                arrowprops=dict(arrowstyle="-", color=ANNOT_COLOR, lw=0.7),
                ha="left" if dx > 0 else "right",
                va="center", zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor=ANNOT_COLOR, linewidth=0.7, alpha=0.9),
            )

    # Tick labels: convert log10 back to CPM
    ticks = [0, 1, 2, 3, 4]
    tick_labels = ["0", "10", "100", "1k", "10k"]
    ax.set_xticks(ticks); ax.set_xticklabels(tick_labels)
    ax.set_yticks(ticks); ax.set_yticklabels(tick_labels)
    ax.set_xlim(-0.1, lim)
    ax.set_ylim(-0.1, lim)
    ax.set_xlabel("mean CPM — Control", fontsize=9)
    ax.set_ylabel("mean CPM — AD", fontsize=9)
    ax.set_title(
        "Transcript expression level: AD vs Control  |  colour = proxy type  |  black edge = junior pass",
        fontsize=10, pad=6,
    )
    ax.tick_params(labelsize=8)

    proxy_patches = [mpatches.Patch(color=PROXY_COLORS[pt], label=pt)
                     for pt in ["C", "D", "NMD", "N"]]
    pass_patch = mpatches.Patch(facecolor="white", edgecolor="black",
                                 label="Junior pass", linewidth=0.8)
    ax.legend(handles=proxy_patches + [pass_patch],
              title="proxy_type", bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)
    save(fig, "22_cpm_ad_vs_ct.png")


# ═══════════════════════════════════════════════════════════════════════════
# SANKEY
# ═══════════════════════════════════════════════════════════════════════════

def plot_sankey(df: pd.DataFrame) -> None:
    """19 — Sankey: candidate group → biotype_class → proxy_type → junior gate."""
    import plotly.graph_objects as go

    def _rgba(hex_color: str, alpha: float = 0.45) -> str:
        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)
        return f"rgba({r},{g},{b},{alpha})"

    group_order = ["trial_failure_candidate", "new_target_candidate"]
    bio_order   = BIO_ORDER
    proxy_order = ["C", "D", "NMD", "N"]
    gate_order  = [True, False]
    gate_labels = {True: "Pass", False: "Drop"}

    n_grp = len(group_order)
    n_bio = len(bio_order)
    n_pxy = len(proxy_order)

    group_idx = {g: i                        for i, g in enumerate(group_order)}
    bio_idx   = {b: n_grp + i                for i, b in enumerate(bio_order)}
    proxy_idx = {p: n_grp + n_bio + i        for i, p in enumerate(proxy_order)}
    gate_idx  = {t: n_grp + n_bio + n_pxy + i for i, t in enumerate(gate_order)}

    group_labels = {
        "trial_failure_candidate": "Trial-failure candidate",
        "new_target_candidate":    "New-target candidate",
    }
    node_labels = (
        [group_labels[g] for g in group_order]
        + bio_order
        + proxy_order
        + [gate_labels[t] for t in gate_order]
    )
    node_colors = (
        [GROUP_COLORS[g]      for g in group_order]
        + [BIO_COLORS[b]        for b in bio_order]
        + [PROXY_COLORS[p]    for p in proxy_order]
        + [GATE_COLORS[t]     for t in gate_order]
    )

    sources, targets, values, link_colors = [], [], [], []

    # candidate group → biotype
    for grp in group_order:
        for bio in bio_order:
            n = int(((df["candidate_group"] == grp) & (df["biotype_class"] == bio)).sum())
            if n:
                sources.append(group_idx[grp])
                targets.append(bio_idx[bio])
                values.append(n)
                link_colors.append(_rgba(GROUP_COLORS[grp], 0.35))

    # biotype → proxy
    for bio in bio_order:
        for pxy in proxy_order:
            n = int(((df["biotype_class"] == bio) & (df["proxy_type"] == pxy)).sum())
            if n:
                sources.append(bio_idx[bio])
                targets.append(proxy_idx[pxy])
                values.append(n)
                link_colors.append(_rgba(BIO_COLORS[bio], 0.35))

    # proxy → junior gate
    for pxy in proxy_order:
        for t in gate_order:
            n = int(((df["proxy_type"] == pxy) & (df["junior_pass"] == t)).sum())
            if n:
                sources.append(proxy_idx[pxy])
                targets.append(gate_idx[t])
                values.append(n)
                link_colors.append(_rgba(PROXY_COLORS[pxy], 0.35))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=18,
            thickness=22,
            line=dict(color="#cccccc", width=0.5),
            label=node_labels,
            color=node_colors,
            hovertemplate="%{label}: %{value}<extra></extra>",
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
        ),
    ))

    fig.update_layout(
        title=dict(
            text="Candidate group → Biotype class → Proxy type → Junior gate",
            font=dict(size=15, family="sans-serif", color=FG),
            x=0.5, xanchor="center",
        ),
        font=dict(family="sans-serif", size=12, color=FG),
        paper_bgcolor=BG,
        plot_bgcolor=BG,
        width=1100,
        height=600,
        margin=dict(l=20, r=20, t=60, b=20),
    )

    out = OUT_DIR / "19_sankey_group_biotype_proxy_gate.png"
    fig.write_image(str(out), scale=2)
    print(f"  saved → {out.name}")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    print("Loading hits_enriched.csv ...")
    df = load()
    print(f"  {len(df):,} hits  |  {df['gene_name'].nunique():,} genes\n")

    print("Rendering plots ...")
    plot_junior_gate_overview(df)
    plot_biotype_by_celltype(df)
    plot_proxy_by_celltype(df)
    plot_usage_direction(df)
    plot_volcano(df)
    plot_psi_scatter(df)
    plot_delta_violin(df)
    plot_ot_vs_effect(df)
    plot_ot_specificity(df)
    plot_ot_label_by_celltype(df)
    plot_ct_count_dist(df)
    plot_multict_heatmap(df)
    plot_cds_diff_hist(df)
    plot_domain_names_tier1(df)
    plot_structural_feat_by_biotype(df)
    plot_junior_pass_lollipop(df)
    plot_junior_pass_feature_heatmap(df)
    plot_log2fc_vs_delta(df)
    plot_ot_genetic_vs_lit(df)
    plot_cpm_scatter(df)
    plot_sankey(df)

    n_saved = len(list(OUT_DIR.glob("*.png")))
    print(f"\nDone — {n_saved} plots in {OUT_DIR}")


if __name__ == "__main__":
    main()

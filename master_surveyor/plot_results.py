"""MASTER_SURVEYOR — preliminary visualization suite for the 228-hit J4 shortlist."""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

from master_surveyor.config import HITS_CSV, OUT_DIR

# ---------------------------------------------------------------------------
# Design system (mirrors assistant_surveyor / junior_surveyor conventions)
# ---------------------------------------------------------------------------
NEON_CMAP = LinearSegmentedColormap.from_list(
    "neon", ["#0a0020", "#5522cc", "#cc1188", "#cc5500", "#bbbb00"]
)

def _neon(t: float) -> str:
    return mcolors.to_hex(NEON_CMAP(t))

CHANGE_COLORS = {
    "frameshift_stop":    _neon(0.95),
    "N_truncation":       "#E5670A",
    "C_truncation":       "#C41B3A",
    "internal_indel":     _neon(0.50),
    "substitution":       _neon(0.42),
    "N_extension":        _neon(0.35),
    "internal_insertion": _neon(0.28),
    "C_extension":        _neon(0.20),
    "identical":          "#cccccc",
    "no_sequence":        "#eeeeee",
}
_CHANGE_ORDER = [
    "frameshift_stop", "N_truncation", "C_truncation",
    "internal_indel", "substitution",
    "N_extension", "internal_insertion", "C_extension",
    "identical", "no_sequence",
]

OT_LABEL_COLORS = {"supported": _neon(0.50), "emerging": _neon(0.75), "novel": "#999999"}

BG         = "#ffffff"
FG         = "#1a1a1a"
GRID       = "#e0e0e0"
ANNOT_NAVY = "#003399"
ANNOT_TEAL = "#007755"
LOSS_COLOR = "#C41B3A"

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _push_apart(ys: list[float], min_gap: float, lo: float, hi: float) -> list[float]:
    """Iteratively separate label positions to avoid overlap."""
    ys = list(ys)
    for _ in range(400):
        moved = False
        order = sorted(range(len(ys)), key=lambda i: ys[i])
        for k in range(len(order) - 1):
            i, j = order[k], order[k + 1]
            gap = ys[j] - ys[i]
            if gap < min_gap:
                push = (min_gap - gap) / 2
                ys[i] -= push
                ys[j] += push
                moved = True
        if not moved:
            break
    return [max(lo, min(hi, y)) for y in ys]


def _style(ax: plt.Axes) -> None:
    ax.set_facecolor(BG)
    ax.tick_params(colors=FG, labelsize=8)
    ax.xaxis.label.set_color(FG)
    ax.yaxis.label.set_color(FG)
    ax.title.set_color(FG)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID)
    ax.set_axisbelow(True)


def _grid(ax: plt.Axes) -> None:
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)


def save(fig: plt.Figure, name: str, extra_artists: list | None = None) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.patch.set_facecolor(BG)
    # bbox_inches="tight" doesn't reliably measure a legend added via
    # ax.add_artist() (needed when a plot has two legends) — pass it
    # explicitly via bbox_extra_artists or it gets clipped off the canvas.
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=BG,
                bbox_extra_artists=extra_artists)
    plt.close(fig)
    print(f"  saved → {path.name}")


def _split_domains(s) -> list[str]:
    if not isinstance(s, str) or s.strip().lower() in ("", "none", "nan"):
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def _size_scale(values: pd.Series, lo: float = 25.0, hi: float = 420.0) -> np.ndarray:
    """Map a raw value column to matplotlib scatter marker *area* (points^2)."""
    v = values.astype(float).values
    if v.max() == v.min():
        return np.full(len(v), (lo + hi) / 2)
    return lo + (v - v.min()) / (v.max() - v.min()) * (hi - lo)


def _size_legend(ax: plt.Axes, raw: pd.Series, sizes: np.ndarray, title: str,
                  loc: str, fmt: str = "{:.2f}", fig: plt.Figure | None = None,
                  bbox_to_anchor: tuple[float, float] | None = None) -> mpatches.Patch:
    """Manual bubble-size legend (low/mid/high) — matplotlib has no automatic one.

    Pass `fig` + `bbox_to_anchor` to place it outside the axes (figure-fraction
    coordinates) instead of floating inside the plot at `loc`.
    """
    import matplotlib.lines as mlines
    order = np.argsort(raw.values)
    picks = [order[0], order[len(order) // 2], order[-1]]
    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", markerfacecolor="none",
                      markeredgecolor=FG, markeredgewidth=0.9,
                      markersize=np.sqrt(sizes[i]), label=fmt.format(raw.values[i]))
        for i in picks
    ]
    target = fig if fig is not None else ax
    kwargs = dict(handles=handles, title=title, loc=loc, fontsize=7,
                  title_fontsize=7.5, frameon=True, framealpha=1, edgecolor=GRID,
                  labelspacing=1.4, borderpad=1.1)
    if bbox_to_anchor is not None:
        kwargs["bbox_to_anchor"] = bbox_to_anchor
    return target.legend(**kwargs)


def load() -> pd.DataFrame:
    """Load junior_surveyor's hits_deep.csv, restricted to selected_for_next_stage."""
    df = pd.read_csv(HITS_CSV)
    df = df[df["selected_for_next_stage"]].copy()
    df["protein_change_type"] = df["protein_change_type"].fillna("no_sequence")
    df["affected_domain"]     = df["affected_domain"].fillna("none")
    df["chembl_max_phase"]    = df["chembl_max_phase"].fillna(0).astype(int)
    df["dgidb_interactions"]  = df["dgidb_interactions"].fillna(0).astype(int)
    return df


# ---------------------------------------------------------------------------
# Plot 01 — Volcano: significance vs effect size
# ---------------------------------------------------------------------------

def plot_volcano(df: pd.DataFrame) -> None:
    """01 — -log10(chi_padj) vs delta_usage, colour = protein_change_type."""
    d = df.copy()
    d["neg_log10_padj"] = -np.log10(d["chi_padj"].clip(lower=1e-300))
    ymax = d["neg_log10_padj"].max() * 1.10

    fig, ax = plt.subplots(figsize=(9, 7))

    present = [ct for ct in _CHANGE_ORDER if ct in d["protein_change_type"].unique()]
    for ct in present:
        sub = d[d["protein_change_type"] == ct]
        ax.scatter(sub["delta_usage"], sub["neg_log10_padj"],
                   c=CHANGE_COLORS[ct], s=32, alpha=0.85, linewidths=0.4,
                   edgecolors="white", zorder=2, rasterized=True)

    ax.axhline(0, color="black", lw=0.6, ls="--", alpha=0.3)
    ax.axvline(0, color="black", lw=0.6, ls="--", alpha=0.3)
    ax.axhline(-np.log10(0.05), color="red", lw=0.8, ls="--", alpha=0.6)

    # Label the most extreme points (rank on significance + |effect|)
    d["_rank"] = (d["neg_log10_padj"].rank(pct=True)
                  + d["delta_usage"].abs().rank(pct=True))
    top = d.nlargest(15, "_rank")
    label_ys = _push_apart(top["neg_log10_padj"].tolist(), min_gap=ymax * 0.035,
                           lo=0.0, hi=ymax)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        dx = 0.05 if row["delta_usage"] >= 0 else -0.05
        ax.annotate(
            row["gene_name"],
            xy=(row["delta_usage"], row["neg_log10_padj"]),
            xytext=(row["delta_usage"] + dx, ly),
            fontsize=7.5, color=ANNOT_NAVY,
            arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
            ha="left" if dx > 0 else "right", va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9),
        )

    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Δ usage (AD − Control)", fontsize=9)
    ax.set_ylabel("−log₁₀ (χ² adjusted p-value)", fontsize=9)
    ax.set_title(
        f"Master shortlist (n={len(d):,}): significance vs effect size  |  colour = protein change type",
        fontsize=10, pad=8,
    )
    ax.tick_params(labelsize=8)
    _style(ax); _grid(ax)

    handles = [mpatches.Patch(color=CHANGE_COLORS[ct], label=ct) for ct in present]
    ax.legend(handles=handles, title="protein_change_type",
              bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, title_fontsize=8)

    save(fig, "01_volcano.png")


# ---------------------------------------------------------------------------
# Plot 02 — Protein change type distribution
# ---------------------------------------------------------------------------

def plot_change_type_bar(df: pd.DataFrame) -> None:
    """02 — Distribution of protein_change_type (deduped per transcript)."""
    d = df.drop_duplicates("transcript_name")
    counts = d["protein_change_type"].value_counts()
    order = [ct for ct in _CHANGE_ORDER if ct in counts.index]
    counts = counts.reindex(order)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ypos = np.arange(len(order))
    ax.barh(ypos, counts.values, color=[CHANGE_COLORS[ct] for ct in order],
            height=0.6, zorder=2)
    for y, v in zip(ypos, counts.values):
        ax.text(v + max(counts.values) * 0.015, y, f"{v}  ({v / len(d) * 100:.0f}%)",
                va="center", fontsize=8, color=FG)

    ax.set_yticks(ypos)
    ax.set_yticklabels(order, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, counts.max() * 1.22)
    ax.set_xlabel("transcripts", fontsize=9)
    ax.set_title(f"Protein change type — master shortlist (n={len(d):,} transcripts)",
                 fontsize=10, pad=8)
    _style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(False)

    save(fig, "02_change_type_bar.png")


# ---------------------------------------------------------------------------
# Plot 03 — Domain gain/loss per transcript
# ---------------------------------------------------------------------------

def plot_domain_gain_loss(df: pd.DataFrame) -> None:
    """03 — Diverging bar: domains lost vs gained, one row per transcript."""
    d_all = df.drop_duplicates("transcript_name").copy()
    d_all["n_lost"]   = d_all["domains_lost"].apply(lambda s: len(_split_domains(s)))
    d_all["n_gained"] = d_all["domains_gained"].apply(lambda s: len(_split_domains(s)))
    d_all = d_all[(d_all["n_lost"] > 0) | (d_all["n_gained"] > 0)].copy()
    n_affected = len(d_all)

    d_all["total"] = d_all["n_lost"] + d_all["n_gained"]
    d = d_all.sort_values("total", ascending=False).head(40)
    d = d.iloc[::-1]  # largest at top

    ypos = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(8, max(4.5, len(d) * 0.24)))
    ax.barh(ypos, -d["n_lost"], color=LOSS_COLOR, height=0.65, zorder=2, label="domain(s) lost")
    ax.barh(ypos, d["n_gained"], color=ANNOT_TEAL, height=0.65, zorder=2, label="domain(s) gained")
    ax.axvline(0, color=FG, lw=0.8)

    ax.set_yticks(ypos)
    ax.set_yticklabels(d["transcript_name"], fontsize=7.5)
    loss_max, gain_max = d["n_lost"].max(), max(d["n_gained"].max(), 1)
    # Loss dominates heavily here — give gain just enough room to be legible
    # rather than mirroring loss's scale and wasting half the plot on blank space.
    ax.set_xlim(-loss_max * 1.15, max(gain_max * 3.5, loss_max * 0.15))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: str(int(abs(x)))))
    ax.set_xlabel("← domains lost                domains gained →", fontsize=8)
    ax.set_title(
        f"Domain gain/loss per transcript  |  top {len(d)} by total events "
        f"({n_affected} affected of {df['transcript_name'].nunique()})",
        fontsize=9.5, pad=8,
    )
    _style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(False)
    ax.legend(loc="lower right", fontsize=8, frameon=True, framealpha=1, edgecolor=GRID)

    save(fig, "03_domain_gain_loss.png")


# ---------------------------------------------------------------------------
# Plot 04 — Drug / target evidence
# ---------------------------------------------------------------------------

def plot_drug_evidence(df: pd.DataFrame) -> None:
    """04 — ChEMBL max phase vs DGIdb interactions, one point per gene."""
    d = df.drop_duplicates("gene_name").copy()
    rng = np.random.default_rng(42)
    d["_x"] = d["chembl_max_phase"].astype(float) + rng.uniform(-0.12, 0.12, size=len(d))

    fig, ax = plt.subplots(figsize=(8, 6.5))
    present = [ct for ct in _CHANGE_ORDER if ct in d["protein_change_type"].unique()]
    for ct in present:
        sub = d[d["protein_change_type"] == ct]
        ax.scatter(sub["_x"], sub["dgidb_interactions"], c=CHANGE_COLORS[ct], s=34,
                   alpha=0.85, linewidths=0.4, edgecolors="white", zorder=2, rasterized=True)

    ymax = d["dgidb_interactions"].max() * 1.12
    top = d.nlargest(10, "dgidb_interactions")
    label_ys = _push_apart(top["dgidb_interactions"].tolist(), min_gap=ymax * 0.045,
                           lo=0.0, hi=ymax)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        ax.annotate(
            row["gene_name"],
            xy=(row["_x"], row["dgidb_interactions"]),
            xytext=(row["_x"] + 0.35, ly),
            fontsize=7.5, color=ANNOT_NAVY,
            arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
            ha="left", va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9),
        )

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlim(-0.4, 5.4)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("ChEMBL max clinical phase  (0 = no compound on record)", fontsize=9)
    ax.set_ylabel("DGIdb interaction count", fontsize=9)
    ax.set_title(f"Drug / target evidence — master shortlist (n={len(d):,} genes)",
                 fontsize=10, pad=8)
    _style(ax); _grid(ax)

    handles = [mpatches.Patch(color=CHANGE_COLORS[ct], label=ct) for ct in present]
    ax.legend(handles=handles, title="protein_change_type",
              bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7.5, title_fontsize=8)

    save(fig, "04_drug_evidence_scatter.png")


# ---------------------------------------------------------------------------
# Plot 05 — AD relevance
# ---------------------------------------------------------------------------

def plot_ad_specificity_vs_score(df: pd.DataFrame) -> None:
    """05 — ad_specificity vs ad_ot_score, colour = ad_ot_label, one point per gene."""
    d = df.drop_duplicates("gene_name").copy()

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for lbl in ["novel", "emerging", "supported"]:
        sub = d[d["ad_ot_label"] == lbl]
        if sub.empty:
            continue
        ax.scatter(sub["ad_ot_score"], sub["ad_specificity"],
                   c=OT_LABEL_COLORS[lbl], s=34, alpha=0.85, linewidths=0.4,
                   edgecolors="white", zorder=2, rasterized=True, label=lbl)

    xmax = d["ad_ot_score"].max() * 1.08
    ymax = d["ad_specificity"].max() * 1.08
    top = d.nlargest(10, "ad_ot_score")
    label_ys = _push_apart(top["ad_specificity"].tolist(), min_gap=ymax * 0.05,
                           lo=0.0, hi=ymax)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        ax.annotate(
            row["gene_name"],
            xy=(row["ad_ot_score"], row["ad_specificity"]),
            xytext=(row["ad_ot_score"] + xmax * 0.03, ly),
            fontsize=7.5, color=ANNOT_NAVY,
            arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
            ha="left", va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9),
        )
    ax.set_xlim(0, xmax * 1.12)

    ax.set_xlabel("AD OpenTargets association score", fontsize=9)
    ax.set_ylabel("AD specificity  (vs. other neurodegenerative disease)", fontsize=9)
    ax.set_title(f"Disease relevance — master shortlist (n={len(d):,} genes)",
                 fontsize=10, pad=8)
    _style(ax); _grid(ax)
    ax.legend(title="ad_ot_label", bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, title_fontsize=8)

    save(fig, "05_ad_specificity_vs_score.png")


# ---------------------------------------------------------------------------
# Plot 06 — Candidate scorecard heatmap
# ---------------------------------------------------------------------------

def plot_scorecard_heatmap(df: pd.DataFrame, top_n: int = 40) -> None:
    """06 — Per-gene normalized heatmap across the key selection metrics."""
    gene_level = df.drop_duplicates("gene_name").set_index("gene_name")

    agg = df.groupby("gene_name").agg(
        max_abs_delta=("delta_usage", lambda s: s.abs().max()),
        min_padj=("chi_padj", "min"),
    )

    mat = pd.DataFrame({
        "|Δ PSI|":          agg["max_abs_delta"],
        "-log10(padj)":     -np.log10(agg["min_padj"].clip(lower=1e-300)),
        "Drug phase":       gene_level["chembl_max_phase"],
        "DGIdb (log1p)":    np.log1p(gene_level["dgidb_interactions"]),
        "AD OT score":      gene_level["ad_ot_score"],
        "AD specificity":   gene_level["ad_specificity"],
    })

    norm = (mat - mat.min()) / (mat.max() - mat.min()).replace(0, 1)
    order = norm.sum(axis=1).sort_values(ascending=False).head(top_n).index
    norm, raw = norm.loc[order], mat.loc[order]

    fig, ax = plt.subplots(figsize=(7.5, max(6.0, len(order) * 0.24)))
    im = ax.imshow(norm.values, aspect="auto", cmap=NEON_CMAP, vmin=0, vmax=1)

    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7.5)

    # NEON_CMAP runs near-black (low norm) -> magenta -> orange -> yellow
    # (high norm), so dark text only reads on the bright high-norm end.
    for i in range(len(order)):
        for j in range(len(mat.columns)):
            v = raw.iloc[i, j]
            txt = f"{v:.2f}" if abs(v) < 10 else f"{v:.0f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                    color=FG if norm.iloc[i, j] > 0.6 else "white")

    cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("normalized per column", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    ax.set_title(f"Master shortlist scorecard — top {len(order)} of {len(gene_level):,} genes",
                 fontsize=10, pad=8)
    _style(ax)
    for spine in ax.spines.values():
        spine.set_visible(False)

    save(fig, "06_scorecard_heatmap.png")


# ---------------------------------------------------------------------------
# 4-channel scatter helpers: one row per gene / per transcript
# ---------------------------------------------------------------------------

def _rep_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    """One row per gene — the transcript×cell_type with the largest |delta_usage|."""
    idx = df.groupby("gene_name")["delta_usage"].apply(lambda s: s.abs().idxmax())
    return df.loc[idx].copy()


def _rep_per_transcript(df: pd.DataFrame) -> pd.DataFrame:
    """One row per transcript — its most significant cell_type (min chi_padj)."""
    idx = df.groupby("transcript_name")["chi_padj"].idxmin()
    return df.loc[idx].copy()


# ---------------------------------------------------------------------------
# Plot 07 — Fused priority: disease relevance x druggability x structure x effect
# ---------------------------------------------------------------------------

def plot_priority_bubble(df: pd.DataFrame) -> None:
    """07 — ad_ot_score vs chembl_max_phase, colour = change type, size = |delta_usage|."""
    d = _rep_per_gene(df)
    rng = np.random.default_rng(7)
    d = d.assign(_y=d["chembl_max_phase"].astype(float) + rng.uniform(-0.12, 0.12, len(d)))
    sizes = _size_scale(d["delta_usage"].abs())

    fig, ax = plt.subplots(figsize=(8.5, 7))
    present = [ct for ct in _CHANGE_ORDER if ct in d["protein_change_type"].unique()]
    for ct in present:
        m = d["protein_change_type"] == ct
        ax.scatter(d.loc[m, "ad_ot_score"], d.loc[m, "_y"], s=sizes[m.values],
                   c=CHANGE_COLORS[ct], alpha=0.8, linewidths=0.4,
                   edgecolors="white", zorder=2, rasterized=True)

    d["_rank"] = d["ad_ot_score"].rank(pct=True) + d["chembl_max_phase"].rank(pct=True)
    top = d.nlargest(10, "_rank")
    ymax = 4.5
    label_ys = _push_apart(top["_y"].tolist(), min_gap=ymax * 0.06, lo=-0.3, hi=ymax)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        ax.annotate(row["gene_name"], xy=(row["ad_ot_score"], row["_y"]),
                    xytext=(row["ad_ot_score"] + 0.02, ly), fontsize=7.5, color=ANNOT_NAVY,
                    arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
                    ha="left", va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9))

    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_ylim(-0.4, ymax)
    ax.set_xlabel("AD OpenTargets association score", fontsize=9)
    ax.set_ylabel("ChEMBL max clinical phase", fontsize=9)
    ax.set_title(f"Priority view — master shortlist (n={len(d):,} genes)  |  "
                 "colour = protein change type  |  size = |Δ PSI|", fontsize=10, pad=8)
    _style(ax); _grid(ax)

    handles = [mpatches.Patch(color=CHANGE_COLORS[ct], label=ct) for ct in present]
    leg1 = ax.legend(handles=handles, title="protein_change_type",
                     bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7.5, title_fontsize=8)
    ax.add_artist(leg1)
    _size_legend(ax, d["delta_usage"].abs(), sizes, "|Δ PSI|", loc="lower left")

    save(fig, "07_priority_bubble.png", extra_artists=[leg1])


# ---------------------------------------------------------------------------
# Plot 08 — Genetic vs literature AD evidence
# ---------------------------------------------------------------------------

def plot_genetic_vs_lit(df: pd.DataFrame) -> None:
    """08 — ad_ot_genetic vs ad_ot_lit, colour = ad_ot_label, size = ad_specificity."""
    d = df.drop_duplicates("gene_name").copy()
    sizes = _size_scale(d["ad_specificity"])

    fig, ax = plt.subplots(figsize=(8, 7))
    lim = max(d["ad_ot_genetic"].max(), d["ad_ot_lit"].max()) * 1.08
    ax.plot([0, lim], [0, lim], color="#aaaaaa", lw=0.8, ls="--", zorder=1)

    for lbl in ["novel", "emerging", "supported"]:
        m = d["ad_ot_label"] == lbl
        if not m.any():
            continue
        ax.scatter(d.loc[m, "ad_ot_genetic"], d.loc[m, "ad_ot_lit"], s=sizes[m.values],
                   c=OT_LABEL_COLORS[lbl], alpha=0.8, linewidths=0.4,
                   edgecolors="white", zorder=2, rasterized=True, label=lbl)

    top = d.nlargest(8, "ad_specificity")
    label_ys = _push_apart(top["ad_ot_lit"].tolist(), min_gap=lim * 0.045, lo=0.0, hi=lim)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        ax.annotate(row["gene_name"], xy=(row["ad_ot_genetic"], row["ad_ot_lit"]),
                    xytext=(row["ad_ot_genetic"] + lim * 0.03, ly), fontsize=7.5, color=ANNOT_NAVY,
                    arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
                    ha="left", va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9))

    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("AD OT genetic-association evidence", fontsize=9)
    ax.set_ylabel("AD OT literature evidence", fontsize=9)
    ax.set_title(f"Genetic vs. literature AD support (n={len(d):,} genes)  |  "
                 "dashed line = equal evidence  |  size = AD specificity", fontsize=10, pad=8)
    _style(ax); _grid(ax)

    leg1 = ax.legend(title="ad_ot_label", bbox_to_anchor=(1.01, 1), loc="upper left",
                     fontsize=8, title_fontsize=8)
    ax.add_artist(leg1)
    _size_legend(ax, d["ad_specificity"], sizes, "AD specificity", loc="lower right")

    save(fig, "08_genetic_vs_lit.png", extra_artists=[leg1])


# ---------------------------------------------------------------------------
# Plot 09 — Volcano, upgraded: + domain count (colour) + edit size (bubble)
# ---------------------------------------------------------------------------

def plot_volcano_upgraded(df: pd.DataFrame) -> None:
    """09 — delta_usage vs -log10(padj), colour = pfam_n_domains, size = aa edit size."""
    d = df.copy()
    d["neg_log10_padj"] = -np.log10(d["chi_padj"].clip(lower=1e-300))
    d["_edit"] = d["mismatch_aa_count"] + d["indel_aa_count"]
    sizes = _size_scale(d["_edit"])
    ymax = d["neg_log10_padj"].max() * 1.10

    fig, ax = plt.subplots(figsize=(9.5, 7))
    sc = ax.scatter(d["delta_usage"], d["neg_log10_padj"], s=sizes,
                    c=d["pfam_n_domains"], cmap=NEON_CMAP, alpha=0.85,
                    linewidths=0.4, edgecolors="white", zorder=2, rasterized=True)

    ax.axhline(0, color="black", lw=0.6, ls="--", alpha=0.3)
    ax.axvline(0, color="black", lw=0.6, ls="--", alpha=0.3)
    ax.axhline(-np.log10(0.05), color="red", lw=0.8, ls="--", alpha=0.6)

    d["_rank"] = d["neg_log10_padj"].rank(pct=True) + d["_edit"].rank(pct=True)
    top = d.nlargest(10, "_rank")
    label_ys = _push_apart(top["neg_log10_padj"].tolist(), min_gap=ymax * 0.05, lo=0.0, hi=ymax)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        dx = 0.05 if row["delta_usage"] >= 0 else -0.05
        ax.annotate(row["gene_name"], xy=(row["delta_usage"], row["neg_log10_padj"]),
                    xytext=(row["delta_usage"] + dx, ly), fontsize=7.5, color=ANNOT_NAVY,
                    arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
                    ha="left" if dx > 0 else "right", va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9))

    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Δ usage (AD − Control)", fontsize=9)
    ax.set_ylabel("−log₁₀ (χ² adjusted p-value)", fontsize=9)
    ax.set_title(f"Volcano, upgraded (n={len(d):,})  |  colour = Pfam domain count  |  "
                 "size = aa mismatch + indel", fontsize=10, pad=8)
    _style(ax); _grid(ax)

    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("Pfam domains on alt isoform", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    # Anchor to the colorbar's own axes (not the figure) so the legend sits
    # directly under it regardless of how much canvas bbox_inches="tight" adds.
    size_leg = _size_legend(cb.ax, d["_edit"], sizes, "aa changed", loc="upper left",
                            fmt="{:.0f}", bbox_to_anchor=(0.0, -0.06))

    save(fig, "09_volcano_upgraded.png", extra_artists=[size_leg])


# ---------------------------------------------------------------------------
# Plot 10 — Drug tractability landscape
# ---------------------------------------------------------------------------

def plot_drug_landscape(df: pd.DataFrame) -> None:
    """10 — chembl_max_phase vs dgidb_interactions (log), colour = ad_ot_score, size = |delta_usage|."""
    d = _rep_per_gene(df)
    rng = np.random.default_rng(10)
    d = d.assign(_x=d["chembl_max_phase"].astype(float) + rng.uniform(-0.12, 0.12, len(d)))
    sizes = _size_scale(d["delta_usage"].abs())

    fig, ax = plt.subplots(figsize=(9, 7))
    sc = ax.scatter(d["_x"], d["dgidb_interactions"], s=sizes,
                    c=d["ad_ot_score"], cmap=NEON_CMAP, alpha=0.85,
                    linewidths=0.4, edgecolors="white", zorder=2, rasterized=True)
    ax.set_yscale("log")

    # Log-scale y: push_apart and label placement both need to work in
    # log10 space, or "equal gaps" on screen come out uneven / labels collide.
    ymin_data, ymax_data = 1.0, d["dgidb_interactions"].max() * 2.2
    top = d.nlargest(8, "ad_ot_score")
    log_ys = np.log10(top["dgidb_interactions"].clip(lower=1))
    log_ys = _push_apart(log_ys.tolist(), min_gap=(np.log10(ymax_data) - 0) * 0.06,
                         lo=0.0, hi=np.log10(ymax_data))
    for (_, row), ly in zip(top.iterrows(), log_ys):
        ax.annotate(
            row["gene_name"], xy=(row["_x"], row["dgidb_interactions"]),
            xytext=(row["_x"] + 0.3, 10 ** ly),
            fontsize=7.5, color=ANNOT_NAVY,
            arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
            ha="left", va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9),
        )

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlim(-0.5, 5.4)
    ax.set_ylim(ymin_data, ymax_data)
    ax.set_xlabel("ChEMBL max clinical phase  (0 = no compound on record)", fontsize=9)
    ax.set_ylabel("DGIdb interaction count (log scale)", fontsize=9)
    ax.set_title(f"Drug tractability landscape (n={len(d):,} genes)  |  "
                 "colour = AD OT score  |  size = |Δ PSI|", fontsize=10, pad=8)
    _style(ax)
    ax.yaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)

    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("AD OpenTargets score", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    _size_legend(ax, d["delta_usage"].abs(), sizes, "|Δ PSI|", loc="upper left")

    save(fig, "10_drug_landscape.png")


# ---------------------------------------------------------------------------
# Plot 11 — Confidence vs consequence
# ---------------------------------------------------------------------------

def plot_confidence_vs_consequence(df: pd.DataFrame) -> None:
    """11 — pct_identity vs -log10(padj), colour = change type, size = domains lost."""
    d = _rep_per_transcript(df)
    d["neg_log10_padj"] = -np.log10(d["chi_padj"].clip(lower=1e-300))
    d["n_lost"] = d["domains_lost"].apply(lambda s: len(_split_domains(s)))
    sizes = _size_scale(d["n_lost"])
    ymax = d["neg_log10_padj"].max() * 1.10

    fig, ax = plt.subplots(figsize=(9, 7))
    present = [ct for ct in _CHANGE_ORDER if ct in d["protein_change_type"].unique()]
    for ct in present:
        m = d["protein_change_type"] == ct
        ax.scatter(d.loc[m, "pct_identity"], d.loc[m, "neg_log10_padj"], s=sizes[m.values],
                   c=CHANGE_COLORS[ct], alpha=0.85, linewidths=0.4,
                   edgecolors="white", zorder=2, rasterized=True)

    d["_rank"] = (1 - d["pct_identity"]).rank(pct=True) + d["neg_log10_padj"].rank(pct=True)
    top = d.nlargest(10, "_rank")
    label_ys = _push_apart(top["neg_log10_padj"].tolist(), min_gap=ymax * 0.05, lo=0.0, hi=ymax)
    for (_, row), ly in zip(top.iterrows(), label_ys):
        ax.annotate(row["gene_name"], xy=(row["pct_identity"], row["neg_log10_padj"]),
                    xytext=(row["pct_identity"] - 0.015, ly), fontsize=7.5, color=ANNOT_NAVY,
                    arrowprops=dict(arrowstyle="-", color=ANNOT_NAVY, lw=0.6),
                    ha="right", va="center", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor=ANNOT_NAVY, linewidth=0.6, alpha=0.9))

    ax.set_xlim(1.02, d["pct_identity"].min() - 0.03)  # reversed: divergence increases rightward
    ax.set_ylim(0, ymax)
    ax.set_xlabel("← canonical/alt protein identity  (1.0 = identical)", fontsize=9)
    ax.set_ylabel("−log₁₀ (χ² adjusted p-value)", fontsize=9)
    ax.set_title(f"Confidence vs. consequence (n={len(d):,} transcripts)  |  "
                 "colour = protein change type  |  size = domains lost", fontsize=10, pad=8)
    _style(ax); _grid(ax)

    handles = [mpatches.Patch(color=CHANGE_COLORS[ct], label=ct) for ct in present]
    leg1 = ax.legend(handles=handles, title="protein_change_type",
                     bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7.5, title_fontsize=8)
    ax.add_artist(leg1)
    _size_legend(ax, d["n_lost"], sizes, "domains lost", loc="lower left", fmt="{:.0f}")

    save(fig, "11_confidence_vs_consequence.png", extra_artists=[leg1])


# ---------------------------------------------------------------------------
# Plot 12 — Trial-failure candidates: PSI change, Control -> AD
# ---------------------------------------------------------------------------

CELL_SHORT = {
    "Excitatory_neuron": "Exc. neuron",
    "Inhibitory_neuron": "Inh. neuron",
    "Astrocyte":         "Astrocyte",
    "Oligodendrocyte":   "Oligodendrocyte",
    "OPC":               "OPC",
    "Microglia":         "Microglia",
    "Vascular_cell":     "Vascular",
    "Lymphocyte":        "Lymphocyte",
}
TF_CONTROL_COLOR = "#2a78d6"   # blue — matches Control/down-in-AD convention used elsewhere
TF_AD_COLOR      = "#e34948"   # red  — matches AD/up-in-AD convention used elsewhere


def plot_tf_psi_dumbbell(df: pd.DataFrame) -> None:
    """12 — Trial-failure candidates: dominant/canonical isoform usage, Control -> AD.

    Dumbbell chart, one line per candidate. All 18 are CT_enriched by
    construction (see initial_filter.py), so every line points the same
    direction (down) -- the point is the magnitude, not the direction.
    """
    tf = df[df["candidate_group"] == "trial_failure_candidate"].copy()
    tf["cell_short"] = tf["cell_type"].map(lambda c: CELL_SHORT.get(c, c))
    tf = tf.reindex(tf["delta_usage"].abs().sort_values(ascending=False).index)  # biggest |Δ PSI| first
    y = np.arange(len(tf))[::-1]        # biggest decline drawn at the top

    fig, ax = plt.subplots(figsize=(7.5, 0.27 * len(tf) + 1.6))

    # Highlight rows with a large PSI swing or near-total loss in AD.
    highlight = (tf["delta_usage"].abs() > 0.50) | (tf["AD"] < 0.10)
    for yi, hl in zip(y, highlight):
        if hl:
            ax.axhspan(yi - 0.5, yi + 0.5, facecolor="yellow", edgecolor="none",
                       alpha=0.35, zorder=-1)

    for yi, ctrl, ad in zip(y, tf["Control"], tf["AD"]):
        ax.annotate(
            "", xy=(ad, yi), xytext=(ctrl, yi),
            arrowprops=dict(arrowstyle="->", color="#666666", lw=1.0,
                             mutation_scale=14, shrinkA=6, shrinkB=6),
            zorder=2,
        )
    ax.scatter(tf["Control"], y, s=55, color=TF_CONTROL_COLOR, zorder=3,
               marker="s", edgecolors="white", linewidths=0.6, label="Control")
    ax.scatter(tf["AD"], y, s=70, color=TF_AD_COLOR, zorder=3,
               edgecolors="white", linewidths=0.6, label="AD")

    MONO = "Liberation Mono"  # tabular-figure alignment for the numeric columns

    for yi, ctrl, ad in zip(y, tf["Control"], tf["AD"]):
        ax.text(ctrl + 0.02, yi, f"{ctrl * 100:.0f}%", ha="left", va="center",
                 fontsize=7.5, color=FG, fontfamily=MONO)
        ax.text(ad - 0.02, yi, f"{ad * 100:.0f}%", ha="right", va="center",
                 fontsize=7.5, color=FG, fontfamily=MONO)

    ax.set_yticks(y)
    ax.set_yticklabels([""] * len(tf))
    ax.set_xlim(-0.06, 1.12)
    ax.set_ylim(-0.7, len(tf) - 0.3)
    ax.xaxis.set_major_locator(plt.FixedLocator([0, 0.2, 0.4, 0.6, 0.8, 1.0]))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v * 100:.0f}%"))
    plt.setp(ax.get_xticklabels(), fontfamily=MONO)
    ax.set_xlabel("Dominant/canonical isoform usage (PSI)", fontsize=9, fontfamily="Liberation Sans")
    ax.set_title(
        f"Trial-failure candidates: isoform usage collapses in AD  (n={len(tf)})",
        fontsize=11, loc="left", pad=10, fontfamily="Liberation Sans",
    )

    _style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(False)

    handles = [
        mpatches.Patch(color=TF_CONTROL_COLOR, label="Control"),
        mpatches.Patch(color=TF_AD_COLOR, label="AD"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=8,
              frameon=True, framealpha=1, edgecolor=GRID)

    _draw_bold_transcript_italic_celltype_labels(fig, ax, y, tf["transcript_name"], tf["cell_short"])

    fig.tight_layout()
    save(fig, "12_tf_psi_dumbbell.png")


WRAP_TX_LEN = 15  # transcript names longer than this wrap to their own line


def _draw_bold_transcript_italic_celltype_labels(fig, ax, y, transcript_names, cell_shorts) -> None:
    """Two-weight row labels (bold transcript name + italic cell type).

    matplotlib tick labels can't mix styles in one string, so each part is
    drawn as its own Text object rather than relying on set_yticklabels.

    Short names (e.g. "GENE-201") stay on one line: the italic cell-type
    suffix is drawn first, measured, and the bold name placed immediately
    to its left (mirrors the group-header technique in
    plot_dominance_x_mane.py). Long names (novel long-read transcripts like
    "transcriptNNNNN.chrN.nic") would collide with the row's PSI value text
    if kept on one line, so they wrap: bold name on its own line, italic
    cell type stacked directly below it, both right-aligned at the same x.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    row_trans = blended_transform_factory(ax.transAxes, ax.transData)
    for yi, tx, cell in zip(y, transcript_names, cell_shorts):
        if len(tx) > WRAP_TX_LEN:
            ax.text(-0.015, yi + 0.24, tx, transform=row_trans, ha="right", va="center",
                    fontsize=8.5, fontweight="bold", fontfamily="Liberation Sans",
                    color=FG, clip_on=False)
            ax.text(-0.015, yi - 0.24, f"({cell})", transform=row_trans, ha="right", va="center",
                    fontsize=7.5, fontstyle="italic", fontfamily="Liberation Sans",
                    color=FG, clip_on=False)
            continue

        t_suffix = ax.text(-0.015, yi, f" ({cell})", transform=row_trans,
                            ha="right", va="center", fontsize=8.5, fontstyle="italic",
                            fontfamily="Liberation Sans", color=FG, clip_on=False)
        fig.canvas.draw()
        bbox = t_suffix.get_window_extent(renderer=renderer)
        x0_axes = ax.transAxes.inverted().transform((bbox.x0, 0))[0]
        ax.text(x0_axes, yi, tx, transform=row_trans, ha="right", va="center",
                fontsize=8.5, fontweight="bold", fontfamily="Liberation Sans",
                color=FG, clip_on=False)


# ---------------------------------------------------------------------------
# Plot 13 — New-target candidates: PSI change, Control -> AD
# ---------------------------------------------------------------------------

def plot_nt_psi_dumbbell(df: pd.DataFrame) -> None:
    """13 — New-target candidates: alternate isoform usage, Control -> AD.

    Dumbbell chart, one line per candidate (n=79, full list -- not top-N).
    All are AD_enriched by construction (see initial_filter.py), so every
    arrow points the same direction (up/right) -- opposite of the
    trial-failure plot, where the dominant/canonical isoform LOSES usage.
    Marker shapes are swapped vs. the trial-failure plot (AD=square,
    Control=circle) so the two plots are never visually confusable at a
    glance, even though they share the same colour convention.
    """
    nt = df[df["candidate_group"] == "new_target_candidate"].copy()
    nt["cell_short"] = nt["cell_type"].map(lambda c: CELL_SHORT.get(c, c))
    nt = nt.reindex(nt["delta_usage"].abs().sort_values(ascending=False).index)  # biggest |Δ PSI| first
    y = np.arange(len(nt))[::-1]        # biggest gain drawn at the top

    fig, ax = plt.subplots(figsize=(7.5, 0.27 * len(nt) + 1.6))

    # Highlight rows with a large PSI swing or near-total absence in Control
    # -- the mirror image of the trial-failure plot's AD<0.10 condition.
    highlight = (nt["delta_usage"].abs() > 0.50) | (nt["Control"] < 0.10)
    for yi, hl in zip(y, highlight):
        if hl:
            ax.axhspan(yi - 0.5, yi + 0.5, facecolor="yellow", edgecolor="none",
                       alpha=0.35, zorder=-1)

    for yi, ctrl, ad in zip(y, nt["Control"], nt["AD"]):
        ax.annotate(
            "", xy=(ad, yi), xytext=(ctrl, yi),
            arrowprops=dict(arrowstyle="->", color="#666666", lw=1.0,
                             mutation_scale=14, shrinkA=6, shrinkB=6),
            zorder=2,
        )
    ax.scatter(nt["Control"], y, s=70, color=TF_CONTROL_COLOR, zorder=3,
               edgecolors="white", linewidths=0.6, label="Control")
    ax.scatter(nt["AD"], y, s=55, color=TF_AD_COLOR, zorder=3,
               marker="s", edgecolors="white", linewidths=0.6, label="AD")

    MONO = "Liberation Mono"  # tabular-figure alignment for the numeric columns

    for yi, ctrl, ad in zip(y, nt["Control"], nt["AD"]):
        ax.text(ctrl - 0.02, yi, f"{ctrl * 100:.0f}%", ha="right", va="center",
                 fontsize=7.5, color=FG, fontfamily=MONO)
        ax.text(ad + 0.02, yi, f"{ad * 100:.0f}%", ha="left", va="center",
                 fontsize=7.5, color=FG, fontfamily=MONO)

    ax.set_yticks(y)
    ax.set_yticklabels([""] * len(nt))
    ax.set_xlim(-0.06, 1.12)
    ax.set_ylim(-0.7, len(nt) - 0.3)
    ax.xaxis.set_major_locator(plt.FixedLocator([0, 0.2, 0.4, 0.6, 0.8, 1.0]))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v * 100:.0f}%"))
    plt.setp(ax.get_xticklabels(), fontfamily=MONO)
    ax.set_xlabel("Alternate/minor isoform usage (PSI)", fontsize=9, fontfamily="Liberation Sans")
    ax.set_title(
        f"New-target candidates: alternate isoform usage rises in AD  (n={len(nt)})",
        fontsize=11, loc="left", pad=10, fontfamily="Liberation Sans",
    )

    _style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.yaxis.grid(False)

    handles = [
        mpatches.Patch(color=TF_CONTROL_COLOR, label="Control"),
        mpatches.Patch(color=TF_AD_COLOR, label="AD"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8,
              frameon=True, framealpha=1, edgecolor=GRID)

    _draw_bold_transcript_italic_celltype_labels(fig, ax, y, nt["transcript_name"], nt["cell_short"])

    fig.tight_layout()
    save(fig, "13_nt_psi_dumbbell.png")


# ---------------------------------------------------------------------------
# Plot 14 — Both candidate groups: Control PSI vs AD PSI scatter
# ---------------------------------------------------------------------------

GROUP_COLORS = {
    "trial_failure_candidate": TF_CONTROL_COLOR,  # blue
    "new_target_candidate":    TF_AD_COLOR,       # red
}
GROUP_LABELS = {
    "trial_failure_candidate": "Trial-failure candidate",
    "new_target_candidate":    "New-target candidate",
}


def plot_candidate_scatter(df: pd.DataFrame) -> None:
    """14 — Control PSI vs AD PSI for both candidate groups, one scatter.

    Trial-failure candidates fall below the y=x line (AD < Control);
    new-target candidates fall above it (AD > Control) -- that's guaranteed
    by construction (initial_filter.py), not a finding. The dotted
    y=x+-50 lines mark the |delta PSI| = 50 boundary used to highlight rows
    in the two dumbbell plots (12/13).
    """
    d = df[df["candidate_group"].isin(GROUP_LABELS)].copy()
    d["ct_pct"] = d["Control"] * 100
    d["ad_pct"] = d["AD"] * 100
    d["abs_aa_change"] = d["mismatch_aa_count"] + d["indel_aa_count"]
    sizes = _size_scale(d["abs_aa_change"], lo=25, hi=350)

    fig, ax = plt.subplots(figsize=(7, 7))
    # Fix the coordinate system (limits + aspect) before any pixel-space text
    # measurement/placement below -- transData must be stable at draw time,
    # or positions computed against it get invalidated by a later aspect change.
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("equal")

    xs = np.array([0, 100])
    ax.plot(xs, xs, color="#666666", lw=1.2, zorder=2, label="_nolegend_")
    ax.plot(xs, xs - 50, color="#d03b3b", lw=1.2, linestyle=":", zorder=2)
    ax.plot(xs, xs + 50, color="#d03b3b", lw=1.2, linestyle=":", zorder=2)

    # Perpendicular-to-diagonal indicator: |delta PSI| grows moving away from
    # y=x in either direction. Anchored at the diagonal's midpoint, which is
    # empty of points by construction (every candidate sits well to one side).
    center = 50.0
    arrow_len = 9.0
    step = arrow_len / np.sqrt(2)
    REGION_GREY = "#999999"
    region_text = {1: "Usage increased in AD", -1: "Usage decreased in AD"}
    for sign in (1, -1):
        tip = (center - sign * step, center + sign * step)
        ax.annotate("", xy=tip, xytext=(center, center),
                    arrowprops=dict(arrowstyle="->", color="#555555", lw=1.3, mutation_scale=12),
                    zorder=4)
        label_xy = (center - sign * step * 1.15, center + sign * step * 1.15)
        ax.text(*label_xy, region_text[sign], ha="center", va="center", rotation=45,
                 fontsize=11, color=REGION_GREY, fontfamily="Liberation Sans", zorder=4)

    for grp, label in GROUP_LABELS.items():
        mask = d["candidate_group"] == grp
        ax.scatter(d.loc[mask, "ct_pct"], d.loc[mask, "ad_pct"], s=sizes[mask.values],
                   color=GROUP_COLORS[grp], alpha=0.85, edgecolors="white",
                   linewidths=0.6, zorder=3, label=label)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_xlabel("Control PSI", fontsize=9, fontfamily="Liberation Sans")
    ax.set_ylabel("AD PSI", fontsize=9, fontfamily="Liberation Sans")
    ax.set_title(
        f"Trial-failure vs. new-target candidates: Control vs. AD isoform usage  (n={len(d)})",
        fontsize=10.5, loc="left", pad=10, fontfamily="Liberation Sans",
    )

    _style(ax)
    _grid(ax)
    plt.setp(ax.get_xticklabels(), fontfamily="Liberation Sans")
    plt.setp(ax.get_yticklabels(), fontfamily="Liberation Sans")
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8,
                     frameon=True, framealpha=1, edgecolor=GRID)
    for text in leg.get_texts():
        text.set_fontfamily("Liberation Sans")
    ax.add_artist(leg)
    size_leg = _size_legend(ax, d["abs_aa_change"], sizes, "abs. AA changed",
                             loc="upper left", fmt="{:.0f}", fig=fig,
                             bbox_to_anchor=(1.02, 0.75))

    fig.tight_layout()
    save(fig, "14_candidate_scatter_ct_ad.png", extra_artists=[leg, size_leg])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading hits_deep.csv (selected_for_next_stage only) …")
    df = load()
    print(f"  {len(df):,} rows  |  {df['gene_name'].nunique():,} genes  |  "
          f"{df['transcript_name'].nunique():,} transcripts\n")

    print("Rendering plots …")
    plot_volcano(df)
    plot_change_type_bar(df)
    plot_domain_gain_loss(df)
    plot_drug_evidence(df)
    plot_ad_specificity_vs_score(df)
    plot_scorecard_heatmap(df)
    plot_priority_bubble(df)
    plot_genetic_vs_lit(df)
    plot_volcano_upgraded(df)
    plot_drug_landscape(df)
    plot_confidence_vs_consequence(df)
    plot_tf_psi_dumbbell(df)
    plot_nt_psi_dumbbell(df)
    plot_candidate_scatter(df)

    n = len(list(OUT_DIR.glob("*.png")))
    print(f"\nDone — {n} plot(s) in {OUT_DIR}")


if __name__ == "__main__":
    main()

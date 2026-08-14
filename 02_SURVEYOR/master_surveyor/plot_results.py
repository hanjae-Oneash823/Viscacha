"""MASTER_SURVEYOR — preliminary visualization suite for the 228-hit J4 shortlist."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

# Ensure 02_SURVEYOR (this package's parent) is on path regardless of
# invocation cwd.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor.config import HITS_CSV, OUT_DIR
from master_surveyor.m0_select import representative_row as _tf_representative_row

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

BG         = "#ffffff"
FG         = "#1a1a1a"
GRID       = "#e0e0e0"
ANNOT_NAVY = "#003399"

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
    """Load junior_surveyor's hits_deep.csv, restricted to selected_for_next_stage.

    trial_failure_candidate hits are long-format (one row per ranked
    alternate isoform, all sharing the same selected_for_next_stage value) --
    collapsed here to one representative row per hit, the gate-driving
    alternate with the largest AD usage gain, so downstream plots (which
    assume one row per hit/transcript) don't show a hit duplicated across
    its N ranked alternates. new_target_candidate hits are already one row
    per hit (alt_rank=0) and pass through unchanged.
    """
    df = pd.read_csv(HITS_CSV)
    df = df[df["selected_for_next_stage"]].copy()

    nt = df[df["candidate_group"] == "new_target_candidate"]
    tf = df[df["candidate_group"] == "trial_failure_candidate"]
    if not tf.empty:
        tf = (tf.groupby(["gene_name", "cell_type"], group_keys=False)
                .apply(_tf_representative_row))
    df = pd.concat([nt, tf], ignore_index=True, sort=False)

    df["protein_change_type"] = df["protein_change_type"].fillna("no_sequence")
    df["affected_domain"]     = df["affected_domain"].fillna("none")
    df["chembl_max_phase"]    = df["chembl_max_phase"].fillna(0).astype(int)
    df["dgidb_interactions"]  = df["dgidb_interactions"].fillna(0).astype(int)
    return df


# ---------------------------------------------------------------------------
# Plot 01 — Protein change type distribution
# ---------------------------------------------------------------------------

def plot_change_type_bar(df: pd.DataFrame) -> None:
    """01 — Distribution of protein_change_type (deduped per transcript)."""
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

    save(fig, "01_change_type_bar.png")


# ---------------------------------------------------------------------------
# Plot 02 — Drug / target evidence
# ---------------------------------------------------------------------------

def plot_drug_evidence(df: pd.DataFrame) -> None:
    """02 — ChEMBL max phase vs DGIdb interactions, one point per gene."""
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

    save(fig, "02_drug_evidence_scatter.png")


# ---------------------------------------------------------------------------
# Plot 03 — Volcano, upgraded: + domain count (colour) + edit size (bubble)
# ---------------------------------------------------------------------------

def plot_volcano_upgraded(df: pd.DataFrame) -> None:
    """03 — delta_usage vs -log10(padj), colour = pfam_n_domains, size = aa edit size."""
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

    save(fig, "03_volcano_upgraded.png", extra_artists=[size_leg])


# ---------------------------------------------------------------------------
# Plot 04 — Trial-failure candidates: PSI change, Control -> AD
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
    """04 — Trial-failure candidates: dominant/canonical isoform usage, Control -> AD.

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
    save(fig, "04_tf_psi_dumbbell.png")


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
# Plot 05 — New-target candidates: PSI change, Control -> AD
# ---------------------------------------------------------------------------

def plot_nt_psi_dumbbell(df: pd.DataFrame) -> None:
    """05 — New-target candidates: alternate isoform usage, Control -> AD.

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
    save(fig, "05_nt_psi_dumbbell.png")


# ---------------------------------------------------------------------------
# Plot 06 — Both candidate groups: Control PSI vs AD PSI scatter
# ---------------------------------------------------------------------------

NOVEL_TARGET_COLOR = "#e8a628"   # amber — distinct from the TF/repurposing red-blue pair

GROUP_COLORS = {
    "trial_failure_candidate":    TF_CONTROL_COLOR,     # blue
    "drug_repurposing_candidate": TF_AD_COLOR,          # red
    "novel_target_candidate":     NOVEL_TARGET_COLOR,   # amber
}
GROUP_LABELS = {
    "trial_failure_candidate":    "Trial-failure candidate",
    "drug_repurposing_candidate": "Drug-repurposing candidate",
    "novel_target_candidate":     "Novel-target candidate",
}


def plot_candidate_scatter(df: pd.DataFrame) -> None:
    """06 — Control PSI vs AD PSI for the three master_group candidates, one scatter.

    Trial-failure candidates fall below the y=x line (AD < Control);
    the two new_target-derived groups (drug_repurposing / novel_target) fall
    above it (AD > Control) -- that's guaranteed by construction
    (initial_filter.py), not a finding. The dotted y=x+-50 lines mark the
    |delta PSI| = 50 boundary used to highlight rows in the two dumbbell
    plots (04/05). drug_repurposing vs novel_target is master_group's
    post-hoc split of new_target_candidate by chemical-tractability evidence
    (see j4_gate.py) -- it doesn't affect PSI position, only colour.
    """
    d = df[df["master_group"].isin(GROUP_LABELS)].copy()
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
        mask = d["master_group"] == grp
        ax.scatter(d.loc[mask, "ct_pct"], d.loc[mask, "ad_pct"], s=sizes[mask.values],
                   color=GROUP_COLORS[grp], alpha=0.85, edgecolors="white",
                   linewidths=0.6, zorder=3, label=label)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_xlabel("Control PSI", fontsize=9, fontfamily="Liberation Sans")
    ax.set_ylabel("AD PSI", fontsize=9, fontfamily="Liberation Sans")
    ax.set_title(
        f"Trial-failure vs. drug-repurposing vs. novel-target candidates: "
        f"Control vs. AD isoform usage  (n={len(d)})",
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
    save(fig, "06_candidate_scatter_ct_ad.png", extra_artists=[leg, size_leg])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading hits_deep.csv (selected_for_next_stage only) …")
    df = load()
    print(f"  {len(df):,} rows  |  {df['gene_name'].nunique():,} genes  |  "
          f"{df['transcript_name'].nunique():,} transcripts")
    for group in ["trial_failure_candidate", "drug_repurposing_candidate", "novel_target_candidate"]:
        n = (df["master_group"] == group).sum()
        print(f"    {group}: {n:,} hits")
    print()

    print("Rendering plots …")
    plot_change_type_bar(df)
    plot_drug_evidence(df)
    plot_volcano_upgraded(df)
    plot_tf_psi_dumbbell(df)
    plot_nt_psi_dumbbell(df)
    plot_candidate_scatter(df)

    n = len(list(OUT_DIR.glob("*.png")))
    print(f"\nDone — {n} plot(s) in {OUT_DIR}")


if __name__ == "__main__":
    main()

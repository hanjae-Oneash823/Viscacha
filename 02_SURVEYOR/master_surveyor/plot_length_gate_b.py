"""10b -- canonical vs alt protein length, standalone Gate B view.

Same data and geometry as the right panel of plot_protein_change.py
(09_protein_change.png), pulled out on its own so it reads independently of
the changed-residue-count panel. Every point is one of the 94 distinct
canonical/alt protein pairs from gate_matrix.build_gate_matrix(); the
diagonal is "no length change", the red dashed line is the Gate B floor
(MIN_KEPT_FRAC, alt >= 50% of canonical length) -- points below it fail
Gate B on length alone, independent of the (now-dropped) CDS-annotation
clause.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from master_surveyor.gate_matrix import MIN_KEPT_FRAC, build_gate_matrix
from master_surveyor.plot_results import BG, CHANGE_COLORS, FG, GRID, _style, save

HILITE = "#C41B3A"
MUTED = "#6b6b6b"
GENE_LABEL_BOX = dict(boxstyle="round,pad=0.16", facecolor=BG,
                      edgecolor="#d7dbe0", linewidth=0.65)

# Extremes worth naming so the axis range is self-explanatory.
LENGTH_CALLOUTS = {
    "SGMS1-206": (1.55, 1.05),
    "CLK4-208": (0.30, 0.55),
    "DCLK1-205": (1.55, 0.62),
    "GABRA2-206": (0.26, 1.55),
    "BIRC6-215": (0.42, 1.30),
}

# Fixed data-space label slots for the continuous linear plot.  The four
# short-protein callouts occupy the same crowded lower-left region, so simple
# multiplicative offsets make their boxes collide.
LINEAR_LABEL_POSITIONS = {
    "GABRA2-206": (100, 220, "right"),
    "CLK4-208": (100, 70, "right"),
    "SGMS1-206": (720, 175, "left"),
    "DCLK1-205": (1180, 85, "left"),
    "BIRC6-215": (2050, 4950, "left"),
}


def plot(pairs: pd.DataFrame, log_scale: bool = True,
         max_length: float | None = None) -> None:
    d = pairs.copy()
    if max_length is not None:
        d = d[(d["can_aa_len"] <= max_length) & (d["alt_aa_len"] <= max_length)].copy()
    order = (d.groupby("protein_change_type")["n_changed"]
             .median().sort_values(ascending=False).index.tolist())

    fig, ax = plt.subplots(figsize=(8.6, 8.0), facecolor=BG)
    _style(ax)

    if log_scale:
        lo = min(d["alt_aa_len"].min(), d["can_aa_len"].min()) * 0.7
        hi = max(d["alt_aa_len"].max(), d["can_aa_len"].max()) * 1.4
    else:
        lo = 0.0
        hi = max_length if max_length is not None else (
            max(d["alt_aa_len"].max(), d["can_aa_len"].max()) * 1.05)
    line = np.array([lo, hi])

    ax.fill_between(line, lo, line * MIN_KEPT_FRAC, color="#f5f6f8", zorder=0,
                    linewidth=0)
    ax.plot(line, line, color=MUTED, linestyle=(0, (5, 4)), linewidth=1.1,
            zorder=1)
    ax.plot(line, line * MIN_KEPT_FRAC, color=HILITE, linestyle=(0, (6, 3)),
            linewidth=1.4, zorder=2)
    ax.text(hi * 0.92, hi * 0.92, "same length  ", fontsize=8.5, color=MUTED,
            ha="right", va="bottom")
    n_below = int((d["alt_aa_len"] / d["can_aa_len"] < MIN_KEPT_FRAC).sum())
    ax.text(hi * 0.92, hi * 0.92 * MIN_KEPT_FRAC,
            f"gate B: {MIN_KEPT_FRAC:.0%} kept  ", fontsize=9, color=HILITE,
            ha="right", va="bottom", fontweight="bold")

    for kind in order:
        sub = d[d["protein_change_type"] == kind]
        ax.scatter(sub["can_aa_len"], sub["alt_aa_len"], s=58,
                   color=CHANGE_COLORS.get(kind, MUTED), alpha=0.85,
                   linewidth=0.6, edgecolor=BG, zorder=3)

    for name, (fx, fy) in LENGTH_CALLOUTS.items():
        hit = d[d["alt_transcript_name"] == name]
        if hit.empty:
            continue
        r = hit.iloc[0]
        # Multiplicative offsets are tuned for the log-log layout; on linear
        # axes the same factors can push labels off-canvas (e.g. BIRC6-215
        # near the top-right corner), so clamp into the visible range.
        if log_scale:
            tx = min(r["can_aa_len"] * fx, hi * 0.97)
            ty = min(r["alt_aa_len"] * fy, hi * 0.97)
            ha = "left" if fx >= 1 else "right"
        else:
            tx, ty, ha = LINEAR_LABEL_POSITIONS.get(
                name, (r["can_aa_len"] * fx, r["alt_aa_len"] * fy,
                       "left" if fx >= 1 else "right"))
        ax.annotate(name, xy=(r["can_aa_len"], r["alt_aa_len"]),
                    xytext=(tx, ty),
                    fontsize=8.5, color=FG, va="center",
                    ha=ha, zorder=5, bbox=GENE_LABEL_BOX,
                    arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.7,
                                    shrinkA=0, shrinkB=3))

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("canonical protein length  (aa)", fontsize=10.5, color=FG)
    ax.set_ylabel("alt protein length  (aa)", fontsize=10.5, color=FG)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    fig.suptitle(
        f"Gate B — canonical vs alt protein length "
        f"({n_below} of {len(d)} pairs fall below the {MIN_KEPT_FRAC:.0%} floor)",
        fontsize=14.5, color=FG, fontweight="bold", x=0.012, ha="left", y=0.995,
    )
    fig.text(0.012, 0.955,
             f"{len(d)} trial-failure protein pairs ({d['gene_name'].nunique()} genes) "
             f"— colored by protein_change_type"
             + (f" — ≤ {max_length:,.0f} aa on both axes" if max_length else ""),
             fontsize=9.5, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=9,
               color=CHANGE_COLORS.get(k, MUTED), markeredgecolor=BG,
               label=k.replace("_", " "))
        for k in order
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=min(len(order), 5), frameon=False, fontsize=9,
               handletextpad=0.4, columnspacing=1.4)

    fig.subplots_adjust(left=0.11, right=0.97, top=0.90, bottom=0.16)
    name = (
        "10b_length_gate_b.png" if log_scale else
        "10b_length_gate_b_linear_0_2000aa.png" if max_length == 2000 else
        "10b_length_gate_b_linear.png"
    )
    save(fig, name)


def plot_linear_broken(pairs: pd.DataFrame,
                        break_lo: float = 2450,
                        x_break_hi: float = 4700,
                        y_break_hi: float = 3700,
                        visual_gap: float = 150) -> None:
    """Linear canonical-vs-alt length on one axes with scale breaks.

    The empty x and y ranges are compressed to a narrow visual gap.  This
    preserves a conventional single-frame scatter plot while marking that
    distances across each gap are not on the same linear scale.
    """
    d = pairs.copy()
    order = (d.groupby("protein_change_type")["n_changed"]
             .median().sort_values(ascending=False).index.tolist())
    lo = 0.0
    hi_max = max(d["alt_aa_len"].max(), d["can_aa_len"].max()) * 1.05
    def compress(values: pd.Series | np.ndarray | float, gap_hi: float) -> np.ndarray:
        values = np.asarray(values, dtype=float)
        return np.where(values <= break_lo, values,
                        break_lo + visual_gap + values - gap_hi)

    fig, ax = plt.subplots(figsize=(8.8, 8.2), facecolor=BG)
    _style(ax)
    x_hi = break_lo + visual_gap + hi_max - x_break_hi
    y_hi = break_lo + visual_gap + hi_max - y_break_hi

    # Draw lines only on their valid pieces.  Connecting them across a break
    # would falsely imply that the elided range is linearly represented.
    low_line = np.array([lo, break_lo])
    high_diagonal = np.array([max(x_break_hi, y_break_hi), hi_max])
    ax.fill_between(low_line, lo, low_line * MIN_KEPT_FRAC,
                    color="#f5f6f8", zorder=0, linewidth=0)
    ax.plot(low_line, low_line, color=MUTED, linestyle=(0, (5, 4)),
            linewidth=1.1, zorder=1)
    ax.plot(compress(high_diagonal, x_break_hi),
            compress(high_diagonal, y_break_hi), color=MUTED,
            linestyle=(0, (5, 4)), linewidth=1.1, zorder=1)
    ax.plot(low_line, low_line * MIN_KEPT_FRAC, color=HILITE,
            linestyle=(0, (6, 3)), linewidth=1.4, zorder=2)
    # This portion of the Gate-B floor remains visible after the x break,
    # until it reaches the omitted y range.
    gate_hi_x = np.array([x_break_hi, break_lo / MIN_KEPT_FRAC])
    gate_hi_y = gate_hi_x * MIN_KEPT_FRAC
    ax.fill_between(compress(gate_hi_x, x_break_hi), lo,
                    compress(gate_hi_y, y_break_hi), color="#f5f6f8",
                    zorder=0, linewidth=0)
    ax.plot(compress(gate_hi_x, x_break_hi),
            compress(gate_hi_y, y_break_hi), color=HILITE,
            linestyle=(0, (6, 3)), linewidth=1.4, zorder=2)

    for kind in order:
        sub = d[d["protein_change_type"] == kind]
        ax.scatter(compress(sub["can_aa_len"], x_break_hi),
                   compress(sub["alt_aa_len"], y_break_hi), s=58,
                   color=CHANGE_COLORS.get(kind, MUTED), alpha=0.85,
                   linewidth=0.6, edgecolor=BG, zorder=3)

    # Pale bands plus paired slashes make the compressed intervals explicit
    # without splitting the chart into diagonal panels.
    gap_color = "#f7f8fa"
    ax.axvspan(break_lo, break_lo + visual_gap, color=gap_color, zorder=4)
    ax.axhspan(break_lo, break_lo + visual_gap, color=gap_color, zorder=4)
    slash = visual_gap * 0.14
    edge_offset = 0.014
    x_mid = break_lo + visual_gap / 2
    y_mid = break_lo + visual_gap / 2
    for off in (-slash, slash):
        ax.plot([x_mid - slash + off, x_mid + slash + off],
                [-edge_offset, edge_offset],
                transform=ax.get_xaxis_transform(), color=MUTED,
                linewidth=1.3, clip_on=False, zorder=5)
        ax.plot([-edge_offset, edge_offset], [y_mid - slash + off, y_mid + slash + off],
                transform=ax.get_yaxis_transform(), color=MUTED,
                linewidth=1.3, clip_on=False, zorder=5)

    x_low_ticks = np.arange(0, break_lo, 500)
    x_high_ticks = np.arange(np.ceil(x_break_hi / 200) * 200, hi_max + 1, 200)
    y_low_ticks = np.arange(0, break_lo, 500)
    y_high_ticks = np.arange(np.ceil(y_break_hi / 200) * 200, hi_max + 1, 200)
    ax.set_xticks(np.r_[x_low_ticks, compress(x_high_ticks, x_break_hi)],
                  [f"{x:.0f}" for x in np.r_[x_low_ticks, x_high_ticks]])
    ax.set_yticks(np.r_[y_low_ticks, compress(y_high_ticks, y_break_hi)],
                  [f"{y:.0f}" for y in np.r_[y_low_ticks, y_high_ticks]])
    ax.set_xlim(lo, x_hi)
    ax.set_ylim(lo, y_hi)
    ax.set_xlabel("canonical protein length  (aa)", fontsize=10.5, color=FG)
    ax.set_ylabel("alt protein length  (aa)", fontsize=10.5, color=FG)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    label_pos = hi_max * 0.97
    ax.text(compress(label_pos, x_break_hi), compress(label_pos, y_break_hi),
            "same length  ", fontsize=8.5, color=MUTED, ha="right", va="bottom")
    n_below = int((d["alt_aa_len"] / d["can_aa_len"] < MIN_KEPT_FRAC).sum())
    ax.text(break_lo * 0.92, break_lo * 0.92 * MIN_KEPT_FRAC,
            f"gate B: {MIN_KEPT_FRAC:.0%} kept  ", fontsize=9, color=HILITE,
            ha="right", va="bottom", fontweight="bold")

    for name, (fx, fy) in LENGTH_CALLOUTS.items():
        hit = d[d["alt_transcript_name"] == name]
        if hit.empty:
            continue
        r = hit.iloc[0]
        tx = min(r["can_aa_len"] * fx, hi_max * 0.97)
        ty = min(r["alt_aa_len"] * fy, hi_max * 0.97)
        ax.annotate(name,
                    xy=(compress(r["can_aa_len"], x_break_hi),
                        compress(r["alt_aa_len"], y_break_hi)),
                    xytext=(compress(tx, x_break_hi), compress(ty, y_break_hi)),
                    fontsize=8.5, color=FG, va="center",
                    ha="left" if fx >= 1 else "right", zorder=5, bbox=GENE_LABEL_BOX,
                    arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.7,
                                    shrinkA=0, shrinkB=3))

    fig.suptitle(
        f"Gate B — canonical vs alt protein length "
        f"({n_below} of {len(d)} pairs fall below the {MIN_KEPT_FRAC:.0%} floor)",
        fontsize=14.5, color=FG, fontweight="bold", x=0.012, ha="left", y=0.995,
    )
    fig.text(0.012, 0.955,
             f"{len(d)} trial-failure protein pairs ({d['gene_name'].nunique()} genes) "
             f"— colored by protein_change_type — scale breaks "
             f"x: {break_lo:.0f}-{x_break_hi:.0f}; "
             f"y: {break_lo:.0f}-{y_break_hi:.0f} aa, no pairs fall in those bands",
             fontsize=9.5, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=9,
               color=CHANGE_COLORS.get(k, MUTED), markeredgecolor=BG,
               label=k.replace("_", " "))
        for k in order
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=min(len(order), 5), frameon=False, fontsize=9,
               handletextpad=0.4, columnspacing=1.4)

    save(fig, "10b_length_gate_b_linear.png")


def main() -> None:
    pairs = build_gate_matrix()
    plot(pairs, log_scale=True)
    plot(pairs, log_scale=False)
    plot(pairs, log_scale=False, max_length=2000)

    kept = pairs["alt_aa_len"] / pairs["can_aa_len"]
    print(f"\npairs={len(pairs)} genes={pairs.gene_name.nunique()}")
    print(f"below gate B length floor ({MIN_KEPT_FRAC:.0%} kept): "
          f"{int((kept < MIN_KEPT_FRAC).sum())}")
    print(f"pass gate B (length + no premature stop): {int(pairs['B'].sum())}")


if __name__ == "__main__":
    main()

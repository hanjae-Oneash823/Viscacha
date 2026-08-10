"""09 -- how much of the protein actually changes in the alt isoform.

One point per docking pair from gate_matrix.build_gate_matrix() (94 distinct
canonical/alt protein pairs over 55 trial_failure genes; the canonical's own
row is excluded, so every pair here is a real comparison).

Left panel  -- fraction of canonical residues altered, by change class. This
is `true_changed_frac` (n_changed / canonical length) from the edlib CIGAR
walk, NOT the changed_aa_start..end envelope.

The classes separate almost perfectly, and that is the finding: the 55
truncation-like pairs change a median 44% of the sequence, while all 39
others (indels, insertions, extensions) top out at 15.8% and sit at a median
of 1.3%. So `protein_change_type` already tells you the magnitude -- reading
the percentage adds little once you know the class. Truncations themselves
span the full 1.6-97% range, so the class is the informative variable and the
percentage is not a usable one-dimensional ranking.

Right panel -- the same thing in absolute terms, canonical vs alt length. The
0.5x line is gate B: below it more than half the protein is gone and the
result is a fragment, not a reshaped target.

The two panels disagree usefully. A short protein losing 40% and a 4857-aa
protein losing 40% look identical on the left and nothing alike on the right.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from master_surveyor.gate_matrix import (
    MECHANISM_GATES, MIN_KEPT_FRAC, build_gate_matrix,
)
from master_surveyor.plot_results import (
    BG, CHANGE_COLORS, FG, GRID, _style, save,
)

HILITE = "#C41B3A"
MUTED = "#6b6b6b"
BAND = "#f5f6f8"

TRUNCATION_LIKE = ("N_truncation", "C_truncation", "frameshift_stop")

# Pairs worth naming on the length panel: the extremes that explain the axis.
# Offsets are hand-placed multipliers on (x, y) -- the callouts cluster in the
# same corner, so an automatic declutter just stacks them on each other.
LENGTH_CALLOUTS = {
    "SGMS1-206": (1.55, 1.05),
    "CLK4-208": (0.30, 0.55),
    "DCLK1-205": (1.55, 0.62),
    "GABRA2-206": (0.26, 1.55),
    "BIRC6-215": (0.42, 1.30),
}


def _strip_panel(ax, d: pd.DataFrame) -> list[str]:
    """Jittered per-pair strip of changed fraction, one row per change class."""
    order = (d.groupby("protein_change_type")["true_changed_frac"]
             .median().sort_values(ascending=False).index.tolist())
    rng = np.random.default_rng(0)

    for i, kind in enumerate(order):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=BAND, zorder=0, linewidth=0)
        sub = d[d["protein_change_type"] == kind]
        vals = sub["true_changed_frac"].to_numpy() * 100
        jitter = rng.uniform(-0.26, 0.26, len(vals))
        ax.scatter(vals, i + jitter, s=52, color=CHANGE_COLORS.get(kind, MUTED),
                   alpha=0.85, linewidth=0.6, edgecolor=BG, zorder=3)
        med = float(np.median(vals))
        ax.plot([med, med], [i - 0.36, i + 0.36], color=FG, linewidth=2.0,
                zorder=4, solid_capstyle="butt")
        ax.text(101.5, i, f"n={len(sub)}", fontsize=8.5, color=MUTED,
                va="center", ha="left")

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([k.replace("_", " ") for k in order], fontsize=10,
                       color=FG)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlim(-2, 100)
    ax.set_xlabel("canonical residues altered in the alt  (%)", fontsize=10.5,
                  color=FG)
    ax.set_title("A   how much of the sequence changes", fontsize=11.5,
                 color=FG, fontweight="bold", loc="left", pad=12)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=1)
    ax.yaxis.grid(False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)

    # Where the classes stop overlapping. Data-derived, not a chosen cutoff.
    split = d.loc[~d["protein_change_type"].isin(TRUNCATION_LIKE),
                  "true_changed_frac"].max() * 100
    ax.axvline(split, color=MUTED, linestyle=(0, (4, 4)), linewidth=1.1, zorder=2)
    # Sits in the block that is empty by construction: no non-truncation class
    # reaches this far right, which is exactly what the caption says.
    ax.text((split + 100) / 2, (len(order) - 1) * 0.72,
            f"no indel, insertion or extension pair changes\n"
            f"more than {split:.0f}% — class alone predicts the magnitude",
            fontsize=8.5, color=MUTED, ha="center", va="center", linespacing=1.6,
            style="italic", zorder=5)
    return order


def _length_panel(ax, d: pd.DataFrame, order: list[str]) -> None:
    lo = min(d["alt_aa_len"].min(), d["can_aa_len"].min()) * 0.7
    hi = max(d["alt_aa_len"].max(), d["can_aa_len"].max()) * 1.4
    line = np.array([lo, hi])

    ax.fill_between(line, lo, line * MIN_KEPT_FRAC, color=BAND, zorder=0,
                    linewidth=0)
    ax.plot(line, line, color=MUTED, linestyle=(0, (5, 4)), linewidth=1.1,
            zorder=1)
    ax.plot(line, line * MIN_KEPT_FRAC, color=HILITE, linestyle=(0, (6, 3)),
            linewidth=1.4, zorder=2)
    ax.text(hi * 0.92, hi * 0.92, "same length  ", fontsize=8.5, color=MUTED,
            ha="right", va="bottom", rotation=45, rotation_mode="anchor")
    n_below = int((d["alt_aa_len"] / d["can_aa_len"] < MIN_KEPT_FRAC).sum())
    ax.text(hi * 0.92, hi * 0.92 * MIN_KEPT_FRAC,
            f"gate B: {MIN_KEPT_FRAC:.0%} kept  ", fontsize=9, color=HILITE,
            ha="right", va="bottom", rotation=45, rotation_mode="anchor",
            fontweight="bold")

    for kind in order:
        sub = d[d["protein_change_type"] == kind]
        ax.scatter(sub["can_aa_len"], sub["alt_aa_len"], s=52,
                   color=CHANGE_COLORS.get(kind, MUTED), alpha=0.85,
                   linewidth=0.6, edgecolor=BG, zorder=3)

    for name, (fx, fy) in LENGTH_CALLOUTS.items():
        hit = d[d["alt_transcript_name"] == name]
        if hit.empty:
            continue
        r = hit.iloc[0]
        ax.annotate(name, xy=(r["can_aa_len"], r["alt_aa_len"]),
                    xytext=(r["can_aa_len"] * fx, r["alt_aa_len"] * fy),
                    fontsize=8.5, color=FG, va="center",
                    ha="left" if fx >= 1 else "right", zorder=5,
                    arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.7,
                                    shrinkA=0, shrinkB=3))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("canonical protein length  (aa)", fontsize=10.5, color=FG)
    ax.set_ylabel("alt protein length  (aa)", fontsize=10.5, color=FG)
    ax.set_title(f"B   …and whether a protein is left  "
                 f"({n_below} of {len(d)} fall below gate B)",
                 fontsize=11.5, color=FG, fontweight="bold", loc="left", pad=12)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def plot(pairs: pd.DataFrame) -> None:
    d = pairs.copy()
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(16.0, 7.6), facecolor=BG,
        gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.20,
                         left=0.105, right=0.985, top=0.845, bottom=0.155),
    )
    _style(ax_a)
    _style(ax_b)

    order = _strip_panel(ax_a, d)
    _length_panel(ax_b, d, order)

    trunc = d["protein_change_type"].isin(TRUNCATION_LIKE)
    fig.text(0.012, 0.985,
             f"Protein-level consequence of the switch — {len(d)} trial-failure "
             f"protein pairs ({d['gene_name'].nunique()} genes)",
             fontsize=15.5, color=FG, fontweight="bold", ha="left", va="top")
    fig.text(0.012, 0.947,
             f"{int(trunc.sum())} truncation-like pairs change a median "
             f"{d.loc[trunc, 'true_changed_frac'].median() * 100:.0f}% of the "
             f"sequence; the other {int((~trunc).sum())} change a median "
             f"{d.loc[~trunc, 'true_changed_frac'].median() * 100:.1f}% and never "
             f"exceed {d.loc[~trunc, 'true_changed_frac'].max() * 100:.0f}%",
             fontsize=10, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=9,
               color=CHANGE_COLORS.get(k, MUTED), markeredgecolor=BG,
               label=k.replace("_", " "))
        for k in order
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=len(order), frameon=False, fontsize=9.5, handletextpad=0.4,
               columnspacing=1.5)

    save(fig, "09_protein_change.png")


def main() -> None:
    pairs = build_gate_matrix()
    plot(pairs)

    d = pairs
    print(f"\npairs={len(d)} genes={d.gene_name.nunique()}")
    print("changed fraction by class:")
    summary = (d.groupby("protein_change_type")["true_changed_frac"]
               .agg(["count", "median", "min", "max"]).sort_values("median",
                                                                   ascending=False))
    print((summary * [1, 100, 100, 100]).round(2).to_string())
    kept = d["alt_aa_len"] / d["can_aa_len"]
    print(f"\nbelow gate B length floor ({MIN_KEPT_FRAC:.0%} kept): "
          f"{int((kept < MIN_KEPT_FRAC).sum())}")
    mech = d[d[list(MECHANISM_GATES)].all(axis=1)]
    print(f"mechanism-only set changed fraction: "
          f"{(mech['true_changed_frac'] * 100).round(1).tolist()}")


if __name__ == "__main__":
    main()

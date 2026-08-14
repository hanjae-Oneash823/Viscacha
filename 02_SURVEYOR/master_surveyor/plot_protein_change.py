"""09 -- how much of the protein actually changes in the alt isoform.

One point per docking pair from gate_matrix.build_gate_matrix() (94 distinct
canonical/alt protein pairs over 55 trial_failure genes; the canonical's own
row is excluded, so every pair here is a real comparison).

Left panel  -- RAW COUNT of canonical residues altered, by change class. This
is `n_changed` from the edlib CIGAR walk, NOT the changed_aa_start..end
envelope, and NOT normalized by canonical length (that normalized version --
true_changed_frac -- is what this panel used to plot; see git history).

Swapping fraction for count breaks the clean class separation that the
fraction gave you: the two ranges now overlap substantially (truncation-like
pairs span 1-1369 residues changed, everything else spans 1-236, and 22 of
the 55 truncation-like pairs fall inside that overlap). That's not a bug --
it's the reason the fraction was used originally. A 14-residue N-truncation
of a small protein and a 236-residue internal indel in a big one land at
opposite ends of "how disruptive is this" despite the count comparison
suggesting otherwise; raw count conflates "how much changed" with "how big
was the protein to begin with," which is exactly what the fraction was
built to factor out. Read this panel for absolute scale (e.g. "will this
change register on a docking pocket sized in residues"), not for a clean
per-class ranking -- use the fraction (or protein_change_type directly) for
that.

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
from matplotlib.transforms import blended_transform_factory

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
    """Jittered per-pair strip of changed-residue COUNT, one row per change class."""
    order = (d.groupby("protein_change_type")["n_changed"]
             .median().sort_values(ascending=False).index.tolist())
    rng = np.random.default_rng(0)
    # x is in axes fraction, y in data coords -- keeps the "n=" labels pinned
    # just past the right edge regardless of the log-scale data range.
    row_label = blended_transform_factory(ax.transAxes, ax.transData)

    for i, kind in enumerate(order):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=BAND, zorder=0, linewidth=0)
        sub = d[d["protein_change_type"] == kind]
        vals = sub["n_changed"].to_numpy()
        names = sub["alt_transcript_name"].to_numpy()
        jitter = rng.uniform(-0.26, 0.26, len(vals))
        ax.scatter(vals, i + jitter, s=52, color=CHANGE_COLORS.get(kind, MUTED),
                   alpha=0.85, linewidth=0.6, edgecolor=BG, zorder=3)
        med = float(np.median(vals))
        ax.plot([med, med], [i - 0.36, i + 0.36], color=FG, linewidth=2.0,
                zorder=4, solid_capstyle="butt")
        ax.text(1.015, i, f"n={len(sub)}", fontsize=8.5, color=MUTED,
                va="center", ha="left", transform=row_label)

        # Name the extreme-end gene at each end of the row -- these are
        # exactly the points that create (or shrink) the cross-class overlap
        # called out below, so they're the ones worth identifying by name.
        if len(vals) > 1:
            lo_pos, hi_pos = int(np.argmin(vals)), int(np.argmax(vals))
            for pos, dx, dy in ((lo_pos, -8, -11), (hi_pos, 8, 11)):
                ax.annotate(
                    names[pos], xy=(vals[pos], i + jitter[pos]),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=7.3, color=FG, va="center",
                    ha="right" if dx < 0 else "left", zorder=5,
                    arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.6,
                                    alpha=0.7, shrinkA=0, shrinkB=3),
                )

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([k.replace("_", " ") for k in order], fontsize=10,
                       color=FG)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xscale("log")
    ax.set_xlim(0.7, d["n_changed"].max() * 1.5)
    ax.set_xlabel("canonical residues altered in the alt  (count)", fontsize=10.5,
                  color=FG)
    ax.set_title("A   how much of the sequence changes", fontsize=11.5,
                 color=FG, fontweight="bold", loc="left", pad=12)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=1)
    ax.yaxis.grid(False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)

    # Unlike the fraction, raw count does NOT cleanly separate the classes --
    # say so instead of drawing a split line that would misrepresent the data.
    other = d.loc[~d["protein_change_type"].isin(TRUNCATION_LIKE), "n_changed"]
    trunc = d.loc[d["protein_change_type"].isin(TRUNCATION_LIKE), "n_changed"]
    n_overlap = int((trunc <= other.max()).sum())
    ax.text(0.97, 0.03,
            f"counts overlap between classes: {n_overlap} of {len(trunc)} "
            f"truncation-like pairs change fewer residues than the largest\n"
            f"indel/insertion/extension pair ({int(other.max())}) — count "
            f"conflates magnitude with protein size, fraction does not",
            fontsize=8, color=MUTED, ha="right", va="bottom", linespacing=1.6,
            style="italic", zorder=5, transform=ax.transAxes)
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
             f"{int(d.loc[trunc, 'n_changed'].median())} residues (as a fraction "
             f"of protein length, a median {d.loc[trunc, 'true_changed_frac'].median() * 100:.0f}%); "
             f"the other {int((~trunc).sum())} change a median "
             f"{int(d.loc[~trunc, 'n_changed'].median())} residues "
             f"({d.loc[~trunc, 'true_changed_frac'].median() * 100:.1f}% of length)",
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

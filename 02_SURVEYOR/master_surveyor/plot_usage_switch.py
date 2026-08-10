"""08 -- the isoform switch itself: canonical loss vs alt usage level in AD.

One point per docking pair from gate_matrix.build_gate_matrix() (94 distinct
canonical/alt protein pairs over 55 trial_failure genes).

  x  PERCENTAGE POINTS of within-gene isoform usage the canonical transcript
     LOSES in AD (control minus AD). Every trial_failure hit is CT_enriched by
     construction, so this is positive for all 94 pairs -- it is a magnitude.
  y  the alt's ABSOLUTE usage level in AD (alt_usage_pct_AD), not the change.
     "How much of this gene's transcript pool is the alt isoform in AD?"

Because y is a level and gate A is defined on the CHANGE, gate A cannot be
drawn as a horizontal line here -- passers span 10.0-88.8% AD usage and
failures span 3.2-50.8%, so the two overlap. The change is instead drawn per
point as a stem running from the alt's control level up (or down) to its AD
level: the DOT is the level, the STEM is the change, and gate A is the stem
being at least +10 pp long. Marker colour carries the gate A verdict.

The 50% reference is the meaningful landmark on a level axis: above it the alt
is the majority isoform in AD, i.e. the "alternative" has become the dominant
species. Only 8 of 94 pairs get there.

Worth noting on this axis and invisible on a delta axis: 6 of the 13
mechanism-only pairs sit at exactly 0% in controls, so the alt is not a shifted
isoform but one that is absent from controls entirely and appears in AD.
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
    MECHANISM_GATES, MIN_USAGE_DELTA, build_gate_matrix,
)
from master_surveyor.plot_results import BG, FG, GRID, _neon, _style, save

BAR = _neon(0.30)
HILITE = "#C41B3A"
NEUTRAL = "#cccccc"
MUTED = "#6b6b6b"

THRESHOLD_PP = MIN_USAGE_DELTA * 100

# Draw each pair's control level as a small dot joined to its AD level by a
# stem. The y axis is the AD LEVEL either way; the stem is what makes gate A
# legible, since two dots at the same height can differ on it. Set False for
# plain level-only markers.
SHOW_CONTROL_STEMS = True


def _declutter(ys: np.ndarray, min_gap: float) -> np.ndarray:
    """Push overlapping label anchors apart, preserving order."""
    order = np.argsort(ys)
    out = ys.astype(float).copy()
    for prev, cur in zip(order, order[1:]):
        if out[cur] - out[prev] < min_gap:
            out[cur] = out[prev] + min_gap
    return out


def plot(pairs: pd.DataFrame) -> None:
    d = pairs.copy()
    d["can_drop"] = -d["delta_usage"] * 100      # CT_enriched => always positive
    d["alt_ad"] = d["alt_usage_pct_AD"] * 100    # LEVEL in AD, not the change
    d["alt_ct"] = d["alt_usage_pct_control"] * 100
    d["mech"] = d[list(MECHANISM_GATES)].all(axis=1)

    fig, ax = plt.subplots(figsize=(12.2, 8.8), facecolor=BG)
    _style(ax)

    xmax = d["can_drop"].max() * 1.10
    ymin, ymax = 0, d["alt_ad"].max() * 1.14

    # Above 50% the alt is the majority isoform in AD.
    ax.axhspan(50, ymax, color="#f5f6f8", zorder=0, linewidth=0)
    ax.axhline(50, color=MUTED, linestyle=(0, (5, 4)), linewidth=1.1, zorder=1)
    n_major = int((d["alt_ad"] > 50).sum())
    ax.text(xmax * 0.995, 51.5,
            f"alt is the majority isoform in AD  ({n_major} of {len(d)})",
            fontsize=9, color=MUTED, ha="right", va="bottom", zorder=1)

    groups = [
        (~d["A"], NEUTRAL, 34, 0.85, "below gate A"),
        (d["A"] & ~d["mech"], BAR, 58, 0.95, "passes gate A"),
        (d["mech"], HILITE, 92, 1.0, f"mechanism-only set ({'+'.join(MECHANISM_GATES)})"),
    ]
    # Stem = the change gate A is defined on (control level -> AD level); dot =
    # the level itself. Gate A cannot be a horizontal line on a level axis, so
    # it lives in the stem length and the marker colour.
    for mask, color, size, alpha, _ in groups:
        sub = d[mask]
        if SHOW_CONTROL_STEMS:
            ax.vlines(sub["can_drop"], sub["alt_ct"], sub["alt_ad"], color=color,
                      linewidth=1.1, alpha=alpha * 0.55, zorder=3)
            ax.scatter(sub["can_drop"], sub["alt_ct"], s=size * 0.30, color=color,
                       alpha=alpha * 0.75, linewidth=0, zorder=3)
        ax.scatter(sub["can_drop"], sub["alt_ad"], s=size, color=color,
                   alpha=alpha, linewidth=0.6, edgecolor=BG, zorder=4)

    # Label only the mechanism-only set -- labelling 94 points is unreadable.
    # Two label columns, split on the x midpoint: labels always run outward
    # toward the sparse edges, so leaders stay short and never cross the dense
    # centre. Declutter each column independently.
    lab = d[d["mech"]]
    mid = (lab["can_drop"].min() + lab["can_drop"].max()) / 2
    columns = [(lab[lab["can_drop"] <= mid], -1, "right"),
               (lab[lab["can_drop"] > mid], +1, "left")]
    for grp, direction, ha in columns:
        if grp.empty:
            continue
        grp = grp.sort_values("alt_ad")
        ly = _declutter(grp["alt_ad"].to_numpy(), min_gap=(ymax - ymin) * 0.052)
        for (_, r), yl in zip(grp.iterrows(), ly):
            ax.annotate(
                r["alt_transcript_name"],
                xy=(r["can_drop"], r["alt_ad"]),
                xytext=(r["can_drop"] + direction * xmax * 0.042, yl),
                fontsize=8.5, color=HILITE, va="center", ha=ha,
                fontweight="bold", zorder=5,
                arrowprops=dict(arrowstyle="-", color=HILITE, linewidth=0.7,
                                alpha=0.55, shrinkA=0, shrinkB=2),
            )

    ax.set_xlim(0, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("canonical isoform usage LOST in AD  (percentage points)",
                  fontsize=10.5, color=FG)
    ax.set_ylabel("alt isoform usage in AD  (% of the gene's transcript pool)",
                  fontsize=10.5, color=FG)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    n_pass = int(d["A"].sum())
    n_absent = int((d["alt_ct"] == 0).sum())
    fig.text(0.012, 0.985,
             f"Alt isoform usage in AD — {len(d)} trial-failure protein pairs "
             f"({d['gene_name'].nunique()} genes)",
             fontsize=15.5, color=FG, fontweight="bold", ha="left", va="top")
    fig.text(0.012, 0.947,
             f"dot = usage level in AD, stem = shift from the control level; "
             f"{n_pass} of {len(d)} alts clear gate A (stem ≥ +{THRESHOLD_PP:.0f} pp) "
             f"and {n_absent} are absent from controls entirely",
             fontsize=10, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=9, color=color,
               markeredgecolor=BG, label=label)
        for _, color, _, _, label in groups
    ]
    if SHOW_CONTROL_STEMS:
        handles.append(Line2D([], [], marker="o", linestyle="-", markersize=4.5,
                              color=MUTED, alpha=0.7,
                              label="small dot = control level, stem = shift"))
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=9.5,
              handletextpad=0.4, labelspacing=0.5)

    save(fig, "08_usage_switch_scatter.png")


def main() -> None:
    pairs = build_gate_matrix()
    plot(pairs)

    d = pairs.assign(can_drop=-pairs["delta_usage"] * 100,
                     alt_ad=pairs["alt_usage_pct_AD"] * 100,
                     alt_ct=pairs["alt_usage_pct_control"] * 100)
    print(f"\npairs={len(d)}  canonical drop {d.can_drop.min():.1f}–"
          f"{d.can_drop.max():.1f} pp  |  alt usage in AD {d.alt_ad.min():.1f}–"
          f"{d.alt_ad.max():.1f}%")
    print(f"clear gate A (≥ +{THRESHOLD_PP:.0f} pp): {int(d['A'].sum())}")
    print(f"alt is the majority isoform in AD (> 50%): {int((d.alt_ad > 50).sum())}")
    print(f"alt absent from controls (0%): {int((d.alt_ct == 0).sum())}")
    print(f"alt usage in AD, gate A passers {d[d.A].alt_ad.min():.1f}–"
          f"{d[d.A].alt_ad.max():.1f}% vs failures {d[~d.A].alt_ad.min():.1f}–"
          f"{d[~d.A].alt_ad.max():.1f}% (overlapping — level is not the gate)")


if __name__ == "__main__":
    main()

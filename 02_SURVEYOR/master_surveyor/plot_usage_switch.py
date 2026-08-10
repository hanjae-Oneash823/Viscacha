"""08 -- the isoform switch itself: canonical loss vs alt gain, per protein pair.

One point per docking pair from gate_matrix.build_gate_matrix() (94 distinct
canonical/alt protein pairs over 55 trial_failure genes).

Both axes are PERCENTAGE POINTS of within-gene isoform usage, AD minus control:
  x  how much usage the canonical transcript LOSES in AD. Every trial_failure
     hit is CT_enriched by construction, so this is positive for all 94 pairs
     and the axis carries no sign information -- it is a magnitude.
  y  how much usage this particular alt GAINS. Free to be negative: most alts
     of a losing canonical lose ground too.

The y = x diagonal is the compositional ceiling for a two-isoform gene: usage
sums to 1 within a gene, so an alt sitting on the diagonal has absorbed the
canonical's entire loss by itself. Points below it mean the loss was shared
out across other isoforms. Points ABOVE it (2 of 94) are not an error -- they
need a third isoform to also be losing.

That constraint is why gate A is far more selective than it looks: at most one
alt per gene can absorb a large share, so raising the bar from 0 to +10 pp
cuts 94 pairs to 40.
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
    d["alt_gain"] = d["alt_usage_delta"] * 100
    d["mech"] = d[list(MECHANISM_GATES)].all(axis=1)

    fig, ax = plt.subplots(figsize=(12.2, 8.8), facecolor=BG)
    _style(ax)

    xmax = d["can_drop"].max() * 1.10
    ymin = min(d["alt_gain"].min() - 5, -12)
    ymax = d["alt_gain"].max() * 1.18

    # Region excluded by gate A, drawn first so every mark sits on top of it.
    ax.axhspan(ymin, THRESHOLD_PP, color="#f5f6f8", zorder=0, linewidth=0)

    # Compositional ceiling.
    lim = max(xmax, ymax)
    ax.plot([0, lim], [0, lim], color=MUTED, linestyle=(0, (5, 4)), linewidth=1.1,
            zorder=1)
    # Anchor the caption ON the diagonal, inside the y range -- the diagonal
    # runs off the top of the axes well before it reaches xmax.
    anchor = ymax * 0.80
    ax.text(anchor + xmax * 0.012, anchor, "alt absorbs the entire canonical loss",
            fontsize=8.5, color=MUTED, ha="left", va="bottom", zorder=1)

    ax.axhline(0, color=GRID, linewidth=1.0, zorder=1)
    ax.axhline(THRESHOLD_PP, color=HILITE, linestyle=(0, (6, 3)), linewidth=1.4,
               zorder=2)
    ax.text(xmax * 0.995, THRESHOLD_PP + 1.2,
            f"gate A: alt gains ≥ +{THRESHOLD_PP:.0f} pp",
            fontsize=9.5, color=HILITE, ha="right", va="bottom", fontweight="bold",
            zorder=3)

    groups = [
        (~d["A"], NEUTRAL, 34, 0.85, "below gate A"),
        (d["A"] & ~d["mech"], BAR, 58, 0.95, "passes gate A"),
        (d["mech"], HILITE, 92, 1.0, f"mechanism-only set ({'+'.join(MECHANISM_GATES)})"),
    ]
    for mask, color, size, alpha, _ in groups:
        sub = d[mask]
        ax.scatter(sub["can_drop"], sub["alt_gain"], s=size, color=color,
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
        grp = grp.sort_values("alt_gain")
        ly = _declutter(grp["alt_gain"].to_numpy(), min_gap=(ymax - ymin) * 0.058)
        for (_, r), yl in zip(grp.iterrows(), ly):
            ax.annotate(
                r["alt_transcript_name"],
                xy=(r["can_drop"], r["alt_gain"]),
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
    ax.set_ylabel("alt isoform usage GAINED in AD  (percentage points)",
                  fontsize=10.5, color=FG)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    n_pass = int(d["A"].sum())
    fig.text(0.012, 0.985,
             f"Isoform switch magnitude — {len(d)} trial-failure protein pairs "
             f"({d['gene_name'].nunique()} genes)",
             fontsize=15.5, color=FG, fontweight="bold", ha="left", va="top")
    fig.text(0.012, 0.947,
             f"{n_pass} of {len(d)} alts clear gate A; "
             f"{int((d['alt_gain'] < 0).sum())} lose usage alongside their canonical",
             fontsize=10, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=9, color=color,
               markeredgecolor=BG, label=label)
        for _, color, _, _, label in groups
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=9.5,
              handletextpad=0.4, labelspacing=0.5)

    save(fig, "08_usage_switch_scatter.png")


def main() -> None:
    pairs = build_gate_matrix()
    plot(pairs)

    d = pairs.assign(can_drop=-pairs["delta_usage"] * 100,
                     alt_gain=pairs["alt_usage_delta"] * 100)
    print(f"\npairs={len(d)}  canonical drop {d.can_drop.min():.1f}–"
          f"{d.can_drop.max():.1f} pp  |  alt gain {d.alt_gain.min():.1f}–"
          f"{d.alt_gain.max():.1f} pp")
    print(f"clear gate A (≥ +{THRESHOLD_PP:.0f} pp): {int(d['A'].sum())}")
    print(f"alt loses usage too (< 0): {int((d.alt_gain < 0).sum())}")
    print(f"above the compositional diagonal: {int((d.alt_gain > d.can_drop).sum())}")


if __name__ == "__main__":
    main()

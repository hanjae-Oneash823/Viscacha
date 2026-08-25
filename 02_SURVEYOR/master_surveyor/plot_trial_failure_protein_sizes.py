"""Protein-size distribution for the trial-failure candidate category.

Uses the uncollapsed gate-matrix universe: one row per distinct canonical / alt
protein pair.  Thus a gene with several genuinely different alternative
proteins contributes each distinct protein comparison, while duplicate
cell-type observations do not inflate the distribution.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from master_surveyor.gate_matrix import build_gate_matrix
from master_surveyor.plot_results import BG, FG, GRID, _style, save

CANONICAL_COLOR = "#2455A4"
ALT_COLOR = "#D05A2A"
MUTED = "#6B6B6B"


def _ecdf(values: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(values.to_numpy(dtype=float))
    return x, np.arange(1, len(x) + 1) / len(x)


def plot(pairs: pd.DataFrame) -> None:
    """Save a distribution and paired-size view for trial-failure proteins."""
    d = pairs[["gene_name", "alt_transcript_name", "can_aa_len", "alt_aa_len"]].copy()
    canonical = d["can_aa_len"]
    alt = d["alt_aa_len"]
    lo = max(1, int(min(canonical.min(), alt.min()) * 0.8))
    hi = int(max(canonical.max(), alt.max()) * 1.25)
    bins = np.geomspace(lo, hi, 22)

    fig, (ax_dist, ax_pairs) = plt.subplots(
        1, 2, figsize=(13.0, 6.4), facecolor=BG,
        gridspec_kw={"width_ratios": [1.05, 0.95], "wspace": 0.25},
    )
    _style(ax_dist)
    _style(ax_pairs)

    # Equal-width bins on a log axis show the full range without allowing the
    # handful of multi-thousand-residue proteins to obscure smaller proteins.
    ax_dist.hist(canonical, bins=bins, color=CANONICAL_COLOR, alpha=0.62,
                 label="canonical", edgecolor=BG, linewidth=0.5)
    ax_dist.hist(alt, bins=bins, color=ALT_COLOR, alpha=0.62,
                 label="alternative", edgecolor=BG, linewidth=0.5)
    for vals, color in ((canonical, CANONICAL_COLOR), (alt, ALT_COLOR)):
        x, y = _ecdf(vals)
        ax_dist.step(x, y * max(ax_dist.get_ylim()[1], 1), where="post",
                     color=color, linewidth=1.7)
    ax_dist.set_xscale("log")
    ax_dist.set_xlim(lo, hi)
    ax_dist.set_xlabel("protein length (aa; log scale)", fontsize=10.5, color=FG)
    ax_dist.set_ylabel("distinct protein pairs", fontsize=10.5, color=FG)
    ax_dist.set_title("A   size distribution", loc="left", fontsize=11.5,
                      fontweight="bold", color=FG, pad=10)
    ax_dist.legend(frameon=False, fontsize=9.5, loc="upper left")
    ax_dist.xaxis.grid(True, color=GRID, linewidth=0.7)
    ax_dist.yaxis.grid(True, color=GRID, linewidth=0.7)
    for side in ("top", "right"):
        ax_dist.spines[side].set_visible(False)

    # Each faint connector is a real canonical/alternative pair.  The box
    # summary makes the typical size immediately legible without hiding pairs.
    rng = np.random.default_rng(7)
    y_can = rng.normal(0, 0.055, len(d))
    y_alt = 1 + rng.normal(0, 0.055, len(d))
    for left, right, yl, yr in zip(canonical, alt, y_can, y_alt):
        ax_pairs.plot([left, right], [yl, yr], color="#B9BDC5", alpha=0.38,
                      linewidth=0.65, zorder=1)
    ax_pairs.scatter(canonical, y_can, s=24, color=CANONICAL_COLOR,
                     edgecolor=BG, linewidth=0.4, alpha=0.85, zorder=2)
    ax_pairs.scatter(alt, y_alt, s=24, color=ALT_COLOR,
                     edgecolor=BG, linewidth=0.4, alpha=0.85, zorder=2)
    ax_pairs.set_xscale("log")
    ax_pairs.set_xlim(lo, hi)
    ax_pairs.set_yticks([0, 1], ["canonical", "alternative"])
    ax_pairs.set_xlabel("protein length (aa; log scale)", fontsize=10.5, color=FG)
    ax_pairs.set_title("B   paired protein sizes", loc="left", fontsize=11.5,
                       fontweight="bold", color=FG, pad=10)
    ax_pairs.xaxis.grid(True, color=GRID, linewidth=0.7)
    ax_pairs.yaxis.grid(False)
    for side in ("top", "right", "left"):
        ax_pairs.spines[side].set_visible(False)

    summary = (
        f"canonical median {canonical.median():,.0f} aa; "
        f"alternative median {alt.median():,.0f} aa"
    )
    fig.text(0.012, 0.975,
             f"Protein sizes in trial-failure candidates — {len(d)} distinct "
             f"canonical/alternative pairs across {d['gene_name'].nunique()} genes",
             fontsize=14.5, fontweight="bold", color=FG, ha="left", va="top")
    fig.text(0.012, 0.935, summary, fontsize=9.5, color=MUTED,
             ha="left", va="top")
    fig.subplots_adjust(top=0.84, left=0.08, right=0.98, bottom=0.13)
    save(fig, "10d_trial_failure_protein_sizes.png")


def main() -> None:
    pairs = build_gate_matrix()
    plot(pairs)
    print(
        f"pairs={len(pairs)} genes={pairs.gene_name.nunique()} | "
        f"canonical median={pairs.can_aa_len.median():.0f} aa | "
        f"alternative median={pairs.alt_aa_len.median():.0f} aa"
    )


if __name__ == "__main__":
    main()

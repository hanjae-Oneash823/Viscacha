"""10 -- how many alt isoforms each trial_failure gene contributes.

The docking unit is a distinct (canonical, alt) PROTEIN pair, so a gene with
three alts that encode two distinct proteins counts as two. Three genes are
affected: IFNAR2-204/-205 and PADI2-203/-204 each encode identical proteins,
and STAT3 has 5 alt transcripts collapsing to 4 proteins -- 97 alt transcripts
become 94 pairs.

This is the distribution m0_select.representative_row() flattens: it keeps one
alt per gene, so the 28 genes with more than one would each lose every alt but
the top-scoring one. That is why the gate analysis works off the uncollapsed
set -- see gate_matrix.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import pandas as pd

from master_surveyor.gate_matrix import build_gate_matrix
from master_surveyor.plot_results import BG, FG, GRID, _neon, _style, save

BAR = _neon(0.30)
MUTED = "#6b6b6b"


def plot(pairs: pd.DataFrame) -> pd.Series:
    per_gene = pairs.groupby("gene_name").size()
    dist = per_gene.value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(8.6, 6.0), facecolor=BG)
    _style(ax)

    ax.bar(dist.index, dist.values, width=0.62, color=BAR, linewidth=0, zorder=3)
    for n_alts, n_genes in dist.items():
        ax.text(n_alts, n_genes + dist.max() * 0.018, str(n_genes), ha="center",
                va="bottom", fontsize=11, color=FG, fontweight="bold", zorder=4)

    ax.set_xticks(dist.index)
    ax.set_xlim(dist.index.min() - 0.6, dist.index.max() + 0.6)
    ax.set_ylim(0, dist.max() * 1.16)
    ax.set_xlabel("alt protein pairs contributed by the gene", fontsize=10.5,
                  color=FG)
    ax.set_ylabel("genes", fontsize=10.5, color=FG)
    ax.tick_params(labelsize=10, colors=MUTED, length=0)
    ax.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.xaxis.grid(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    n_multi = int((per_gene > 1).sum())
    fig.text(0.012, 0.977,
             f"Alt isoforms per trial-failure gene — {int(per_gene.sum())} pairs "
             f"across {len(per_gene)} genes",
             fontsize=14, color=FG, fontweight="bold", ha="left", va="top")
    fig.text(0.012, 0.930,
             f"{n_multi} genes contribute more than one alt; collapsing to one "
             f"per gene would drop {int(per_gene.sum() - len(per_gene))} pairs",
             fontsize=9.5, color=MUTED, ha="left", va="top")
    fig.subplots_adjust(top=0.855, left=0.095, right=0.975, bottom=0.105)

    save(fig, "10_alts_per_gene.png")
    return per_gene


def main() -> None:
    per_gene = plot(build_gate_matrix())
    dist = per_gene.value_counts().sort_index()
    print(f"\npairs={int(per_gene.sum())} genes={len(per_gene)}")
    for n_alts, n_genes in dist.items():
        print(f"  {n_alts} alt pair(s): {n_genes} genes")
    print("genes with the most:",
          ", ".join(f"{g} ({n})" for g, n in per_gene.nlargest(4).items()))


if __name__ == "__main__":
    main()

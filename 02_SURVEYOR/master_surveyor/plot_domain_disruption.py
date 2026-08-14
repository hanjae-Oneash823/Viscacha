"""09b -- which domain gets hit, and how badly.

09 (plot_protein_change.py) answers "how much of the SEQUENCE changed" and
"is a protein left at all" -- both are magnitude/survival questions, blind to
function. This is the functional follow-up: of the 56 (of 94) pairs where a
changed residue actually lands in a Pfam domain, which domain, and
is it clipped at the edge or gutted?

One row per domain-hit pair, from gate_matrix.build_gate_matrix()
(`domain_disruption` / `top_domain` / `domains_touched_kept`, computed in
_domain_overlap via the real edlib changed-residue set intersected with Pfam
coordinates -- not the changed_aa_start..end envelope).

  x         domain_disruption: fraction of the TOP touched domain's own
            residues that changed (0 = edge nick, 1 = the domain is gone).
  marker    filled = Pfam no longer detects that domain accession anywhere
            in the alt; ring = Pfam still calls a (possibly degenerate) hit
            for it. These disagree more than you'd expect -- TUT7-201 is
            98% disrupted by residue count but Pfam still calls PAP_assoc,
            because HMM hits don't require full-length coverage.
  color     protein_change_type, same palette as 09 and plot_results.

38 pairs never touch an annotated domain at all (indels/insertions landing
in linkers, or genes with sparse Pfam coverage) and are out of scope here --
see 09 for the sequence-level view that covers all 94.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from master_surveyor.gate_matrix import build_gate_matrix
from master_surveyor.plot_results import BG, CHANGE_COLORS, FG, GRID, _style, save

MUTED = "#6b6b6b"
BAND = "#f5f6f8"


def _kept(row: pd.Series) -> bool:
    kept = row["domains_touched_kept"]
    return isinstance(kept, str) and row["top_domain"] in [
        d.strip() for d in kept.split(",")
    ]


def plot(pairs: pd.DataFrame) -> None:
    # Domain overlap is descriptive context, not an active gate.
    d = pairs[pairs["n_touched"] > 0].copy()
    d["kept"] = d.apply(_kept, axis=1)
    d["disr_pct"] = d["domain_disruption"] * 100
    d = d.sort_values("disr_pct", ascending=False).reset_index(drop=True)
    n = len(d)

    order = [k for k in CHANGE_COLORS if k in set(d["protein_change_type"])]

    # Fixed-inch header/footer bands regardless of row count, so the title
    # and legend keep the same size and spacing whether n is 20 or 90 --
    # only the plot body (proportional to n) grows.
    top_in, bottom_in, body_per_row_in = 0.95, 1.55, 0.215
    fig_h = body_per_row_in * n + top_in + bottom_in
    fig, ax = plt.subplots(figsize=(10.5, fig_h), facecolor=BG)
    _style(ax)

    y = range(n)
    for yi in y:
        if yi % 2 == 0:
            ax.axhspan(yi - 0.5, yi + 0.5, color=BAND, zorder=0, linewidth=0)

    colors = d["protein_change_type"].map(lambda k: CHANGE_COLORS.get(k, MUTED))
    ax.hlines(y, 0, d["disr_pct"], color=colors, linewidth=1.3, alpha=0.7, zorder=2)

    filled = ~d["kept"]
    ax.scatter(d.loc[filled, "disr_pct"], d.index[filled], s=62,
               color=colors[filled], edgecolor=BG, linewidth=0.6, zorder=3)
    ax.scatter(d.loc[~filled, "disr_pct"], d.index[~filled], s=62,
               facecolor=BG, edgecolor=colors[~filled], linewidth=1.6, zorder=3)

    for yi, row in d.iterrows():
        if row["n_touched"] > 1:
            ax.text(row["disr_pct"] + 2.2, yi, f"+{row['n_touched'] - 1}",
                    fontsize=7, color=MUTED, va="center", ha="left", zorder=4)

    labels = [f"{r.alt_transcript_name}  ·  {r.top_domain}" for r in d.itertuples()]
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8.5, color=FG)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(-2, 108)
    ax.set_xlabel("residues of the top touched domain that changed  (%)",
                  fontsize=10.5, color=FG, labelpad=10)
    ax.tick_params(labelsize=9, colors=MUTED, length=0)
    ax.xaxis.grid(True, color=GRID, lw=0.7, zorder=1)
    ax.yaxis.grid(False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)

    n_lost = int(filled.sum())
    fig.suptitle(
        f"Which domain gets hit — {n} of {len(pairs)} trial-failure pairs "
        f"that touch a Pfam domain",
        fontsize=15, color=FG, fontweight="bold", x=0.012, ha="left", y=0.995,
    )
    fig.text(0.012, 1 - 0.40 / fig_h,
             f"filled = Pfam no longer detects the domain in the alt; "
             f"ring = Pfam still calls a (possibly degenerate) hit  —  "
             f"{n_lost} lost, {n - n_lost} still detected  —  "
             f"+N = other domains also touched",
             fontsize=9.5, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=8,
               color=CHANGE_COLORS.get(k, MUTED), markeredgecolor=BG,
               label=k.replace("_", " "))
        for k in order
    ] + [
        Line2D([], [], marker="o", linestyle="", markersize=8, color=MUTED,
               markeredgecolor=BG, label="domain lost (filled)"),
        Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=BG,
               markeredgecolor=MUTED, markeredgewidth=1.6, label="domain still detected (ring)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.18 / fig_h),
               ncol=min(len(handles), 5), frameon=False, fontsize=9,
               handletextpad=0.4, columnspacing=1.4)

    fig.subplots_adjust(left=0.28, right=0.97,
                        top=1 - top_in / fig_h, bottom=bottom_in / fig_h)
    save(fig, "09b_domain_disruption.png")


def main() -> None:
    pairs = build_gate_matrix()
    plot(pairs)

    d = pairs[pairs["n_touched"] > 0]
    d = d.assign(kept=d.apply(_kept, axis=1))
    print(f"\ndomain-hit pairs: {len(d)} / {len(pairs)}")
    print(f"domain lost from Pfam detection: {int((~d['kept']).sum())}")
    print(f"domain still Pfam-detected despite the hit: {int(d['kept'].sum())}")
    print(f"median disruption: {d['domain_disruption'].median() * 100:.1f}%")
    print(f"multi-domain pairs (n_touched > 1): {int((d['n_touched'] > 1).sum())}")


if __name__ == "__main__":
    main()

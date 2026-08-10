"""07 -- five-gate overlap (UpSet-style) for trial_failure_candidate docking pairs.

Gate definitions, the uncollapsed pair substrate, and the rationale for every
audit fix live in gate_matrix.py -- read that first. This module only renders.

Bars are SUPERSET counts: each column is "pairs passing AT LEAST the marked
gates", so columns overlap and do NOT sum to the total. Exclusive-intersection
counts (the usual UpSet convention) answer a different question and are much
harder to read off when the point is a selection funnel.

Runs the gate matrix itself, so there is no run-this-first step. Needs edlib
via gate_matrix -- use the oneash_dtu env.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from master_surveyor.config import OUT_DIR
from master_surveyor.gate_matrix import (
    GATE_LABEL, GATES, MECHANISM_GATES, build_gate_matrix, funnel,
)
from master_surveyor.plot_results import BG, FG, GRID, _neon, _style, save

# House design system (plot_results): neon purple for marks, the CHANGE_COLORS
# C_truncation crimson as the single accent, "identical" grey as the neutral.
BAR = _neon(0.30)
HILITE = "#C41B3A"
NEUTRAL = "#cccccc"
MUTED = "#6b6b6b"
DOT_OFF = "#dedede"
BAND = "#f5f6f8"


def _tint(hex_color: str, amount: float) -> str:
    """Blend toward white -- a secondary that stays on the house hue."""
    r, g, b = mcolors.to_rgb(hex_color)
    return mcolors.to_hex(tuple(c + (1.0 - c) * amount for c in (r, g, b)))


SET_BAR = _tint(BAR, 0.62)


def superset_counts(gates: pd.DataFrame) -> pd.DataFrame:
    """Every non-empty gate subset, counted as 'passes AT LEAST these gates'."""
    recs = []
    for r in range(1, len(GATES) + 1):
        for sub in itertools.combinations(GATES, r):
            hit = gates[list(sub)].all(axis=1)
            recs.append({
                "key": "".join("1" if g in sub else "0" for g in GATES),
                "gates": "".join(sub),
                "n": int(hit.sum()),
                "ngates": r,
                "genes": sorted(set(gates.loc[hit, "gene_name"])),
            })
    return (pd.DataFrame(recs).set_index("key")
            .sort_values(["n", "ngates"], ascending=[False, True]))


def plot_upset(gates: pd.DataFrame, combos: pd.DataFrame) -> None:
    n_pairs, n_genes = len(gates), gates["gene_name"].nunique()
    k = len(combos)
    all_key = "1" * len(GATES)
    mech_key = "".join("1" if g in MECHANISM_GATES else "0" for g in GATES)

    fig = plt.figure(figsize=(15.4, 9.2), facecolor=BG)
    gs = fig.add_gridspec(
        3, 2, width_ratios=[0.95, 4.5], height_ratios=[3.0, 2.5, 1.7],
        left=0.145, right=0.99, top=0.885, bottom=0.075, wspace=0.02, hspace=0.10,
    )
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_mat)
    ax_fun = fig.add_subplot(gs[2, :])
    for ax in (ax_bar, ax_mat, ax_set, ax_fun):
        _style(ax)

    x = np.arange(k)
    y = np.arange(len(GATES))
    top = combos["n"].max()

    # ---- panel 1: superset sizes ------------------------------------------
    ax_bar.bar(x, combos["n"], width=0.52, linewidth=0,
               color=[HILITE if key == all_key else BAR for key in combos.index])
    for xi, (key, n) in enumerate(zip(combos.index, combos["n"])):
        ax_bar.text(xi, n + top * 0.014, str(n), ha="center", va="bottom",
                    fontsize=7.5, color=HILITE if key == all_key else FG,
                    fontweight="bold" if key == all_key else "normal")
    ax_bar.set_ylabel("protein pairs passing\nAT LEAST these gates", fontsize=10,
                      color=FG, linespacing=1.5)
    ax_bar.set_ylim(0, top * 1.30)
    ax_bar.set_yticks([0, 30, 60, 90, 120])
    ax_bar.tick_params(axis="y", labelsize=8.5, colors=MUTED, length=0)
    ax_bar.tick_params(axis="x", length=0, labelbottom=False)
    ax_bar.yaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax_bar.xaxis.grid(False)
    for side in ("top", "right", "bottom"):
        ax_bar.spines[side].set_visible(False)

    ai = list(combos.index).index(all_key)
    n_all = int(combos.loc[all_key, "n"])
    ax_bar.annotate(
        f"all five gates — {n_all} pair{'s' if n_all != 1 else ''}\n"
        + (", ".join(combos.loc[all_key, "genes"]) or "(none)"),
        xy=(ai, n_all + top * 0.03), xytext=(ai - 2.5, top * 0.30),
        fontsize=9.5, color=HILITE, fontweight="bold", ha="right", va="center",
        linespacing=1.45,
        arrowprops=dict(arrowstyle="-", color=HILITE, linewidth=1.2,
                        shrinkA=0, shrinkB=3, connectionstyle="arc3,rad=0.2"),
    )

    # Outline the mechanism-only bar: neighbours share its height, so a leader
    # line alone does not identify which column it points at.
    mi = list(combos.index).index(mech_key)
    n_mech = int(combos.loc[mech_key, "n"])
    ax_bar.bar(mi, n_mech, width=0.52, facecolor="none", edgecolor=FG,
               linewidth=1.3, zorder=5)
    ax_bar.annotate(
        f"mechanism-only working set (no drug gate)\n"
        f"{'+'.join(MECHANISM_GATES)} — {n_mech} pairs",
        xy=(mi, n_mech + top * 0.02), xytext=(mi - 2.5, top * 0.55),
        fontsize=9.5, color=FG, ha="right", va="center", linespacing=1.45,
        arrowprops=dict(arrowstyle="-", color=FG, linewidth=1.1,
                        shrinkA=0, shrinkB=3, connectionstyle="arc3,rad=0.2"),
    )

    # ---- panel 2: membership matrix ---------------------------------------
    for yi in range(len(GATES)):
        if yi % 2 == 0:
            ax_mat.add_patch(plt.Rectangle((-0.7, yi - 0.5), k, 1.0, color=BAND,
                                           zorder=0, linewidth=0))
    for xi, key in enumerate(combos.index):
        members = [i for i, c in enumerate(key) if c == "1"]
        on = HILITE if key == all_key else FG
        ax_mat.scatter([xi] * len(GATES), y, s=26, color=DOT_OFF, zorder=2,
                       linewidth=0)
        ax_mat.plot([xi, xi], [min(members), max(members)], color=on,
                    linewidth=1.4, zorder=3, solid_capstyle="round")
        ax_mat.scatter([xi] * len(members), members, s=26, color=on, zorder=4,
                       linewidth=0)
    ax_mat.set_xlim(-0.7, k - 0.3)
    ax_mat.set_ylim(len(GATES) - 0.5, -0.5)
    ax_mat.set_yticks(y)
    ax_mat.set_yticklabels([])
    # No combo strings under the columns: the dots already encode them, and a
    # rotated label row steals the vertical space the left panel's axis label
    # and the funnel need. Read them off 07_gate_overlap_table.csv instead.
    ax_mat.tick_params(length=0, labelbottom=False)
    ax_mat.xaxis.grid(False)
    ax_mat.yaxis.grid(False)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)

    # ---- panel 3: per-gate totals -----------------------------------------
    totals = [int(gates[g].sum()) for g in GATES]
    ax_set.barh(y, totals, height=0.38, color=SET_BAR, linewidth=0)
    for yi, t in enumerate(totals):
        ax_set.text(t + n_pairs * 0.025, yi, str(t), va="center", ha="left",
                    fontsize=9, color=FG, fontweight="bold")
    ax_set.set_xlim(n_pairs * 1.16, 0)
    ax_set.invert_xaxis()
    # No axis label: it renders under the funnel row and gets swallowed. Each
    # bar is value-labelled and the total is in the title.
    ax_set.set_xticks([0, 50, 100])
    ax_set.tick_params(axis="x", labelsize=8.5, colors=MUTED, length=0)
    ax_set.tick_params(axis="y", length=0)
    ax_set.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax_set.yaxis.grid(False)
    for spine in ax_set.spines.values():
        spine.set_visible(False)

    # Gate row labels ride on ax_set.transAxes, NOT ax_mat.transData -- drawn
    # in data space at negative x they land under the sibling axes and vanish.
    for yi, g in enumerate(GATES):
        letter, name, detail = GATE_LABEL[g]
        yf = 1.0 - (yi + 0.5) / len(GATES)
        ax_set.text(-0.055, yf + 0.035, f"{letter}   {name}",
                    transform=ax_set.transAxes, ha="right", va="center",
                    fontsize=9.5, color=FG, fontweight="bold", clip_on=False)
        ax_set.text(-0.055, yf - 0.045, detail, transform=ax_set.transAxes,
                    ha="right", va="center", fontsize=8, color=MUTED,
                    clip_on=False)

    # ---- panel 4: cumulative funnel ---------------------------------------
    labels, counts = funnel(gates)
    fy = np.arange(len(counts))
    ax_fun.barh(fy, counts, height=0.42, linewidth=0,
                color=[NEUTRAL] + [HILITE if n <= counts[-1] else BAR
                                   for n in counts[1:]])
    for yi, n in enumerate(counts):
        ax_fun.text(n + n_pairs * 0.008, yi, str(n), va="center", ha="left",
                    fontsize=9.5, color=FG, fontweight="bold")
    ax_fun.set_yticks(fy)
    ax_fun.set_yticklabels(labels, fontsize=9, color=FG)
    ax_fun.invert_yaxis()
    ax_fun.set_xlim(0, n_pairs * 1.08)
    ax_fun.set_xlabel("protein pairs surviving cumulative gates", fontsize=9,
                      color=MUTED)
    ax_fun.tick_params(axis="x", labelsize=8.5, colors=MUTED, length=0)
    ax_fun.tick_params(axis="y", length=0)
    ax_fun.xaxis.grid(True, color=GRID, lw=0.7, zorder=0)
    ax_fun.yaxis.grid(False)
    for spine in ax_fun.spines.values():
        spine.set_visible(False)

    # Call out gate C specifically, NOT the steepest drop. Gate A is always the
    # largest cut simply because it goes first against the ungated total, which
    # says nothing. Gate C is the one that decides the shortlist, and it is a
    # literature-coverage filter rather than a property of the protein.
    di = labels.index("A–C")
    ax_fun.annotate(
        f"gate C ({GATE_LABEL['C'][1]}) cuts {counts[di - 1]} → {counts[di]}",
        xy=(counts[di] + n_pairs * 0.05, di), xytext=(n_pairs * 0.30, di + 0.15),
        fontsize=9.5, color=HILITE, va="center", ha="left",
        arrowprops=dict(arrowstyle="-", color=HILITE, linewidth=1.2,
                        shrinkA=0, shrinkB=4),
    )

    # ---- titles & legend ---------------------------------------------------
    fig.text(0.023, 0.962,
             f"Five selection gates over the {n_pairs} trial-failure protein "
             f"pairs ({n_genes} genes)",
             fontsize=16, color=FG, fontweight="bold", ha="left", va="top")
    fig.text(0.023, 0.928,
             f"bars = pairs passing AT LEAST the marked gates "
             f"({k} combinations, overlapping)",
             fontsize=10, color=MUTED, ha="left", va="top")

    handles = [
        Line2D([], [], marker="s", linestyle="", markersize=8, color=BAR,
               label="gate combination"),
        Line2D([], [], marker="s", linestyle="", markersize=8, color=HILITE,
               label="all five gates"),
        Line2D([], [], marker="s", linestyle="", markersize=8, color=SET_BAR,
               label="single-gate total (left panel)"),
    ]
    # Legend rides in the header band, top right. Inside the funnel panel it
    # lands on top of the full-width ungated bar.
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               ncol=3, frameon=False, fontsize=9.5, handletextpad=0.5,
               columnspacing=1.6)

    save(fig, "07_gate_overlap_upset.png")


def write_tables(gates: pd.DataFrame, combos: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    combo_path = OUT_DIR / "07_gate_overlap_table.csv"
    (combos.reset_index()
           .assign(genes=lambda d: [", ".join(v) for v in d["genes"]])
           .rename(columns={"n": "n_pairs_passing_at_least"})
           [["gates", "n_pairs_passing_at_least", "genes"]]
           .to_csv(combo_path, index=False))
    print(f"  saved → {combo_path.name}  ({len(combos)} combinations)")

    pair_path = OUT_DIR / "07_gate_pairs.csv"
    gates.to_csv(pair_path, index=False)
    print(f"  saved → {pair_path.name}  ({len(gates)} pairs)")


def main() -> None:
    gates = build_gate_matrix()
    combos = superset_counts(gates)
    plot_upset(gates, combos)
    write_tables(gates, combos)

    labels, counts = funnel(gates)
    print(f"\npairs={len(gates)} genes={gates.gene_name.nunique()} "
          f"| funnel={' → '.join(str(c) for c in counts)}")
    mech = gates[gates[list(MECHANISM_GATES)].all(axis=1)]
    print(f"{'+'.join(MECHANISM_GATES)} (mechanism-only):",
          ", ".join(sorted(set(mech.gene_name))) or "(none)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Diverging bar chart of the isoform_role (Dominant/Minor, empirical Control-PSI)
x tx_role_mane (Canonical/Alternate, MANE Select) 2x2 table, split by
usage_direction (up/down in AD) — the 3-way breakdown of the 2,514 MANE-covered
DTU hits computed in classify_hit_scenarios.py / classify_hit_scenarios_mane.py.

One row per (isoform_role, tx_role_mane) combo; bar extends right for
AD_enriched (up in AD) count, left for CT_enriched (down in AD) count.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python plot_dominance_x_mane.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.transforms import blended_transform_factory

plt.rcParams["font.family"] = "Liberation Sans"

REPO_ROOT = Path(__file__).resolve().parent
IN_CSV    = REPO_ROOT / "outputs/scenario_analysis/hit_scenarios_dominance_x_mane.csv"
OUT_DIR   = REPO_ROOT / "outputs/scenario_analysis/plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Reference palette (dataviz skill) — light mode, diverging blue/red pair.
COLOR_AD_ENRICHED = "#e34948"   # red  — up in AD
COLOR_CT_ENRICHED = "#2a78d6"   # blue — down in AD
INK_PRIMARY   = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED     = "#898781"
GRIDLINE      = "#e1e0d9"
BASELINE      = "#c3c2b7"
SURFACE       = "#fcfcfb"

ROW_ORDER = [
    ("Dominant", "Canonical"),
    ("Dominant", "Alternate"),
    ("Minor",    "Canonical"),
    ("Minor",    "Alternate"),
]


def main() -> None:
    df = pd.read_csv(IN_CSV)
    covered = df[df["tx_role_mane"] != "no_MANE_coverage"]

    rows = []
    for isoform_role, mane_role in ROW_ORDER:
        sub = covered[
            (covered["isoform_role"] == isoform_role)
            & (covered["tx_role_mane"] == mane_role)
        ]
        vc = sub["usage_direction"].value_counts()
        rows.append({
            "isoform_role": isoform_role,
            "label": mane_role,
            "AD_enriched": int(vc.get("AD_enriched", 0)),
            "CT_enriched": int(vc.get("CT_enriched", 0)),
        })
    plot_df = pd.DataFrame(rows)
    plot_df["total"] = plot_df["AD_enriched"] + plot_df["CT_enriched"]

    fig, ax = plt.subplots(figsize=(7.5, 4.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    y_pos = range(len(plot_df))[::-1]  # top row = first entry
    bar_height = 0.6

    ax.barh(y_pos, -plot_df["CT_enriched"], height=bar_height,
            color=COLOR_CT_ENRICHED, label="Down in AD (CT_enriched)", zorder=3)
    ax.barh(y_pos, plot_df["AD_enriched"], height=bar_height,
            color=COLOR_AD_ENRICHED, label="Up in AD (AD_enriched)", zorder=3)

    # Value labels at each segment end. Kept at a minimum standoff from the
    # center badge so they don't collide with it on rows with a tiny segment
    # (the badge's on-screen width is fixed regardless of the data scale).
    BADGE_CLEARANCE = 150
    for y, ct, ad, total in zip(y_pos, plot_df["CT_enriched"], plot_df["AD_enriched"], plot_df["total"]):
        ax.text(-max(ct + 15, BADGE_CLEARANCE), y, f"{ct}", ha="right", va="center",
                 fontsize=9.5, color=INK_PRIMARY)
        ax.text(max(ad + 15, BADGE_CLEARANCE), y, f"{ad}", ha="left", va="center",
                 fontsize=9.5, color=INK_PRIMARY)
        ax.text(0, y, f"n={total}", ha="center", va="center",
                 fontsize=8, color="white", fontweight="normal", zorder=5,
                 bbox=dict(boxstyle="round,pad=0.3", facecolor=INK_SECONDARY,
                           edgecolor="none", alpha=0.65))

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(plot_df["label"], fontsize=10.5, color=INK_PRIMARY)
    ax.tick_params(axis="y", pad=8)

    # Group headers ("Dominant" / "Minor") to the left of the per-row labels,
    # vertically centered across each group's two rows.
    group_trans = blended_transform_factory(ax.transAxes, ax.transData)
    for isoform_role in ("Dominant", "Minor"):
        y_vals = [y for y, role in zip(y_pos, plot_df["isoform_role"]) if role == isoform_role]
        y_center = sum(y_vals) / len(y_vals)
        ax.text(-0.20, y_center, isoform_role, transform=group_trans,
                 ha="right", va="center", fontsize=12, fontweight="bold",
                 color=INK_PRIMARY)

    max_val = max(plot_df["AD_enriched"].max(), plot_df["CT_enriched"].max())
    xlim = max_val * 1.28
    ax.set_xlim(-xlim, xlim)
    ax.axvline(0, color=BASELINE, linewidth=1.2, zorder=2)

    xticks = ax.get_xticks()
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{abs(int(t))}" for t in xticks], color=INK_MUTED, fontsize=9)
    ax.set_xlabel("Number of hits (n = 2,514 MANE-covered)", color=INK_SECONDARY, fontsize=10)

    ax.grid(axis="x", color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(BASELINE)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", color=INK_MUTED)

    # Group separator + labels between the Dominant and Minor row pairs.
    ax.axhline(1.5, color=GRIDLINE, linewidth=0.8, linestyle="--", zorder=1)

    ax.set_title(
        "Isoform role × transcript role, by AD usage direction",
        fontsize=13, color=INK_PRIMARY, loc="left", pad=28,
    )

    # Direction headers in place of a legend — colored to match the bars,
    # sitting directly above the halves of the plot they describe.
    ax.text(0.02, 1.06, "← Down in AD", transform=ax.transAxes,
             ha="left", va="bottom", fontsize=8.5, fontweight="bold",
             color=COLOR_CT_ENRICHED)
    ax.text(0.98, 1.06, "Up in AD →", transform=ax.transAxes,
             ha="right", va="bottom", fontsize=8.5, fontweight="bold",
             color=COLOR_AD_ENRICHED)

    ax.set_xlabel("Number of hits (n = 2,514 MANE-covered)", color=INK_SECONDARY, fontsize=10)
    fig.text(0.5, 0.015,
             "Isoform role = empirical Control-PSI dominance rank per (gene, cell type)  ·  "
             "Transcript role = MANE Select canonical vs. alternate",
             ha="center", va="bottom", fontsize=8, color=INK_MUTED)

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out_png = OUT_DIR / "dominance_x_mane_direction.png"
    out_pdf = OUT_DIR / "dominance_x_mane_direction.pdf"
    fig.savefig(out_png, dpi=300, facecolor=SURFACE)
    fig.savefig(out_pdf, facecolor=SURFACE)
    print(f"Saved -> {out_png}")
    print(f"Saved -> {out_pdf}")


if __name__ == "__main__":
    main()

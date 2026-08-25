#!/usr/bin/env python3
"""Create standalone, reusable docking plots and scientific schematics."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
OUT = CAMPAIGN / "standalone_figures"

INK = "#17212B"
SLATE = "#637485"
GRID = "#DCE3E8"
TEAL = "#00A896"
ORANGE = "#F4A261"
PURPLE = "#7146C7"
RED = "#E0524D"
BLUE = "#3478A8"


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.titlesize": 17,
        "axes.titleweight": "bold",
        "axes.labelsize": 12,
        "axes.labelcolor": INK,
        "axes.edgecolor": GRID,
        "axes.linewidth": 1.0,
        "xtick.color": SLATE,
        "ytick.color": SLATE,
        "text.color": INK,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def finish(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        kwargs = {"dpi": 320} if suffix == "png" else {}
        fig.savefig(
            OUT / f"{stem}.{suffix}",
            bbox_inches="tight",
            pad_inches=0.10,
            facecolor="white",
            **kwargs,
        )
    plt.close(fig)


def strip_spines(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)


def fyn_rmsd() -> None:
    fyn = CAMPAIGN / "FYN_saracatinib"
    summary = json.loads((fyn / "analysis" / "summary.json").read_text())
    rows = list(csv.DictReader((fyn / "analysis" / "all_poses.csv").open()))
    top = [r for r in rows if int(r["pose_rank"]) == 1]
    rmsd = np.asarray([float(r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"]) for r in top])
    order = np.argsort(rmsd)
    rmsd = rmsd[order]
    x = np.arange(1, len(rmsd) + 1)

    fig, ax = plt.subplots(figsize=(9.2, 5.6), facecolor="white")
    ax.vlines(x, 1.65, rmsd, color="#B8E5DF", lw=3, zorder=1)
    ax.scatter(x, rmsd, s=105, color=TEAL, edgecolor="white", linewidth=1.5, zorder=3)
    ax.axhline(2.0, color=RED, lw=2, ls=(0, (5, 4)), zorder=0)
    ax.text(9.42, 2.012, "validation threshold  2.0 Å", color=RED, ha="right", va="bottom", fontsize=10.5)
    ax.axhline(np.median(rmsd), color=INK, lw=1.4, alpha=.75)
    ax.text(9.42, np.median(rmsd) - .012, f"median  {np.median(rmsd):.2f} Å", color=INK, ha="right", va="top", fontsize=10.5)
    ax.set(xlim=(.55, 9.45), ylim=(1.65, 2.08), xlabel="Independent Vina run (ordered by RMSD)", ylabel="Top-pose RMSD to crystal ligand (Å)")
    ax.set_xticks(x)
    ax.set_yticks([1.7, 1.8, 1.9, 2.0])
    ax.set_title("FYN–saracatinib redocking is reproducible", loc="left", pad=16)
    ax.text(0, 1.01, f"{summary['top_pose_recovered_under_2A']}/{summary['completed_runs']} independent runs recovered the crystallographic pose below 2 Å", transform=ax.transAxes, color=SLATE, fontsize=11, va="bottom")
    ax.grid(axis="y", color=GRID, lw=.8, alpha=.7)
    strip_spines(ax)
    fig.tight_layout()
    finish(fig, "FYN_redocking_RMSD_replicates")


def fyn_domain_architecture() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 3.7), facecolor="white")
    ax.set_xlim(0, 540); ax.set_ylim(-.05, 1.25); ax.axis("off")
    ax.set_title("Predicted alternate FYN product lacks the validated kinase pocket", loc="left", pad=10)
    ax.text(0, .82, "Canonical FYN", fontsize=11.5, color=SLATE, va="center")
    ax.text(0, .31, "Alternate product", fontsize=11.5, color=SLATE, va="center")

    x0 = 105
    ax.add_patch(Rectangle((x0, .69), 415, .26, fc="#EDF1F4", ec="none"))
    domains = [
        (65, 145, "SH3", "#9BD8D0"),
        (150, 255, "SH2", "#ABC9DF"),
        (270, 520, "kinase", "#F3C69D"),
    ]
    for start, end, label, color in domains:
        left = x0 + (start / 537) * 415
        width = ((end - start) / 537) * 415
        ax.add_patch(FancyBboxPatch((left, .69), width, .26, boxstyle="round,pad=.008,rounding_size=.025", fc=color, ec="white", lw=1.2))
        ax.text(left + width / 2, .82, label, ha="center", va="center", fontsize=11, weight="bold")

    end_x = x0 + (115 / 537) * 415
    ax.add_patch(FancyBboxPatch((x0, .18), end_x - x0, .26, boxstyle="round,pad=.008,rounding_size=.025", fc=RED, ec="none", alpha=.88))
    ax.axvline(end_x, ymin=.22, ymax=.48, color=RED, lw=2)
    ax.text(end_x + 10, .31, "ends at residue 115", color=RED, fontsize=12, weight="bold", va="center")
    ax.text(x0, .02, "No complete SH3, SH2, or kinase domain is retained.", color=SLATE, fontsize=10.5)
    ax.text(520, .60, "537 aa", ha="right", color=SLATE, fontsize=9.5)
    fig.tight_layout()
    finish(fig, "FYN_domain_architecture")


def kit_scores() -> None:
    kit = CAMPAIGN / "KIT_masitinib"
    values = np.asarray([float(r["vina_affinity_kcal_mol"]) for r in csv.DictReader((kit / "analysis" / "masitinib_1T46_replicates.csv").open())])
    values = np.sort(values)
    rng = np.random.default_rng(715)
    jitter = rng.uniform(-.055, .055, len(values))

    fig, ax = plt.subplots(figsize=(9.2, 4.9), facecolor="white")
    ax.scatter(values, jitter, s=115, color=PURPLE, edgecolor="white", linewidth=1.5, zorder=3)
    ax.axvline(values.mean(), color=INK, lw=1.8)
    ax.axvspan(values.mean() - values.std(ddof=1), values.mean() + values.std(ddof=1), color=PURPLE, alpha=.10, lw=0)
    ax.text(values.mean(), .15, f"mean  {values.mean():.3f}", ha="center", va="bottom", fontsize=10.5, weight="bold")
    ax.set_xlim(-12.90, -12.72); ax.set_ylim(-.17, .22); ax.set_yticks([])
    ax.set_xlabel("Top Vina score (kcal/mol)")
    ax.set_title("Canonical KIT–masitinib scores are tightly reproducible", loc="left", pad=16)
    ax.text(0, 1.01, f"9 independent runs  •  SD {values.std(ddof=1):.3f} kcal/mol  •  range {values.min():.3f} to {values.max():.3f}", transform=ax.transAxes, color=SLATE, fontsize=11, va="bottom")
    ax.grid(axis="x", color=GRID, lw=.8, alpha=.7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    fig.tight_layout()
    finish(fig, "KIT_masitinib_score_replicates")


def kit_coverage() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 3.7), facecolor="white")
    ax.set_xlim(620, 835); ax.set_ylim(-.1, 1.15); ax.axis("off")
    ax.set_title("PDB 1T46 cannot support a direct KIT-223 comparison", loc="left", pad=10)
    y = .58
    ax.plot([630, 689], [y, y], color=BLUE, lw=17, solid_capstyle="round")
    ax.plot([762, 825], [y, y], color=BLUE, lw=17, solid_capstyle="round")
    ax.plot([689, 762], [y, y], color=GRID, lw=17, ls=(0, (2, 2)), solid_capstyle="butt")
    ax.text(659, .74, "resolved", color=BLUE, fontsize=10.5, ha="center")
    ax.text(793, .74, "resolved", color=BLUE, fontsize=10.5, ha="center")
    ax.text(725.5, .80, "unresolved in 1T46\nresidues 690–761", color=SLATE, fontsize=11, ha="center", linespacing=1.25)
    ax.axvline(715, ymin=.19, ymax=.69, color=RED, lw=2.4)
    ax.scatter([715], [y], s=95, color=RED, edgecolor="white", lw=1.2, zorder=4)
    ax.text(715, .17, "KIT-223 deletion site\nresidue 715", color=RED, fontsize=11.5, weight="bold", ha="center", va="top")
    ax.text(620, -.06, "A mutant score would require a separately rebuilt and validated loop model.", color=SLATE, fontsize=10.5)
    fig.tight_layout()
    finish(fig, "KIT_1T46_structure_coverage")


if __name__ == "__main__":
    fyn_rmsd()
    fyn_domain_architecture()
    kit_scores()
    kit_coverage()

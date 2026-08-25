#!/usr/bin/env python3
"""Create restrained, presentation-grade docking figures from real outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
ASSETS = CAMPAIGN / "presentation_v2" / "assets"
FINAL = CAMPAIGN / "presentation_v2" / "final"

BG = "#F7F8FA"
WHITE = "#FFFFFF"
INK = "#17212B"
SLATE = "#647381"
LIGHT = "#DDE3E8"
TEAL = "#00A896"
ORANGE = "#F4A261"
PURPLE = "#6F42C1"
RED = "#D9534F"
BLUE = "#2F6B9A"


def base_figure() -> plt.Figure:
    fig = plt.figure(figsize=(16, 9), dpi=120, facecolor=BG)
    fig.subplots_adjust(0, 0, 1, 1)
    return fig


def title(fig: plt.Figure, headline: str, subtitle: str) -> None:
    fig.text(.055, .925, headline, fontsize=27, weight="bold", color=INK, va="top")
    fig.text(.055, .865, subtitle, fontsize=12.5, color=SLATE, va="top")
    fig.add_artist(plt.Line2D([.055, .945], [.825, .825], color=LIGHT, lw=1.2, transform=fig.transFigure))


def cropped(path: Path) -> Image.Image:
    im = Image.open(path).convert("RGBA")
    alpha = np.asarray(im.getchannel("A"))
    ys, xs = np.where(alpha > 8)
    if len(xs):
        margin = 25
        box = (max(0, int(xs.min()) - margin), max(0, int(ys.min()) - margin), min(im.width, int(xs.max()) + margin), min(im.height, int(ys.max()) + margin))
        im = im.crop(box)
    return im


def image_axis(fig: plt.Figure, rect: list[float], image: Image.Image) -> plt.Axes:
    ax = fig.add_axes(rect)
    ax.imshow(image)
    ax.axis("off")
    return ax


def save(fig: plt.Figure, stem: str) -> None:
    FINAL.mkdir(parents=True, exist_ok=True)
    fig.savefig(FINAL / f"{stem}.png", dpi=240, facecolor=BG)
    fig.savefig(FINAL / f"{stem}.pdf", dpi=240, facecolor=BG)
    plt.close(fig)


def fyn_result_slide() -> None:
    fyn = CAMPAIGN / "FYN_saracatinib"
    summary = json.loads((fyn / "analysis" / "summary.json").read_text())
    rows = list(csv.DictReader((fyn / "analysis" / "all_poses.csv").open()))
    top = [r for r in rows if int(r["pose_rank"]) == 1]
    rmsd = np.asarray([float(r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"]) for r in top])

    fig = base_figure()
    title(fig, "The docking protocol reliably recovers saracatinib in FYN", "Crystal-structure control (PDB 10DJ) • nine independently seeded AutoDock Vina runs")
    image_axis(fig, [.045, .12, .52, .68], cropped(ASSETS / "fyn_overview_ray.png"))
    fig.text(.075, .135, "Protein cartoon", color=SLATE, fontsize=10)
    fig.text(.20, .135, "●", color=TEAL, fontsize=15, va="center")
    fig.text(.218, .135, "crystal pose", color=INK, fontsize=10)
    fig.text(.32, .135, "●", color=ORANGE, fontsize=15, va="center")
    fig.text(.338, .135, "best Vina pose", color=INK, fontsize=10)

    # Large metrics, deliberately unboxed.
    fig.text(.61, .735, "9/9", fontsize=35, weight="bold", color=TEAL)
    fig.text(.61, .695, "top poses below\nthe 2 Å threshold", fontsize=10.5, color=SLATE, linespacing=1.25)
    fig.text(.805, .735, f"{summary['top_pose_rmsd_median_A']:.2f} Å", fontsize=35, weight="bold", color=INK)
    fig.text(.805, .695, "median top-pose\nRMSD", fontsize=10.5, color=SLATE, linespacing=1.25)

    ax = fig.add_axes([.61, .48, .32, .16], facecolor="none")
    x = np.arange(1, len(rmsd) + 1)
    ax.scatter(x, rmsd, s=58, color=TEAL, edgecolor=WHITE, linewidth=.8, zorder=3)
    ax.axhline(2.0, color=RED, lw=1.5, ls=(0, (4, 3)))
    ax.text(9.15, 2.035, "2 Å validation threshold", color=RED, fontsize=9, va="bottom", ha="right")
    ax.set_xlim(.5, 9.5); ax.set_ylim(1.68, 2.08)
    ax.set_xticks(x); ax.set_xlabel("independent run", fontsize=9, color=SLATE)
    ax.set_ylabel("top-pose RMSD (Å)", fontsize=9, color=SLATE)
    ax.tick_params(colors=SLATE, labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(LIGHT)

    # Isoform/domain diagram.
    axd = fig.add_axes([.61, .17, .32, .22])
    axd.set_xlim(0, 537); axd.set_ylim(0, 1); axd.axis("off")
    axd.text(0, .93, "Structural interpretation of the alternate transcript", color=INK, fontsize=12, weight="bold")
    axd.text(0, .65, "canonical", color=SLATE, fontsize=9, va="center")
    axd.text(0, .28, "alternate", color=SLATE, fontsize=9, va="center")
    domains = [(65, 145, "SH3", "#9FD9D1"), (150, 255, "SH2", "#AFC9DF"), (270, 520, "kinase", "#F5C9A2")]
    axd.add_patch(Rectangle((60, .56), 470, .18, color="#E8ECEF", ec="none"))
    for start, end, label, color in domains:
        axd.add_patch(FancyBboxPatch((start, .55), end-start, .20, boxstyle="round,pad=.01,rounding_size=4", fc=color, ec="none"))
        axd.text((start+end)/2, .65, label, ha="center", va="center", fontsize=9, color=INK, weight="bold")
    axd.add_patch(Rectangle((60, .20), 55, .16, color=RED, alpha=.85, ec="none"))
    axd.axvline(115, ymin=.12, ymax=.49, color=RED, lw=2)
    axd.text(122, .28, "ends at residue 115", color=RED, fontsize=10, weight="bold", va="center")
    axd.text(60, .02, "No complete SH3, SH2, or kinase domain is retained.", color=SLATE, fontsize=9)

    fig.text(.055, .045, "Conclusion", fontsize=10, weight="bold", color=INK)
    fig.text(.122, .045, "The validated canonical pocket is absent from the predicted truncated product.", fontsize=10, color=INK)
    fig.text(.945, .045, "RMSD is fixed-frame and not symmetry-corrected.", fontsize=8.5, color=SLATE, ha="right")
    save(fig, "01_FYN_validated_result")


def fyn_pocket_slide() -> None:
    fig = base_figure()
    title(fig, "Crystal and docked saracatinib occupy the same FYN pocket", "Proper bond topology and ray-traced molecular rendering • best seeded Vina pose shown in orange")
    image_axis(fig, [.045, .10, .65, .70], cropped(ASSETS / "fyn_pocket_ray.png"))

    ax = fig.add_axes([.73, .18, .22, .54]); ax.axis("off")
    ax.text(0, .93, "Pose overlay", fontsize=16, weight="bold", color=INK, transform=ax.transAxes)
    ax.text(0, .80, "●", fontsize=20, color=TEAL, transform=ax.transAxes, va="center")
    ax.text(.09, .80, "crystal saracatinib", fontsize=11, color=INK, transform=ax.transAxes, va="center")
    ax.text(0, .70, "●", fontsize=20, color=ORANGE, transform=ax.transAxes, va="center")
    ax.text(.09, .70, "best Vina pose", fontsize=11, color=INK, transform=ax.transAxes, va="center")
    ax.plot([0, 1], [.60, .60], color=LIGHT, lw=1, transform=ax.transAxes)
    ax.text(0, .47, "1.75 Å", fontsize=34, weight="bold", color=TEAL, transform=ax.transAxes)
    ax.text(0, .40, "best top-pose RMSD", fontsize=11, color=SLATE, transform=ax.transAxes)
    ax.text(0, .24, "What this establishes", fontsize=12, weight="bold", color=INK, transform=ax.transAxes)
    ax.text(0, .07, "The receptor preparation, search box,\nand Vina settings can recover the\nknown canonical binding mode.", fontsize=10.5, color=SLATE, linespacing=1.45, transform=ax.transAxes)
    fig.text(.055, .045, "This validates pose recovery for canonical FYN; it does not measure cellular efficacy or prove disease causality.", fontsize=9.5, color=SLATE)
    save(fig, "02_FYN_pocket_overlay")


def kit_slide() -> None:
    kit = CAMPAIGN / "KIT_masitinib"
    values = np.asarray([float(r["vina_affinity_kcal_mol"]) for r in csv.DictReader((kit / "analysis" / "masitinib_1T46_replicates.csv").open())])
    fig = base_figure()
    title(fig, "KIT–masitinib provides a canonical baseline—not an isoform comparison", "Experimental c-KIT template PDB 1T46 • nine independently seeded canonical-pocket runs")
    image_axis(fig, [.045, .20, .40, .56], cropped(ASSETS / "kit_overview_ray.png"))
    fig.text(.075, .18, "canonical c-KIT cartoon", fontsize=9.5, color=SLATE)
    fig.text(.245, .18, "●", fontsize=15, color=PURPLE, va="center")
    fig.text(.263, .18, "docked masitinib", fontsize=9.5, color=INK)

    ax = fig.add_axes([.52, .59, .40, .14], facecolor="none")
    y = np.ones_like(values)
    ax.scatter(values, y, s=62, color=PURPLE, edgecolor=WHITE, lw=.8, zorder=3)
    ax.axvline(values.mean(), color=INK, lw=1.4)
    ax.text(values.mean(), 1.15, f"mean {values.mean():.2f}", ha="center", color=INK, fontsize=10, weight="bold")
    ax.set_xlim(-12.91, -12.70); ax.set_ylim(.75, 1.28); ax.set_yticks([])
    ax.set_xlabel("top Vina score (kcal/mol)", color=SLATE, fontsize=9)
    ax.tick_params(axis="x", colors=SLATE, labelsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False); ax.spines["bottom"].set_color(LIGHT)
    fig.text(.52, .755, "Canonical result is numerically reproducible", fontsize=13, weight="bold", color=INK)
    fig.text(.52, .72, f"9 runs • SD {values.std(ddof=1):.03f} kcal/mol • range {values.min():.3f} to {values.max():.3f}", fontsize=10, color=SLATE)

    axg = fig.add_axes([.52, .22, .40, .27]); axg.set_xlim(540, 950); axg.set_ylim(0, 1); axg.axis("off")
    axg.text(540, .92, "Why no KIT-223 docking score is reported", fontsize=13, weight="bold", color=INK)
    axg.plot([565, 689], [.51, .51], color=BLUE, lw=12, solid_capstyle="round")
    axg.plot([762, 933], [.51, .51], color=BLUE, lw=12, solid_capstyle="round")
    axg.plot([689, 762], [.51, .51], color=LIGHT, lw=12, ls=(0, (2, 2)), solid_capstyle="butt")
    axg.axvline(715, ymin=.24, ymax=.78, color=RED, lw=2)
    axg.text(715, .18, "KIT-223 deletion\nresidue 715", ha="center", color=RED, fontsize=10, weight="bold")
    axg.text(725, .68, "unresolved in 1T46\n(residues 690–761)", ha="center", color=SLATE, fontsize=9.5)
    axg.text(625, .66, "resolved", ha="center", color=BLUE, fontsize=9)
    axg.text(845, .66, "resolved", ha="center", color=BLUE, fontsize=9)
    axg.text(540, .02, "A comparative score would depend on an unvalidated rebuilt loop model.", color=SLATE, fontsize=9.5)

    fig.add_artist(FancyBboxPatch((.52, .075), .40, .085, boxstyle="round,pad=.01,rounding_size=.012", fc="#FDEDEC", ec="none", transform=fig.transFigure))
    fig.text(.54, .122, "No KIT-223 score", color=RED, fontsize=13, weight="bold", va="center")
    fig.text(.69, .122, "The experimental structure does not cover the altered residue.", color=INK, fontsize=10.5, va="center")
    fig.text(.055, .045, "Interpretation: stable canonical baseline only; no claim about isoform-dependent affinity.", fontsize=9.5, color=SLATE)
    save(fig, "03_KIT_canonical_baseline")


if __name__ == "__main__":
    fyn_result_slide()
    fyn_pocket_slide()
    kit_slide()

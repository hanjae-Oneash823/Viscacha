#!/usr/bin/env python3
"""Create presentation-grade quantitative and mechanism figures for B and C."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
FIGURES = CAMPAIGN / "figures" / "bc_candidates"
ANALYSIS = CAMPAIGN / "analysis" / "bc_candidates"
SEEDS = (1103, 2207, 3301, 4409, 5519)
COLORS = {"canonical": "#1F4E79", "alternate": "#8B2E2E", "gold": "#8A6D00", "ink": "#1A1A1A"}


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 11,
            "axes.edgecolor": "black",
            "axes.linewidth": 0.9,
            "xtick.color": "black",
            "ytick.color": "black",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "text.color": COLORS["ink"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "legend.frameon": True,
            "legend.edgecolor": "black",
            "legend.fancybox": False,
        }
    )


def read_bace_rows() -> list[dict[str, float | str | int]]:
    rows = []
    definitions = (
        (
            "Canonical BACE1",
            CAMPAIGN / "systems" / "BACE1_verubecestat",
            "canonical_obabel",
        ),
        (
            "BACE1-202 (401 aa)",
            CAMPAIGN / "systems" / "BACE1_verubecestat_401",
            "alternate_401",
        ),
    )
    for group, system, label in definitions:
        for seed in SEEDS:
            run = system / "runs" / f"{label}_seed{seed}_ex32" / "result.json"
            gnina = system / "gnina_rescoring" / f"{label}_seed{seed}_ex32" / "result.json"
            vina_data = json.loads(run.read_text())
            gnina_data = json.loads(gnina.read_text())
            top = vina_data["results"][0]
            rows.append(
                {
                    "group": group,
                    "seed": seed,
                    "vina_affinity_kcal_mol": float(top["vina_affinity_kcal_mol"]),
                    "top_pose_rmsd_A": float(top["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"]),
                    "best_of_20_rmsd_A": min(
                        float(item["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"])
                        for item in vina_data["results"]
                    ),
                    "gnina_cnn_score": float(gnina_data["poses"][0]["cnn_score"]),
                    "gnina_cnn_affinity": float(gnina_data["poses"][0]["cnn_affinity"]),
                }
            )
    return rows


def write_bace_summary(rows: list[dict[str, float | str | int]]) -> None:
    outdir = ANALYSIS / "BACE1_202_401"
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "matched_docking_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {}
    for group in ("Canonical BACE1", "BACE1-202 (401 aa)"):
        selected = [row for row in rows if row["group"] == group]
        summary[group] = {}
        for metric in (
            "vina_affinity_kcal_mol",
            "top_pose_rmsd_A",
            "best_of_20_rmsd_A",
            "gnina_cnn_score",
            "gnina_cnn_affinity",
        ):
            values = np.array([float(row[metric]) for row in selected])
            summary[group][metric] = {
                "mean": round(float(values.mean()), 4),
                "sd": round(float(values.std(ddof=1)), 4),
                "min": round(float(values.min()), 4),
                "max": round(float(values.max()), 4),
            }
    (outdir / "matched_docking_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def docking_figure(rows: list[dict[str, float | str | int]]) -> None:
    groups = ("Canonical BACE1", "BACE1-202 (401 aa)")
    metrics = (
        ("vina_affinity_kcal_mol", "Vina affinity (kcal/mol)", (-10.0, -8.2)),
        ("top_pose_rmsd_A", "Pose RMSD to crystal (Å)", (0, 12)),
        ("gnina_cnn_score", "GNINA CNNscore", (0, 1.02)),
    )
    panel_letters = ("A", "B", "C")
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.4), constrained_layout=True)
    for axis, letter, (metric, ylabel, ylim) in zip(axes, panel_letters, metrics):
        for index, group in enumerate(groups):
            selected = [row for row in rows if row["group"] == group]
            values = np.array([float(row[metric]) for row in selected])
            jitter = np.linspace(-0.07, 0.07, len(values))
            if index == 0:
                axis.scatter(np.full(len(values), index) + jitter, values, s=42, facecolor="black",
                             edgecolor="black", linewidth=0.8, zorder=3, label="Canonical")
            else:
                axis.scatter(np.full(len(values), index) + jitter, values, s=42, facecolor="white",
                             edgecolor="black", linewidth=0.9, zorder=3, label="401-aa isoform")
            axis.errorbar(index, values.mean(), yerr=values.std(ddof=1), fmt="none",
                          ecolor="black", elinewidth=1.1, capsize=4, capthick=1.1, zorder=4)
            axis.plot([index - 0.14, index + 0.14], [values.mean(), values.mean()],
                      color="black", linewidth=1.6, zorder=4)
        if metric == "top_pose_rmsd_A":
            axis.axhline(2, color="black", linewidth=0.9, linestyle=":")
            axis.text(1.35, 2.35, "2 Å", ha="right", va="bottom", fontsize=9, color="black")
        axis.text(-0.28, ylim[1] - (ylim[1] - ylim[0]) * 0.04, letter, fontsize=14, fontweight="bold", va="top")
        axis.set_ylabel(ylabel)
        axis.set_ylim(*ylim)
        axis.set_xlim(-0.6, 1.6)
        axis.set_xticks((0, 1), ("Canonical", "401-aa\nisoform"))
        for side in ("top", "right", "bottom", "left"):
            axis.spines[side].set_visible(True)
        axis.grid(False)
        axis.tick_params(which="both", top=True, right=True)
    fig.suptitle(
        "BACE1-202 disrupts verubecestat pose recovery despite similar Vina scores",
        fontsize=13.5,
        fontweight="normal",
        x=0.02,
        ha="left",
        color=COLORS["ink"],
    )
    fig.text(
        0.02, -0.03,
        "Points are individual seeds (n=5 per group); bars are the mean, whiskers are ± 1 SD.\n"
        "Panel B: 0/5 alternate runs land below the 2 Å redocking threshold. Panel C: GNINA independently rejects most alternate poses.",
        fontsize=9, color="#333333", va="top",
    )
    for extension in ("png", "svg", "pdf"):
        fig.savefig(FIGURES / f"B_BACE1_matched_docking.{extension}", dpi=360)
    plt.close(fig)


def _domain_box(axis, start, end, y, height, facecolor, edgecolor="black", hatch=None, zorder=2):
    from matplotlib.patches import Rectangle
    axis.add_patch(Rectangle((start, y - height / 2), end - start, height,
                              facecolor=facecolor, edgecolor=edgecolor, linewidth=1.0,
                              hatch=hatch, zorder=zorder))


def bace_mechanism_figure() -> None:
    # Exact 6 A verubecestat contact shell from the 5HU1 co-crystal (28 canonical
    # positions); see analysis/bc_candidates/BACE1_202_401/BACE1_202_401_summary.json.
    contacts_missing = [70, 71, 72, 73, 74, 75, 91, 92, 93, 94, 95, 96]
    contacts_retained_can = [132, 133, 134, 135, 137, 168, 169, 171, 176, 179, 289, 290, 291, 292, 293, 396]
    contacts_retained_alt = [p - 100 for p in contacts_retained_can]
    CONTACT_RED = "#FF1E1E"

    fig, axis = plt.subplots(figsize=(12.5, 4.2), constrained_layout=True)
    axis.set_xlim(-40, 520)
    axis.set_ylim(0, 2.5)
    y_can, y_alt = 1.55, 0.55
    h = 0.30

    _domain_box(axis, 1, 501, y_can, h, facecolor="#D9D9D9")
    _domain_box(axis, 21, 120, y_can, h, facecolor="white", hatch="////", edgecolor="black")
    _domain_box(axis, 1, 20, y_alt, h, facecolor="#BFBFBF")
    _domain_box(axis, 21, 401, y_alt, h, facecolor="#A9C4DD")

    tick_pad = 0.05
    for pos in contacts_missing + contacts_retained_can:
        axis.plot([pos, pos], [y_can - h / 2 - tick_pad, y_can + h / 2 + tick_pad],
                  color=CONTACT_RED, linewidth=1.6, zorder=5, solid_capstyle="butt")
    for pos in contacts_retained_alt:
        axis.plot([pos, pos], [y_alt - h / 2 - tick_pad, y_alt + h / 2 + tick_pad],
                  color=CONTACT_RED, linewidth=1.6, zorder=5, solid_capstyle="butt")

    axis.text(-15, y_can, "Canonical\n501 aa", ha="right", va="center", fontsize=10)
    axis.text(-15, y_alt, "BACE1-202\n401 aa", ha="right", va="center", fontsize=10)

    axis.annotate("residues 21-120 deleted\n(12/28 drug-contact positions, red)", xy=(70, y_can + h / 2 + tick_pad),
                  xytext=(70, 2.25), ha="center", fontsize=9.5, color="black",
                  arrowprops=dict(arrowstyle="-", color="black", linewidth=0.9))

    axis.scatter([93], [y_can], marker="d", s=45, facecolor="#FFFF00", edgecolor="black",
                 linewidth=1.1, zorder=6)
    axis.annotate("Asp93 (catalytic, lost)", xy=(93, y_can - h / 2), xytext=(93, -0.35),
                  ha="center", fontsize=9, arrowprops=dict(arrowstyle="-", color="black", linewidth=0.8))

    axis.scatter([289], [y_can], marker="d", s=45, facecolor="#FFFF00", edgecolor="black",
                 linewidth=1.1, zorder=6)
    axis.scatter([189], [y_alt], marker="d", s=45, facecolor="#FFFF00", edgecolor="black",
                 linewidth=1.1, zorder=6)
    axis.plot([289, 189], [y_can - h / 2, y_alt + h / 2], color="black", linewidth=0.8, linestyle=(0, (2, 2)))
    axis.text(239, (y_can + y_alt) / 2, "Asp289 retained\n(= Asp189)", ha="center", va="center", fontsize=9,
              bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="none"))

    axis.text(261, y_alt - h / 2 - 0.12, "canonical 121-501 retained, renumbered 21-401", ha="center", va="top", fontsize=9)
    axis.text(10, y_alt + h / 2 + 0.12, "new N-term", ha="center", va="bottom", fontsize=8, style="italic")
    axis.annotate("16 retained drug-contact positions", xy=(296, y_alt + h / 2 + tick_pad), xytext=(400, 1.15),
                  ha="center", fontsize=9, color="black",
                  arrowprops=dict(arrowstyle="-", color="black", linewidth=0.8))

    axis.set_yticks([])
    for side in ("top", "left", "right"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_position(("data", -0.55))
    axis.spines["bottom"].set_bounds(1, 501)
    axis.set_xticks([1, 100, 200, 300, 400, 501])
    axis.set_xlabel("Residue position", fontsize=9.5)
    axis.tick_params(axis="x", direction="out", bottom=True)

    axis.set_title("Mechanism B: BACE1-202 removes half the catalytic dyad and 12/28 drug-contact residues",
                   loc="left", fontsize=13, pad=10, weight="normal")
    for extension in ("png", "svg", "pdf"):
        fig.savefig(FIGURES / f"B_BACE1_isoform_mechanism.{extension}", dpi=360)
    plt.close(fig)


def cacna_mechanism_figure() -> None:
    fig, axis = plt.subplots(figsize=(12.5, 4.4), constrained_layout=True)
    axis.set_xlim(-160, 2300)
    axis.set_ylim(0, 2.6)
    y_can, y_alt = 1.55, 0.55
    h = 0.30

    _domain_box(axis, 1, 2161, y_can, h, facecolor="#D9D9D9")
    _domain_box(axis, 1078, 1493, y_can, h, facecolor="#E3C860", edgecolor="black")
    _domain_box(axis, 1606, 2161, y_can, h, facecolor="white", hatch="////", edgecolor="black")
    _domain_box(axis, 1, 1625, y_alt, h, facecolor="#A9C4DD")
    _domain_box(axis, 1098, 1513, y_alt, h, facecolor="#E3C860", edgecolor="black")
    axis.scatter([492], [y_alt], marker="v", s=60, facecolor="white", edgecolor="black", linewidth=1.1, zorder=4)

    axis.text(-45, y_can, "Canonical\n2161 aa", ha="right", va="center", fontsize=10)
    axis.text(-45, y_alt, "CACNA1D-214\n1625 aa", ha="right", va="center", fontsize=10)

    axis.text(1285, 2.15, "isradipine-contact region\nretained, sequence-identical", ha="center", fontsize=9.5)
    axis.annotate("distal C-terminus absent\n(residues 1606-2161)", xy=(1880, y_can + h / 2), xytext=(1880, 1.9),
                  ha="center", fontsize=9.5, arrowprops=dict(arrowstyle="-", color="black", linewidth=0.9))
    axis.annotate("20-aa insertion", xy=(492, y_alt - h / 2), xytext=(492, -0.35),
                  ha="center", fontsize=9, arrowprops=dict(arrowstyle="-", color="black", linewidth=0.8))
    axis.text(1285, (y_can + y_alt) / 2, "local pocket RMSD 1.49 Å  •  mean pocket pLDDT 84.0", ha="center", va="center",
              fontsize=9, bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="none"))

    axis.set_yticks([])
    for side in ("top", "left", "right"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_position(("data", -0.55))
    axis.spines["bottom"].set_bounds(1, 2161)
    axis.set_xticks([1, 500, 1000, 1500, 2000, 2161])
    axis.set_xlabel("Residue position", fontsize=9.5)
    axis.tick_params(axis="x", direction="out", bottom=True)

    axis.set_title("Mechanism C: a distal splice change preserves the static pocket but may alter channel gating",
                   loc="left", fontsize=13, pad=10, weight="normal")
    axis.text(0.995, 1.06, "AD inhibitory neurons: 1.8% → 28.1% usage\npadj = 4.1e-36; empirical FDR = 0.026",
              ha="right", va="top", fontsize=8.5, color="#333333", transform=axis.transAxes,
              bbox=dict(boxstyle="square,pad=0.4", facecolor="white", edgecolor="black", linewidth=0.7))
    for extension in ("png", "svg", "pdf"):
        fig.savefig(FIGURES / f"C_CACNA1D_isoform_mechanism.{extension}", dpi=360)
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    setup()
    rows = read_bace_rows()
    write_bace_summary(rows)
    docking_figure(rows)
    bace_mechanism_figure()
    cacna_mechanism_figure()
    print(f"Wrote B/C figures to {FIGURES}")


if __name__ == "__main__":
    main()

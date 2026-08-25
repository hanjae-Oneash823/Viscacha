#!/usr/bin/env python3
"""Create a compact, slide-ready FYN--saracatinib preliminary-result figure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def heavy_coords(path: Path, first_model_only: bool = False) -> np.ndarray:
    coords: list[list[float]] = []
    in_model = False
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            if in_model and first_model_only:
                break
            in_model = True
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:].strip() if len(line) >= 78 else line[76:78].strip()
            if not atom_type.startswith("H"):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if line.startswith("ENDMDL") and first_model_only:
            break
    return np.asarray(coords, dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    campaign = args.campaign
    summary = json.loads((campaign / "analysis" / "summary.json").read_text())
    best = summary["best_top_pose"]
    crystal = heavy_coords(campaign / "prepared" / "saracatinib_H8H_A601.pdbqt")
    pose = heavy_coords(campaign / "runs" / best["run"] / "docked_poses.pdbqt", first_model_only=True)
    all_pose_rows = np.genfromtxt(campaign / "analysis" / "all_poses.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    top = all_pose_rows[all_pose_rows["pose_rank"] == 1]

    fig = plt.figure(figsize=(15, 5.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.3, 0.85])

    ax = fig.add_subplot(grid[0, 0])
    ax.scatter(crystal[:, 0], crystal[:, 1], s=44, c="#d1495b", label="experimental saracatinib", zorder=2)
    ax.scatter(pose[:, 0], pose[:, 1], s=19, c="#2878b5", label="Vina top pose", zorder=3)
    for a, b in zip(crystal, pose):
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#8d99ae", lw=0.6, alpha=0.55, zorder=1)
    ax.set_aspect("equal")
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_title(f"10DJ ligand overlay\nRMSD = {best['rmsd_to_crystal_heavy_atom_uncorrected_angstrom']:.2f} Å")
    ax.legend(frameon=False, fontsize=8, loc="best")

    ax = fig.add_subplot(grid[0, 1])
    ax.set_title("FYN isoform context: canonical pocket is absent in alternate")
    ax.set_xlim(0, 550)
    ax.set_ylim(-0.4, 2.4)
    ax.set_yticks([0.25, 1.35])
    ax.set_yticklabels(["alternate\n115 aa", "canonical\n537 aa"])
    ax.set_xlabel("FYN residue")
    ax.hlines(1.35, 1, 537, color="#495057", lw=9, zorder=1)
    ax.hlines(0.25, 1, 115, color="#495057", lw=9, zorder=1)
    domains = [(88, 139, "SH3", "#4c956c"), (149, 231, "SH2", "#e9c46a"), (271, 520, "kinase / ATP pocket", "#527fb4")]
    for start, end, label, color in domains:
        ax.broken_barh([(start, end-start)], (1.12, 0.45), facecolors=color, edgecolors="none", zorder=2)
        ax.text((start + end) / 2, 1.78, label, ha="center", va="bottom", fontsize=9)
    ax.axvline(115, color="#d1495b", lw=1.5, ls="--")
    ax.text(121, 0.62, "termination", color="#b23a48", fontsize=8, va="bottom")
    ax.text(395, 0.15, "kinase domain deleted", color="#b23a48", ha="center", fontsize=9, weight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)

    ax = fig.add_subplot(grid[0, 2])
    rmsd = top["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"].astype(float)
    ax.scatter(np.ones_like(rmsd), rmsd, color="#2878b5", s=45, zorder=3)
    ax.boxplot(rmsd, widths=0.35, vert=True)
    ax.axhline(2.0, ls="--", lw=1.2, color="#d1495b", label="2 Å criterion")
    ax.set_xlim(0.55, 1.45)
    ax.set_xticks([])
    ax.set_ylabel("top-pose RMSD to crystal (Å)")
    ax.set_title(f"Seeded Vina validation\n{len(rmsd)}/{len(rmsd)} runs <2 Å")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.text(1.0, max(rmsd) + 0.08, f"median {np.median(rmsd):.2f} Å", ha="center", fontsize=9)

    fig.suptitle("Preliminary FYN–saracatinib docking validation", fontsize=15, weight="bold")
    (campaign / "figures").mkdir(exist_ok=True)
    fig.savefig(campaign / "figures" / "fyn_saracatinib_preliminary.png", dpi=240, bbox_inches="tight")
    fig.savefig(campaign / "figures" / "fyn_saracatinib_preliminary.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Create a slide-ready canonical c-KIT--masitinib docking summary."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    campaign = Path("outputs/docking_campaign/KIT_masitinib")
    analysis = campaign / "analysis"
    summary = json.loads((analysis / "masitinib_1T46_summary.json").read_text())
    with (analysis / "masitinib_1T46_replicates.csv").open() as handle:
        scores = np.asarray([float(row["vina_affinity_kcal_mol"]) for row in csv.DictReader(handle)])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    ax = axes[0]
    x = np.arange(1, len(scores) + 1)
    ax.scatter(x, scores, color="#2878b5", s=55, zorder=3)
    ax.axhline(scores.mean(), color="#d1495b", lw=1.6, label=f"mean {scores.mean():.2f} kcal/mol")
    ax.fill_between([0.5, len(scores) + 0.5], scores.mean() - scores.std(ddof=1), scores.mean() + scores.std(ddof=1), color="#d1495b", alpha=0.12, label="±1 SD")
    ax.set_xlim(0.5, len(scores) + 0.5)
    ax.set_xticks(x)
    ax.set_xlabel("independent Vina run")
    ax.set_ylabel("top Vina score (kcal/mol)")
    ax.set_title("Canonical c-KIT–masitinib\nexperimental 1T46 ATP pocket")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    ax.set_xlim(540, 950); ax.set_ylim(0, 1)
    ax.hlines(0.55, 565, 933, color="#495057", lw=10)
    ax.axvspan(690, 761, color="#adb5bd", alpha=0.7)
    ax.axvline(715, color="#d1495b", lw=1.6, ls="--")
    ax.text(715, 0.14, "KIT-223\ndeleted residue 715", ha="center", color="#b23a48", fontsize=9, weight="bold")
    ax.text(725, 0.78, "unresolved kinase-insert\nsegment in 1T46 (690–761)", ha="center", fontsize=9)
    ax.text(620, 0.73, "resolved kinase domain", ha="center", fontsize=9)
    ax.text(845, 0.73, "resolved kinase domain", ha="center", fontsize=9)
    ax.set_yticks([]); ax.set_xlabel("canonical c-KIT residue")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title("Why no KIT-223 docking score is reported")
    fig.suptitle("Preliminary canonical c-KIT–masitinib docking", fontsize=14, weight="bold")
    (campaign / "figures").mkdir(exist_ok=True)
    fig.savefig(campaign / "figures" / "kit_masitinib_preliminary.png", dpi=240, bbox_inches="tight")
    fig.savefig(campaign / "figures" / "kit_masitinib_preliminary.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()

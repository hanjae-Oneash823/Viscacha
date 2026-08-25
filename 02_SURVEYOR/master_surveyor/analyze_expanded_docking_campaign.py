#!/usr/bin/env python3
"""Aggregate, quality-label, report, and plot the all-candidate campaign."""

from __future__ import annotations

import csv
import json
import math
import statistics
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
SYSTEMS = CAMPAIGN / "systems"
ANALYSIS = CAMPAIGN / "analysis" / "aggregate"
FIGURES = CAMPAIGN / "figures" / "expanded_campaign"
ANALYSIS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


def read_results(paths: list[Path]) -> list[dict[str, object]]:
    rows = []
    seen = set()
    for path in paths:
        if path in seen or not path.exists():
            continue
        seen.add(path)
        data = json.loads(path.read_text())
        top = data["results"][0]
        rows.append({
            "path": str(path.relative_to(ROOT)),
            "seed": int(data["seed"]),
            "score": float(top["vina_affinity_kcal_mol"]),
            "rmsd": top.get("rmsd_to_crystal_heavy_atom_uncorrected_angstrom"),
        })
    return sorted(rows, key=lambda x: int(x["seed"]))


def group(candidate: str, exact: list[str], globs: list[str]) -> list[dict[str, object]]:
    base = SYSTEMS / candidate / "runs"
    paths = [base / p / "result.json" for p in exact]
    for pattern in globs:
        paths.extend(base.glob(f"{pattern}/result.json"))
    return read_results(paths)


GROUPS = {
    "BACE1 canonical": group("BACE1_verubecestat", ["canonical_obabel_redock_seed20260825_ex32"], ["canonical_obabel_seed*_ex32"]),
    "BACE1-476": group("BACE1_verubecestat", [], ["alternate_476_seed*_ex32"]),
    "BACE1-457": group("BACE1_verubecestat", [], ["alternate_457_seed*_ex32"]),
    "CHRNA7 canonical": group("CHRNA7_encenicline", ["canonical_obabel_redock_seed20260825_ex32"], ["canonical_obabel_seed*_ex32"]),
    "CHRFAM7A at A face": group("CHRNA7_encenicline", [], ["hybrid_A_face_seed*_ex32"]),
    "CHRFAM7A at B face": group("CHRNA7_encenicline", [], ["hybrid_B_face_seed*_ex32"]),
    "GABRA2 canonical": group("GABRA2_AZD7325", [], ["canonical_crossdock_seed*_ex32"]),
    "CACNA1D shared pocket": group("CACNA1D_isradipine", ["canonical_isradipine_seed20260825_ex32"], ["shared_pocket_crossdock_seed*_ex32"]),
    "PDE9A canonical": group("PDE9A_BI409306", [], ["canonical_crossdock_seed*_ex32"]),
}


def stats(rows: list[dict[str, object]]) -> dict[str, object]:
    scores = [float(x["score"]) for x in rows]
    rmsds = [float(x["rmsd"]) for x in rows if x["rmsd"] is not None]
    return {
        "n": len(scores),
        "mean_score": round(statistics.mean(scores), 3),
        "sd_score": round(statistics.stdev(scores), 3) if len(scores) > 1 else 0.0,
        "min_score": round(min(scores), 3),
        "max_score": round(max(scores), 3),
        "mean_top_pose_rmsd": round(statistics.mean(rmsds), 3) if rmsds else None,
        "sd_top_pose_rmsd": round(statistics.stdev(rmsds), 3) if len(rmsds) > 1 else (0.0 if rmsds else None),
        "top_pose_success_lt2A": sum(x < 2.0 for x in rmsds) if rmsds else None,
    }


SUMMARIES = {name: stats(rows) for name, rows in GROUPS.items()}


def paired_delta(reference: str, alternate: str) -> dict[str, object]:
    ref = {int(x["seed"]): float(x["score"]) for x in GROUPS[reference]}
    alt = {int(x["seed"]): float(x["score"]) for x in GROUPS[alternate]}
    seeds = sorted(set(ref) & set(alt))
    deltas = [alt[s] - ref[s] for s in seeds]
    return {
        "alternate_minus_canonical_kcal_mol": round(statistics.mean(deltas), 3),
        "sd": round(statistics.stdev(deltas), 3) if len(deltas) > 1 else 0.0,
        "paired_seeds": seeds,
        "interpretation": "Positive means a less favorable Vina score for the alternate under the matched protocol.",
    }


DELTAS = {
    "BACE1-476": paired_delta("BACE1 canonical", "BACE1-476"),
    "BACE1-457": paired_delta("BACE1 canonical", "BACE1-457"),
    "CHRFAM7A at A face": paired_delta("CHRNA7 canonical", "CHRFAM7A at A face"),
    "CHRFAM7A at B face": paired_delta("CHRNA7 canonical", "CHRFAM7A at B face"),
}


STATUS_ROWS = [
    {"priority": 1, "gene": "FYN", "drug": "saracatinib", "status": "validated + pocket absent", "numeric": "canonical only", "confidence": "high", "reason": "Canonical redocking validated; 115-aa alternate lacks the kinase domain."},
    {"priority": 2, "gene": "KIT", "drug": "masitinib", "status": "attempted; alternate unresolved", "numeric": "canonical baseline", "confidence": "limited", "reason": "Ser715 lies in the unresolved 1T46 segment; predicted-model docking showed clashes and was rejected."},
    {"priority": 3, "gene": "GABRA2", "drug": "AZD7325", "status": "canonical cross-dock + pocket absent", "numeric": "canonical only", "confidence": "exploratory", "reason": "9CSB alpha2/gamma2 pocket passed transfer checks; 73-aa alternate cannot form the receptor."},
    {"priority": 4, "gene": "CACNA1D", "drug": "isradipine", "status": "shared-pocket control", "numeric": "identical coordinates", "confidence": "exploratory", "reason": "The deposited human pocket ends before the alternate C-terminal truncation."},
    {"priority": 5, "gene": "BACE1-476", "drug": "verubecestat", "status": "modeled comparative docking", "numeric": "canonical vs alternate", "confidence": "moderate", "reason": "Canonical cognate redocking validated; alternate pocket has high local pLDDT."},
    {"priority": 6, "gene": "BACE1-457", "drug": "verubecestat", "status": "modeled comparative docking", "numeric": "canonical vs alternate", "confidence": "moderate", "reason": "Deletion removes 5 of 28 experimental pocket-contact residues; model passed local quality checks."},
    {"priority": 7, "gene": "CHRNA7/CHRFAM7A", "drug": "encenicline", "status": "two topology hypotheses", "numeric": "canonical vs two hybrids", "confidence": "hypothesis only", "reason": "Canonical redocking validated, but sample stoichiometry/genotype support is absent."},
    {"priority": 8, "gene": "PDE9A", "drug": "BI 409306", "status": "canonical cross-dock only", "numeric": "alternate undefined", "confidence": "incomplete pair", "reason": "The candidate file does not name the coding-altered PDE9A transcript."},
]


summary_json = {
    "cpu_policy": "All new CPU-bound jobs ran serially with Vina cpu=16; AlphaFold jobs were taskset-affined to CPUs 0-15.",
    "groups": SUMMARIES,
    "paired_deltas": DELTAS,
    "candidate_status": STATUS_ROWS,
    "noncomparability_warning": "Absolute Vina scores must not be compared across different proteins or drugs.",
}
(ANALYSIS / "expanded_summary.json").write_text(json.dumps(summary_json, indent=2) + "\n")

with (ANALYSIS / "candidate_status.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(STATUS_ROWS[0]), lineterminator="\n")
    writer.writeheader(); writer.writerows(STATUS_ROWS)

with (ANALYSIS / "replicate_results.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["group", "seed", "score", "rmsd", "path"],
        lineterminator="\n",
    )
    writer.writeheader()
    for name, rows in GROUPS.items():
        for row in rows:
            writer.writerow({"group": name, **row})


# Presentation-quality paired comparison panels.
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False,
})
palette = ["#213B5C", "#00A6A6", "#E07A5F"]
fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.7))

def dot_summary(ax, names, ylabel, title, metric):
    for i, name in enumerate(names):
        values = [float(r[metric]) for r in GROUPS[name] if r[metric] is not None]
        jitter = np.linspace(-0.11, 0.11, len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, s=52, color=palette[i], alpha=0.78,
                   edgecolor="white", linewidth=0.8, zorder=3)
        mean = np.mean(values); sd = np.std(values, ddof=1) if len(values) > 1 else 0
        ax.errorbar(i, mean, yerr=sd, fmt="D", color="#111827", mfc="white", mec="#111827",
                    markersize=7, capsize=5, linewidth=1.8, zorder=4)
    ax.set_xticks(range(len(names)), [n.replace(" canonical", "\ncanonical").replace("CHRFAM7A at ", "fusion at\n") for n in names])
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.7, alpha=0.65)

dot_summary(axes[0, 0], ["BACE1 canonical", "BACE1-476", "BACE1-457"], "Top Vina score (kcal/mol)", "A  BACE1–verubecestat score", "score")
dot_summary(axes[0, 1], ["BACE1 canonical", "BACE1-476", "BACE1-457"], "Top-pose RMSD to 5HU1 (Å)", "B  BACE1 pose recovery / shift", "rmsd")
axes[0, 1].axhline(2.0, color="#B91C1C", linestyle="--", linewidth=1.2)
axes[0, 1].text(2.35, 2.12, "2 Å validation line", color="#B91C1C", ha="right", fontsize=9)
dot_summary(axes[1, 0], ["CHRNA7 canonical", "CHRFAM7A at A face", "CHRFAM7A at B face"], "Top Vina score (kcal/mol)", "C  α7–encenicline score", "score")
dot_summary(axes[1, 1], ["CHRNA7 canonical", "CHRFAM7A at A face", "CHRFAM7A at B face"], "Top-pose RMSD to 7EKP (Å)", "D  α7 site pose displacement", "rmsd")
axes[1, 1].axhline(2.0, color="#B91C1C", linestyle="--", linewidth=1.2)
fig.suptitle("Matched docking comparisons: replicate-level stability and pose behavior", fontsize=16, fontweight="bold")
fig.tight_layout(rect=[0, 0.065, 1, 0.955], h_pad=2.2, w_pad=2.0)
fig.text(0.5, 0.022, "Diamonds show mean ± SD. Scores are comparable only within each protein–drug system.", ha="center", color="#4B5563", fontsize=10)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(FIGURES / f"matched_comparative_docking.{ext}", dpi=320, bbox_inches="tight")
plt.close(fig)


# Campaign decision/status graphic.
status_colors = {"high": "#128C7E", "moderate": "#2E86AB", "exploratory": "#E9A23B", "limited": "#D97706", "hypothesis only": "#9B5DE5", "incomplete pair": "#6B7280"}
fig, ax = plt.subplots(figsize=(13.5, 8.2))
ax.set_xlim(0, 10); ax.set_ylim(-0.2, 8.0); ax.axis("off")
for idx, row in enumerate(STATUS_ROWS):
    y = 6.85 - idx * 0.91
    color = status_colors[row["confidence"]]
    ax.add_patch(plt.Rectangle((0.15, y - 0.34), 0.12, 0.68, color=color, lw=0))
    ax.text(0.42, y + 0.12, f"{row['priority']}. {row['gene']}  ·  {row['drug']}", fontsize=12, fontweight="bold", va="center")
    ax.text(0.42, y - 0.16, row["status"], fontsize=9.8, color="#374151", va="center")
    ax.text(4.08, y + 0.02, textwrap.fill(row["reason"], width=58), fontsize=8.9, color="#111827", va="center", linespacing=1.15)
    ax.text(9.78, y + 0.02, row["confidence"].upper(), fontsize=8.6, color=color, fontweight="bold", ha="right", va="center")
    if idx < 7:
        ax.plot([0.15, 9.85], [y - 0.46, y - 0.46], color="#E5E7EB", lw=0.8)
ax.text(0.15, 7.72, "All eight candidate rows were attempted", fontsize=17, fontweight="bold", color="#111827")
ax.text(0.15, 7.39, "A numerical score was retained only when the receptor pocket and alternate definition supported it.", fontsize=10.5, color="#4B5563")
for ext in ["png", "pdf", "svg"]:
    fig.savefig(FIGURES / f"all_candidate_status.{ext}", dpi=320, bbox_inches="tight")
plt.close(fig)


# Stable canonical cross-docking scores, explicitly faceted to discourage comparison.
names = ["GABRA2 canonical", "CACNA1D shared pocket", "PDE9A canonical"]
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
for ax, name, color in zip(axes, names, ["#2E86AB", "#E07A5F", "#00A6A6"]):
    values = [float(r["score"]) for r in GROUPS[name]]
    ax.scatter(np.arange(1, len(values) + 1), values, s=58, color=color, edgecolor="white", linewidth=0.8)
    ax.axhline(np.mean(values), color="#111827", lw=1.5)
    ax.fill_between([0.6, len(values) + 0.4], np.mean(values) - np.std(values, ddof=1), np.mean(values) + np.std(values, ddof=1), color=color, alpha=0.15)
    ax.set_title(name.replace(" canonical", "\ncanonical"), fontsize=11)
    ax.set_xlabel("Independent seed"); ax.set_ylabel("Top Vina score (kcal/mol)")
    ax.grid(axis="y", color="#E5E7EB", lw=0.7)
fig.suptitle("Exploratory canonical cross-docking: within-system reproducibility", fontsize=15, fontweight="bold")
fig.text(0.5, -0.02, "Do not compare absolute scores between panels; each uses a different receptor–drug system.", ha="center", color="#4B5563", fontsize=9.5)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(FIGURES / f"canonical_crossdock_stability.{ext}", dpi=320, bbox_inches="tight")
plt.close(fig)


print(json.dumps(summary_json, indent=2))

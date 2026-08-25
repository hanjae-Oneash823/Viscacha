#!/usr/bin/env python3
"""Aggregate GNINA rescoring outputs and render the ggplot2 figure suite."""

from __future__ import annotations

import csv
import json
import math
import os
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
SYSTEMS = CAMPAIGN / "systems"
ANALYSIS = CAMPAIGN / "analysis" / "gnina"
FIGURES = CAMPAIGN / "figures" / "gnina_comparison"

GROUP_ORDER = [
    "FYN canonical",
    "BACE1 canonical",
    "BACE1-476",
    "BACE1-457",
    "CHRNA7 canonical",
    "CHRFAM7A at A face",
    "CHRFAM7A at B face",
]
CANONICAL_GROUPS = ["FYN canonical", "BACE1 canonical", "CHRNA7 canonical"]
COLORS = {
    "FYN canonical": "#274C77",
    "BACE1 canonical": "#1B4965",
    "BACE1-476": "#2A9D8F",
    "BACE1-457": "#E76F51",
    "CHRNA7 canonical": "#3A506B",
    "CHRFAM7A at A face": "#8E6CFF",
    "CHRFAM7A at B face": "#D1495B",
}


def mean(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def sample_sd(values: list[float]) -> float | None:
    return statistics.stdev(values) if len(values) > 1 else (0.0 if values else None)


def round_or_none(value: float | None, digits: int = 4) -> float | None:
    return round(value, digits) if value is not None else None


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 2 or len(x) != len(y):
        return None
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    x_centered = x_array - x_array.mean()
    y_centered = y_array - y_array.mean()
    denominator = math.sqrt(float(np.dot(x_centered, x_centered) * np.dot(y_centered, y_centered)))
    return float(np.dot(x_centered, y_centered) / denominator) if denominator else None


def spearman(x: list[float], y: list[float]) -> float | None:
    return pearson(average_ranks(x), average_ranks(y))


def load_results() -> list[dict[str, Any]]:
    results = []
    for path in sorted(SYSTEMS.glob("*/gnina_rescoring/*/result.json")):
        data = json.loads(path.read_text())
        if data.get("status") != "complete":
            raise RuntimeError(f"Incomplete GNINA result: {path}")
        if not data.get("qc", {}).get("coordinates_preserved"):
            raise RuntimeError(f"Coordinate QC failed: {path}")
        data["result_path"] = str(path.relative_to(ROOT))
        results.append(data)
    if len(results) != 45:
        raise RuntimeError(f"Expected 45 complete GNINA results, found {len(results)}")
    return results


def select_pose(poses: list[dict[str, Any]], metric: str) -> dict[str, Any]:
    return max(poses, key=lambda pose: (float(pose[metric]), -int(pose["pose_rank"])))


def flatten_results(results: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    pose_rows: list[dict[str, Any]] = []
    run_rows: list[dict[str, Any]] = []
    for result in results:
        poses = result["poses"]
        cnn_selected = select_pose(poses, "cnn_score")
        affinity_selected = select_pose(poses, "cnn_affinity")
        vina_top = min(poses, key=lambda pose: int(pose["pose_rank"]))
        vina_quality = [-float(pose["vina_affinity_kcal_mol"]) for pose in poses]
        cnn_scores = [float(pose["cnn_score"]) for pose in poses]
        cnn_affinities = [float(pose["cnn_affinity"]) for pose in poses]
        common = {
            "group": result["group"],
            "system": result["system"],
            "seed": int(result["seed"]),
            "source_run": result["source_run"],
            "result_path": result["result_path"],
        }
        for pose in poses:
            pose_rows.append({**common, **pose})
        run_rows.append(
            {
                **common,
                "n_poses": len(poses),
                "vina_top_affinity_kcal_mol": float(vina_top["vina_affinity_kcal_mol"]),
                "vina_top_rmsd_angstrom": vina_top["vina_rmsd_to_crystal_angstrom"],
                "vina_rank1_cnn_score": float(vina_top["cnn_score"]),
                "vina_rank1_cnn_affinity": float(vina_top["cnn_affinity"]),
                "cnn_selected_vina_rank": int(cnn_selected["pose_rank"]),
                "cnn_selected_cnn_score": float(cnn_selected["cnn_score"]),
                "cnn_selected_cnn_affinity": float(cnn_selected["cnn_affinity"]),
                "cnn_selected_vina_affinity_kcal_mol": float(cnn_selected["vina_affinity_kcal_mol"]),
                "cnn_selected_rmsd_angstrom": cnn_selected["vina_rmsd_to_crystal_angstrom"],
                "affinity_selected_vina_rank": int(affinity_selected["pose_rank"]),
                "affinity_selected_cnn_affinity": float(affinity_selected["cnn_affinity"]),
                "affinity_selected_rmsd_angstrom": affinity_selected["vina_rmsd_to_crystal_angstrom"],
                "cnn_vina_rank1_agree": int(cnn_selected["pose_rank"]) == 1,
                "spearman_vina_quality_vs_cnn_score": spearman(vina_quality, cnn_scores),
                "spearman_vina_quality_vs_cnn_affinity": spearman(vina_quality, cnn_affinities),
            }
        )
    return pose_rows, run_rows


def summarize_group(rows: list[dict[str, Any]]) -> dict[str, Any]:
    selected_rmsd = [float(row["cnn_selected_rmsd_angstrom"]) for row in rows if row["cnn_selected_rmsd_angstrom"] is not None]
    vina_rmsd = [float(row["vina_top_rmsd_angstrom"]) for row in rows if row["vina_top_rmsd_angstrom"] is not None]
    selected_scores = [float(row["cnn_selected_cnn_score"]) for row in rows]
    selected_affinity = [float(row["cnn_selected_cnn_affinity"]) for row in rows]
    rank1_scores = [float(row["vina_rank1_cnn_score"]) for row in rows]
    rank1_affinity = [float(row["vina_rank1_cnn_affinity"]) for row in rows]
    selected_ranks = [int(row["cnn_selected_vina_rank"]) for row in rows]
    score_correlations = [
        float(row["spearman_vina_quality_vs_cnn_score"])
        for row in rows
        if row["spearman_vina_quality_vs_cnn_score"] is not None
    ]
    affinity_correlations = [
        float(row["spearman_vina_quality_vs_cnn_affinity"])
        for row in rows
        if row["spearman_vina_quality_vs_cnn_affinity"] is not None
    ]
    return {
        "n_runs": len(rows),
        "n_poses": sum(int(row["n_poses"]) for row in rows),
        "mean_vina_top_rmsd_angstrom": round_or_none(mean(vina_rmsd)),
        "sd_vina_top_rmsd_angstrom": round_or_none(sample_sd(vina_rmsd)),
        "mean_cnn_selected_rmsd_angstrom": round_or_none(mean(selected_rmsd)),
        "sd_cnn_selected_rmsd_angstrom": round_or_none(sample_sd(selected_rmsd)),
        "cnn_selected_success_lt2A": sum(value < 2.0 for value in selected_rmsd) if selected_rmsd else None,
        "mean_vina_rank1_cnn_score": round_or_none(mean(rank1_scores)),
        "sd_vina_rank1_cnn_score": round_or_none(sample_sd(rank1_scores)),
        "mean_vina_rank1_cnn_affinity": round_or_none(mean(rank1_affinity)),
        "sd_vina_rank1_cnn_affinity": round_or_none(sample_sd(rank1_affinity)),
        "mean_cnn_selected_score": round_or_none(mean(selected_scores)),
        "sd_cnn_selected_score": round_or_none(sample_sd(selected_scores)),
        "mean_cnn_selected_affinity": round_or_none(mean(selected_affinity)),
        "sd_cnn_selected_affinity": round_or_none(sample_sd(selected_affinity)),
        "mean_vina_rank_selected_by_cnn": round_or_none(mean([float(value) for value in selected_ranks]), 3),
        "cnn_vina_rank1_agreement_count": sum(bool(row["cnn_vina_rank1_agree"]) for row in rows),
        "cnn_vina_rank1_agreement_fraction": round_or_none(
            mean([float(bool(row["cnn_vina_rank1_agree"])) for row in rows])
        ),
        "mean_within_run_spearman_vina_vs_cnn_score": round_or_none(mean(score_correlations)),
        "mean_within_run_spearman_vina_vs_cnn_affinity": round_or_none(mean(affinity_correlations)),
    }


def paired_delta(
    grouped_runs: dict[str, list[dict[str, Any]]], reference: str, alternate: str
) -> dict[str, Any]:
    reference_by_seed = {int(row["seed"]): row for row in grouped_runs[reference]}
    alternate_by_seed = {int(row["seed"]): row for row in grouped_runs[alternate]}
    seeds = sorted(set(reference_by_seed) & set(alternate_by_seed))
    rank1_affinity_deltas = [
        float(alternate_by_seed[seed]["vina_rank1_cnn_affinity"])
        - float(reference_by_seed[seed]["vina_rank1_cnn_affinity"])
        for seed in seeds
    ]
    rank1_score_deltas = [
        float(alternate_by_seed[seed]["vina_rank1_cnn_score"])
        - float(reference_by_seed[seed]["vina_rank1_cnn_score"])
        for seed in seeds
    ]
    selected_affinity_deltas = [
        float(alternate_by_seed[seed]["cnn_selected_cnn_affinity"])
        - float(reference_by_seed[seed]["cnn_selected_cnn_affinity"])
        for seed in seeds
    ]
    selected_score_deltas = [
        float(alternate_by_seed[seed]["cnn_selected_cnn_score"])
        - float(reference_by_seed[seed]["cnn_selected_cnn_score"])
        for seed in seeds
    ]
    return {
        "paired_seeds": seeds,
        "primary_rank1_mean_alternate_minus_canonical_cnn_affinity": round_or_none(mean(rank1_affinity_deltas)),
        "primary_rank1_sd_alternate_minus_canonical_cnn_affinity": round_or_none(sample_sd(rank1_affinity_deltas)),
        "primary_rank1_mean_alternate_minus_canonical_cnn_score": round_or_none(mean(rank1_score_deltas)),
        "primary_rank1_sd_alternate_minus_canonical_cnn_score": round_or_none(sample_sd(rank1_score_deltas)),
        "secondary_cnn_selected_mean_alternate_minus_canonical_cnn_affinity": round_or_none(mean(selected_affinity_deltas)),
        "secondary_cnn_selected_sd_alternate_minus_canonical_cnn_affinity": round_or_none(sample_sd(selected_affinity_deltas)),
        "secondary_cnn_selected_mean_alternate_minus_canonical_cnn_score": round_or_none(mean(selected_score_deltas)),
        "secondary_cnn_selected_sd_alternate_minus_canonical_cnn_score": round_or_none(sample_sd(selected_score_deltas)),
        "interpretation": "Primary deltas rescore the matched Vina rank-1 pose from each seed, avoiding bias from unequal pose counts. Negative CNNaffinity means a lower GNINA-predicted affinity for the alternate; CNNscore is pose confidence, not affinity.",
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titleweight": "bold",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.facecolor": "#FBFCFE",
            "figure.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, name: str) -> None:
    for extension in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=360, bbox_inches="tight")
    plt.close(fig)


def paired_validation_figure(grouped_runs: dict[str, list[dict[str, Any]]]) -> None:
    style()
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.7), sharey=True)
    for axis, group, letter in zip(axes, CANONICAL_GROUPS, "ABC"):
        rows = sorted(grouped_runs[group], key=lambda row: int(row["seed"]))
        vina_values = [float(row["vina_top_rmsd_angstrom"]) for row in rows]
        cnn_values = [float(row["cnn_selected_rmsd_angstrom"]) for row in rows]
        for vina_value, cnn_value in zip(vina_values, cnn_values):
            axis.plot([0, 1], [vina_value, cnn_value], color="#AAB4C3", lw=1.1, alpha=0.8, zorder=1)
        axis.scatter(np.zeros(len(rows)), vina_values, s=60, color="#667085", edgecolor="white", lw=0.8, zorder=3)
        axis.scatter(np.ones(len(rows)), cnn_values, s=67, color=COLORS[group], edgecolor="white", lw=0.8, zorder=3)
        axis.axhline(2.0, color="#B42318", linestyle="--", lw=1.1)
        axis.set_xticks([0, 1], ["Vina rank 1", "GNINA\nCNN-selected"])
        axis.set_title(f"{letter}  {group.replace(' canonical', '')}", loc="left")
        axis.grid(axis="y", color="#DCE2EA", lw=0.7)
        axis.text(
            0.98,
            0.96,
            f"GNINA <2 Å: {sum(value < 2 for value in cnn_values)}/{len(cnn_values)}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            color=COLORS[group],
            fontsize=9,
            fontweight="bold",
        )
    axes[0].set_ylabel("RMSD to crystallographic pose (Å)")
    fig.suptitle("Does CNN reranking preserve cognate pose recovery?", fontsize=16, fontweight="bold")
    fig.text(
        0.5,
        0.015,
        "Lines connect selections from the same Vina search. GNINA scored the existing poses without coordinate optimization.",
        ha="center",
        color="#475467",
        fontsize=9.5,
    )
    fig.tight_layout(rect=[0, 0.07, 1, 0.93], w_pad=2.0)
    save_figure(fig, "gnina_pose_selection_validation")


def comparative_figure(grouped_runs: dict[str, list[dict[str, Any]]]) -> None:
    style()
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2))

    def panel(axis: plt.Axes, groups: list[str], metric: str, ylabel: str, title: str) -> None:
        for index, group in enumerate(groups):
            values = [float(row[metric]) for row in grouped_runs[group] if row[metric] is not None]
            jitter = np.linspace(-0.12, 0.12, len(values))
            axis.scatter(
                np.full(len(values), index) + jitter,
                values,
                s=58,
                color=COLORS[group],
                edgecolor="white",
                linewidth=0.8,
                alpha=0.88,
                zorder=3,
            )
            value_mean = statistics.mean(values)
            value_sd = statistics.stdev(values) if len(values) > 1 else 0.0
            axis.errorbar(
                index,
                value_mean,
                yerr=value_sd,
                fmt="D",
                color="#101828",
                mfc="white",
                mec="#101828",
                markersize=7,
                capsize=5,
                linewidth=1.7,
                zorder=4,
            )
        labels = [
            group.replace(" canonical", "\ncanonical")
            .replace("CHRFAM7A at A face", "fusion at\nA face")
            .replace("CHRFAM7A at B face", "fusion at\nB face")
            for group in groups
        ]
        axis.set_xticks(range(len(groups)), labels)
        axis.set_ylabel(ylabel)
        axis.set_title(title, loc="left")
        axis.grid(axis="y", color="#DCE2EA", lw=0.7)

    bace = ["BACE1 canonical", "BACE1-476", "BACE1-457"]
    chrna = ["CHRNA7 canonical", "CHRFAM7A at A face", "CHRFAM7A at B face"]
    panel(axes[0, 0], bace, "vina_rank1_cnn_affinity", "CNNaffinity", "A  BACE1 Vina-rank-1 rescoring")
    panel(axes[0, 1], bace, "cnn_selected_rmsd_angstrom", "RMSD to 5HU1 pose (Å)", "B  BACE1 CNN-selected geometry")
    panel(axes[1, 0], chrna, "vina_rank1_cnn_affinity", "CNNaffinity", "C  α7 Vina-rank-1 rescoring")
    panel(axes[1, 1], chrna, "cnn_selected_rmsd_angstrom", "RMSD to 7EKP pose (Å)", "D  α7 CNN-selected geometry")
    axes[0, 1].axhline(2.0, color="#B42318", linestyle="--", lw=1.1)
    axes[1, 1].axhline(2.0, color="#B42318", linestyle="--", lw=1.1)
    fig.suptitle("GNINA rescoring of matched Vina pose ensembles", fontsize=16, fontweight="bold")
    fig.text(
        0.5,
        0.018,
        "Affinity panels rescore each Vina rank-1 pose; geometry panels show each highest-CNNscore pose. Diamonds show mean ± sample SD.",
        ha="center",
        color="#475467",
        fontsize=9.5,
    )
    fig.tight_layout(rect=[0, 0.055, 1, 0.95], h_pad=2.2, w_pad=2.0)
    save_figure(fig, "gnina_matched_comparison")


def rank_agreement_figure(grouped_runs: dict[str, list[dict[str, Any]]], summaries: dict[str, Any]) -> None:
    style()
    fig, (axis_rank, axis_corr) = plt.subplots(1, 2, figsize=(13.8, 5.5), gridspec_kw={"width_ratios": [1.2, 1]})
    for index, group in enumerate(GROUP_ORDER):
        ranks = [int(row["cnn_selected_vina_rank"]) for row in grouped_runs[group]]
        jitter = np.linspace(-0.11, 0.11, len(ranks))
        axis_rank.scatter(
            np.full(len(ranks), index) + jitter,
            ranks,
            s=54,
            color=COLORS[group],
            edgecolor="white",
            lw=0.7,
            alpha=0.9,
        )
    axis_rank.set_xticks(
        range(len(GROUP_ORDER)),
        [
            name.replace(" canonical", "\ncan.")
            .replace("CHRFAM7A at A face", "fusion\nA face")
            .replace("CHRFAM7A at B face", "fusion\nB face")
            for name in GROUP_ORDER
        ],
        rotation=20,
        ha="right",
    )
    axis_rank.set_ylabel("Original Vina rank selected by CNNscore")
    axis_rank.set_title("A  Pose-selection agreement", loc="left")
    axis_rank.axhline(1, color="#344054", lw=1.0, linestyle="--")
    axis_rank.grid(axis="y", color="#DCE2EA", lw=0.7)

    correlations = [summaries[group]["mean_within_run_spearman_vina_vs_cnn_score"] for group in GROUP_ORDER]
    y_positions = np.arange(len(GROUP_ORDER))
    axis_corr.barh(y_positions, correlations, color=[COLORS[group] for group in GROUP_ORDER], alpha=0.9)
    axis_corr.axvline(0, color="#344054", lw=0.9)
    axis_corr.set_yticks(y_positions, GROUP_ORDER)
    axis_corr.invert_yaxis()
    axis_corr.set_xlim(-1, 1)
    axis_corr.set_xlabel("Mean within-run Spearman ρ")
    axis_corr.set_title("B  Vina-quality vs CNNscore ranking", loc="left")
    axis_corr.grid(axis="x", color="#DCE2EA", lw=0.7)
    for y_position, value in zip(y_positions, correlations):
        axis_corr.text(
            value + (0.025 if value >= 0 else -0.025),
            y_position,
            f"{value:.2f}",
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=9,
        )
    fig.suptitle("Vina and GNINA provide related but non-identical pose rankings", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93], w_pad=2.8)
    save_figure(fig, "vina_gnina_rank_agreement")


def render_publication_figures() -> None:
    """Render the presentation suite with the repository's ggplot2 backend."""
    plotting_script = ROOT / "02_SURVEYOR" / "master_surveyor" / "plot_gnina_comparison.R"
    if not plotting_script.is_file():
        raise FileNotFoundError(plotting_script)
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "16",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    subprocess.run(
        ["taskset", "-c", "0-15", "Rscript", str(plotting_script)],
        check=True,
        cwd=ROOT,
        env=environment,
    )


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    results = load_results()
    pose_rows, run_rows = flatten_results(results)
    grouped_runs: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in run_rows:
        grouped_runs[str(row["group"])].append(row)
    for rows in grouped_runs.values():
        rows.sort(key=lambda row: int(row["seed"]))

    missing_groups = [group for group in GROUP_ORDER if group not in grouped_runs]
    if missing_groups:
        raise RuntimeError(f"Missing groups: {missing_groups}")
    expected_counts = {"FYN canonical": 9, **{group: 6 for group in GROUP_ORDER[1:]}}
    for group, expected in expected_counts.items():
        if len(grouped_runs[group]) != expected:
            raise RuntimeError(f"{group}: expected {expected} runs, found {len(grouped_runs[group])}")

    summaries = {group: summarize_group(grouped_runs[group]) for group in GROUP_ORDER}
    paired_deltas = {
        "BACE1-476": paired_delta(grouped_runs, "BACE1 canonical", "BACE1-476"),
        "BACE1-457": paired_delta(grouped_runs, "BACE1 canonical", "BACE1-457"),
        "CHRFAM7A at A face": paired_delta(grouped_runs, "CHRNA7 canonical", "CHRFAM7A at A face"),
        "CHRFAM7A at B face": paired_delta(grouped_runs, "CHRNA7 canonical", "CHRFAM7A at B face"),
    }
    coordinate_failures = sum(
        int(result["qc"]["heavy_atom_coordinate_mismatches_at_0.001_angstrom"])
        for result in results
    )
    summary = {
        "analysis_type": "Pose-preserving GNINA rescoring of retained AutoDock Vina pose ensembles",
        "run_count": len(run_rows),
        "pose_count": len(pose_rows),
        "coordinate_qc_failures": coordinate_failures,
        "groups": summaries,
        "paired_deltas": paired_deltas,
        "interpretation_constraints": [
            "CNNscore is a pose-confidence/ranking output, not binding affinity.",
            "CNNaffinity is an ML model prediction and must not be subtracted from Vina kcal/mol scores.",
            "Only within-system, preparation-matched canonical/alternate comparisons are interpreted.",
            "Rescoring the same poses provides scoring consensus, not independent conformational sampling.",
        ],
    }

    write_csv(ANALYSIS / "pose_scores.csv", pose_rows)
    write_csv(ANALYSIS / "run_summary.csv", run_rows)
    (ANALYSIS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    render_publication_figures()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

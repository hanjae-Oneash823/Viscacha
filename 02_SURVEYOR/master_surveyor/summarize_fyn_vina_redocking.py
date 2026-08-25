#!/usr/bin/env python3
"""Summarize independently seeded FYN--saracatinib Vina re-docking runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean, median


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for path in sorted(args.runs_dir.glob("*/result.json")):
        result = json.loads(path.read_text())
        for pose in result["results"]:
            rows.append({"run": path.parent.name, **result | pose})
    if not rows:
        raise SystemExit("No completed result.json files found")

    fieldnames = ["run", "seed", "exhaustiveness", "pose_rank", "vina_affinity_kcal_mol", "rmsd_to_crystal_heavy_atom_uncorrected_angstrom"]
    with (args.outdir / "all_poses.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    top = [r for r in rows if r["pose_rank"] == 1 and r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"] is not None]
    recovered = [r for r in top if r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"] < 2.0]
    summary = {
        "completed_runs": len(top),
        "top_pose_recovered_under_2A": len(recovered),
        "top_pose_recovery_rate_under_2A": round(len(recovered) / len(top), 4),
        "top_pose_rmsd_mean_A": round(mean(r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"] for r in top), 4),
        "top_pose_rmsd_median_A": round(median(r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"] for r in top), 4),
        "best_top_pose": min(top, key=lambda r: r["rmsd_to_crystal_heavy_atom_uncorrected_angstrom"]),
        "metric_note": "Fixed-frame heavy-atom RMSD derived from preserved PDBQT atom order; not symmetry corrected.",
        "scientific_note": "A successful re-docking validates this receptor/ligand/grid protocol for the canonical FYN kinase domain only. It does not establish binding affinity or response in cells.",
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

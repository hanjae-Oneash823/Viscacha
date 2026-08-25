#!/usr/bin/env python3
"""Run a reproducible AutoDock Vina re-docking calculation and report RMSD.

RMSD is calculated in the fixed experimental receptor frame using the atom
order preserved from the Meeko-prepared input PDBQT.  It is therefore an
initial, non-symmetry-corrected heavy-atom pose recovery metric.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from vina import Vina


def pdbqt_heavy_models(path: Path) -> list[np.ndarray]:
    models: list[list[list[float]]] = []
    current: list[list[float]] = []
    seen_model = False
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            if current:
                models.append(current)
                current = []
            seen_model = True
        elif line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:].strip() if len(line) >= 78 else ""
            if not atom_type.startswith("H"):
                current.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("ENDMDL") and current:
            models.append(current)
            current = []
    if current:
        models.append(current)
    if not seen_model and not models:
        raise ValueError(f"No coordinates parsed from {path}")
    return [np.asarray(m, dtype=float) for m in models]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--receptor", type=Path, required=True)
    p.add_argument("--ligand", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--center", nargs=3, type=float, required=True)
    p.add_argument("--box-size", nargs=3, type=float, required=True)
    p.add_argument("--seed", type=int, default=20260824)
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--exhaustiveness", type=int, default=16)
    p.add_argument("--n-poses", type=int, default=10)
    p.add_argument(
        "--skip-rmsd",
        action="store_true",
        help="Do not calculate fixed-frame RMSD (required for cross-docking a ligand generated outside the experimental frame).",
    )
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    vina = Vina(sf_name="vina", cpu=args.cpu, seed=args.seed, verbosity=1)
    vina.set_receptor(str(args.receptor))
    vina.set_ligand_from_file(str(args.ligand))
    vina.compute_vina_maps(center=args.center, box_size=args.box_size)
    vina.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)
    pose_path = args.outdir / "docked_poses.pdbqt"
    vina.write_poses(str(pose_path), n_poses=args.n_poses, overwrite=True)

    reference = None if args.skip_rmsd else pdbqt_heavy_models(args.ligand)[0]
    poses = pdbqt_heavy_models(pose_path)
    rmsds: list[float | None] = []
    for pose in poses:
        if reference is None or pose.shape != reference.shape:
            rmsds.append(None)
        else:
            rmsds.append(round(float(np.sqrt(np.mean(np.sum((pose - reference) ** 2, axis=1)))), 4))
    energies = vina.energies(n_poses=args.n_poses)
    rows = []
    for i, energy in enumerate(energies):
        rows.append({
            "pose_rank": i + 1,
            "vina_affinity_kcal_mol": round(float(energy[0]), 4),
            "rmsd_to_crystal_heavy_atom_uncorrected_angstrom": rmsds[i] if i < len(rmsds) else None,
        })
    output = {
        "receptor": str(args.receptor), "ligand": str(args.ligand),
        "center": args.center, "box_size": args.box_size, "seed": args.seed,
        "cpu": args.cpu, "exhaustiveness": args.exhaustiveness,
        "n_poses_requested": args.n_poses, "results": rows,
        "rmsd_note": (
            "Not calculated: cross-docked ligand input did not begin in the experimental receptor frame."
            if args.skip_rmsd
            else "Fixed-frame heavy-atom RMSD; atom-order based and not symmetry corrected."
        ),
    }
    (args.outdir / "result.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

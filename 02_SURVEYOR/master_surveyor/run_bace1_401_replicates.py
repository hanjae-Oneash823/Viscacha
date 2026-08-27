#!/usr/bin/env python3
"""Run matched BACE1-401/verubecestat Vina replicates with a 16-CPU cap."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OLD_SYSTEM = ROOT / "outputs" / "docking_campaign" / "systems" / "BACE1_verubecestat"
SYSTEM = ROOT / "outputs" / "docking_campaign" / "systems" / "BACE1_verubecestat_401"
RUNNER = ROOT / "02_SURVEYOR" / "master_surveyor" / "run_vina_redock.py"
PYTHON = Path("/home/welcome3/anaconda3/envs/pocket_dock/bin/python")
SEEDS = (1103, 2207, 3301, 4409, 5519)


def main() -> None:
    allowed = sorted(os.sched_getaffinity(0))[:16]
    os.sched_setaffinity(0, allowed)
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "16",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    receptor = SYSTEM / "prepared" / "bace1_401_aligned_to_5HU1_obabel.pdbqt"
    ligand = OLD_SYSTEM / "prepared" / "verubecestat_66F_A501.pdbqt"
    if not receptor.exists() or not ligand.exists():
        raise FileNotFoundError(f"Missing receptor or ligand: {receptor}, {ligand}")

    for seed in SEEDS:
        label = f"alternate_401_seed{seed}_ex32"
        output = SYSTEM / "runs" / label
        command = [
            str(PYTHON),
            str(RUNNER),
            "--receptor",
            str(receptor),
            "--ligand",
            str(ligand),
            "--outdir",
            str(output),
            "--center",
            "24.887",
            "10.422",
            "21.858",
            "--box-size",
            "23.129",
            "16.370",
            "18.460",
            "--seed",
            str(seed),
            "--cpu",
            "16",
            "--exhaustiveness",
            "32",
            "--n-poses",
            "20",
        ]
        print(f"START {label}; CPU affinity={allowed}", flush=True)
        completed = subprocess.run(command, env=environment, text=True, capture_output=True)
        log_dir = SYSTEM / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        (log_dir / f"{label}.log").write_text(completed.stdout + completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"Failed {label}; see {log_dir / f'{label}.log'}")
        print(f"DONE  {label}", flush=True)


if __name__ == "__main__":
    main()

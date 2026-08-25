#!/usr/bin/env python3
"""Run expanded-campaign Vina replicates serially with a hard 16-CPU cap."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
SYSTEMS = CAMPAIGN / "systems"
PYTHON = Path("/home/welcome3/anaconda3/envs/pocket_dock/bin/python")
RUNNER = ROOT / "02_SURVEYOR" / "master_surveyor" / "run_vina_redock.py"
SEEDS = [1103, 2207, 3301, 4409, 5519]

CONFIGS = [
    {
        "candidate": "BACE1_verubecestat", "label": "canonical_obabel",
        "receptor": "prepared/bace1_5HU1_chain_A.pdbqt", "ligand": "prepared/verubecestat_66F_A501.pdbqt",
        "center": [24.887, 10.422, 21.858], "box": [23.129, 16.370, 18.460], "rmsd": True,
    },
    {
        "candidate": "BACE1_verubecestat", "label": "alternate_476",
        "receptor": "prepared/bace1_476_alphafold_aligned_to_5HU1_obabel.pdbqt", "ligand": "prepared/verubecestat_66F_A501.pdbqt",
        "center": [24.887, 10.422, 21.858], "box": [23.129, 16.370, 18.460], "rmsd": True,
    },
    {
        "candidate": "BACE1_verubecestat", "label": "alternate_457",
        "receptor": "prepared/bace1_457_alphafold_aligned_to_5HU1_obabel.pdbqt", "ligand": "prepared/verubecestat_66F_A501.pdbqt",
        "center": [24.887, 10.422, 21.858], "box": [23.129, 16.370, 18.460], "rmsd": True,
    },
    {
        "candidate": "CHRNA7_encenicline", "label": "canonical_obabel",
        "receptor": "prepared/chrna7_7EKP_pentamer.pdbqt", "ligand": "prepared/encenicline_I33_A601.pdbqt",
        "center": [142.818, 138.419, 90.348], "box": [22.027, 19.752, 16.264], "rmsd": True,
    },
    {
        "candidate": "CHRNA7_encenicline", "label": "hybrid_A_face",
        "receptor": "prepared/chrna7_chrFam7a_fusion_at_A_face_obabel.pdbqt", "ligand": "prepared/encenicline_I33_A601.pdbqt",
        "center": [142.818, 138.419, 90.348], "box": [22.027, 19.752, 16.264], "rmsd": True,
    },
    {
        "candidate": "CHRNA7_encenicline", "label": "hybrid_B_face",
        "receptor": "prepared/chrna7_chrFam7a_fusion_at_B_face_obabel.pdbqt", "ligand": "prepared/encenicline_I33_A601.pdbqt",
        "center": [142.818, 138.419, 90.348], "box": [22.027, 19.752, 16.264], "rmsd": True,
    },
    {
        "candidate": "GABRA2_AZD7325", "label": "canonical_crossdock",
        "receptor": "prepared/gabra_9CSB_native_pentamer_meeko.pdbqt", "ligand": "prepared/AZD7325.pdbqt",
        "center": [143.009, 97.405, 140.440], "box": [19.879, 15.894, 19.575], "rmsd": False,
    },
    {
        "candidate": "CACNA1D_isradipine", "label": "shared_pocket_crossdock",
        "receptor": "prepared/cacna1d_8E59_chain_A.pdbqt", "ligand": "prepared/isradipine.pdbqt",
        "center": [151.334, 167.442, 149.793], "box": [18.851, 27.189, 19.189], "rmsd": False,
    },
    {
        "candidate": "PDE9A_BI409306", "label": "canonical_crossdock",
        "receptor": "prepared/pde9a_4GH6_chain_A_meeko.pdbqt", "ligand": "prepared/BI409306.pdbqt",
        "center": [77.260, 51.202, 39.870], "box": [19.790, 18.581, 23.684], "rmsd": False,
    },
]


def main() -> None:
    allowed = sorted(os.sched_getaffinity(0))[:16]
    os.sched_setaffinity(0, allowed)
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    })
    print(f"CPU affinity: {allowed}; jobs are serial; Vina cpu=16", flush=True)
    for config in CONFIGS:
        base = SYSTEMS / config["candidate"]
        for seed in SEEDS:
            outdir = base / "runs" / f"{config['label']}_seed{seed}_ex32"
            command = [
                str(PYTHON), str(RUNNER),
                "--receptor", str(base / config["receptor"]),
                "--ligand", str(base / config["ligand"]),
                "--outdir", str(outdir),
                "--center", *(str(x) for x in config["center"]),
                "--box-size", *(str(x) for x in config["box"]),
                "--seed", str(seed), "--cpu", "16", "--exhaustiveness", "32", "--n-poses", "20",
            ]
            if not config["rmsd"]:
                command.append("--skip-rmsd")
            print(f"START {config['candidate']} {config['label']} seed={seed}", flush=True)
            completed = subprocess.run(command, env=env, text=True, capture_output=True)
            log_dir = base / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log = log_dir / f"{config['label']}_seed{seed}_ex32.log"
            log.write_text(completed.stdout + completed.stderr)
            if completed.returncode != 0:
                raise RuntimeError(f"Failed: {' '.join(command)}; see {log}")
            print(f"DONE  {config['candidate']} {config['label']} seed={seed}", flush=True)


if __name__ == "__main__":
    main()

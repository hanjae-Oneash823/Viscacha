#!/usr/bin/env python3
"""Rescore retained Vina poses with GNINA while preserving pose coordinates.

The workflow deliberately uses GNINA ``--score_only`` rather than a second
docking search.  Every Vina pose is converted from the multi-model PDBQT to a
multi-record SDF, scored with the GNINA default CNN ensemble, and checked to
ensure that the heavy-atom coordinate set is unchanged.  Jobs run serially and
the Docker container is restricted to at most 16 logical CPUs.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit import RDLogger


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
SYSTEMS = CAMPAIGN / "systems"
ANALYSIS = CAMPAIGN / "analysis" / "gnina"

DEFAULT_BINARY = Path("/home/welcome3/.local/bin/gnina-1.3.3.cuda12.8.static")
DEFAULT_IMAGE = "nvidia/cuda:12.8.1-cudnn-runtime-ubuntu22.04"
EXPECTED_BINARY_SHA256 = "3340c1f49cd3c7c84d8699182a1c6af13c7fa2a22448d1204640446106f72172"
CPU_LIMIT = 16
CPUSET = "0-15"


@dataclass(frozen=True)
class GroupConfig:
    group: str
    system: str
    exact_runs: tuple[str, ...] = ()
    run_patterns: tuple[str, ...] = ()


GROUP_CONFIGS = (
    GroupConfig(
        "FYN canonical",
        "FYN_saracatinib",
        run_patterns=("vina_initial_seed*_ex*", "vina_seed*_ex*"),
    ),
    GroupConfig(
        "BACE1 canonical",
        "BACE1_verubecestat",
        exact_runs=("canonical_obabel_redock_seed20260825_ex32",),
        run_patterns=("canonical_obabel_seed*_ex32",),
    ),
    GroupConfig(
        "BACE1-476",
        "BACE1_verubecestat",
        run_patterns=("alternate_476_seed*_ex32",),
    ),
    GroupConfig(
        "BACE1-457",
        "BACE1_verubecestat",
        run_patterns=("alternate_457_seed*_ex32",),
    ),
    GroupConfig(
        "BACE1-401 observed AD isoform",
        "BACE1_verubecestat_401",
        run_patterns=("alternate_401_seed*_ex32",),
    ),
    GroupConfig(
        "CHRNA7 canonical",
        "CHRNA7_encenicline",
        exact_runs=("canonical_obabel_redock_seed20260825_ex32",),
        run_patterns=("canonical_obabel_seed*_ex32",),
    ),
    GroupConfig(
        "CHRFAM7A at A face",
        "CHRNA7_encenicline",
        run_patterns=("hybrid_A_face_seed*_ex32",),
    ),
    GroupConfig(
        "CHRFAM7A at B face",
        "CHRNA7_encenicline",
        run_patterns=("hybrid_B_face_seed*_ex32",),
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=check)


def relative_to_root(path_string: str) -> Path:
    path = Path(path_string)
    if not path.is_absolute():
        path = ROOT / path
    try:
        return path.resolve().relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError(f"Input is outside repository root: {path}") from error


def discover_jobs(pilot: bool, group_filter: str | None = None) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for config in GROUP_CONFIGS:
        runs_dir = SYSTEMS / config.system / "runs"
        candidates = [runs_dir / name for name in config.exact_runs]
        for pattern in config.run_patterns:
            candidates.extend(sorted(runs_dir.glob(pattern)))
        for source_run in candidates:
            result_path = source_run / "result.json"
            pose_path = source_run / "docked_poses.pdbqt"
            if source_run in seen or not result_path.exists() or not pose_path.exists():
                continue
            seen.add(source_run)
            data = json.loads(result_path.read_text())
            jobs.append(
                {
                    "group": config.group,
                    "system": config.system,
                    "source_run": source_run,
                    "source_result": result_path,
                    "source_poses": pose_path,
                    "receptor": ROOT / relative_to_root(str(data["receptor"])),
                    "seed": int(data["seed"]),
                    "source_data": data,
                }
            )
    jobs.sort(key=lambda item: (str(item["system"]), str(item["group"]), int(item["seed"])))
    if pilot:
        jobs = [
            job
            for job in jobs
            if job["group"] == "BACE1 canonical" and job["seed"] == 1103
        ]
    if group_filter:
        jobs = [job for job in jobs if job["group"] == group_filter]
    return jobs


def write_pose_sdf(job: dict[str, Any], output_path: Path) -> int:
    pdbqt = PDBQTMolecule.from_file(str(job["source_poses"]), skip_typing=True)
    molecules = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    if len(molecules) != 1 or molecules[0] is None:
        raise RuntimeError(f"Expected one ligand molecule in {job['source_poses']}; got {len(molecules)}")
    molecule = molecules[0]
    results = job["source_data"]["results"]
    if molecule.GetNumConformers() != len(results):
        raise RuntimeError(
            f"Pose count mismatch for {job['source_run']}: "
            f"{molecule.GetNumConformers()} conformers versus {len(results)} JSON records"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_path))
    for index, result in enumerate(results):
        pose_rank = int(result["pose_rank"])
        molecule.SetProp(
            "_Name",
            f"{job['system']}__{job['source_run'].name}__pose{pose_rank:02d}",
        )
        molecule.SetProp("source_group", str(job["group"]))
        molecule.SetProp("source_system", str(job["system"]))
        molecule.SetProp("source_run", job["source_run"].name)
        molecule.SetIntProp("source_seed", int(job["seed"]))
        molecule.SetIntProp("vina_pose_rank", pose_rank)
        molecule.SetDoubleProp("vina_affinity_kcal_mol", float(result["vina_affinity_kcal_mol"]))
        rmsd = result.get("rmsd_to_crystal_heavy_atom_uncorrected_angstrom")
        if rmsd is not None:
            molecule.SetDoubleProp("vina_rmsd_to_crystal_angstrom", float(rmsd))
        writer.write(molecule, confId=index)
    writer.close()
    return len(results)


def heavy_atom_coordinate_signature(molecule: Chem.Mol) -> list[tuple[Any, ...]]:
    positions = molecule.GetConformer().GetPositions()
    return sorted(
        (
            atom.GetSymbol(),
            round(float(positions[atom.GetIdx()][0]), 3),
            round(float(positions[atom.GetIdx()][1]), 3),
            round(float(positions[atom.GetIdx()][2]), 3),
        )
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1
    )


def read_sdf(path: Path) -> list[Chem.Mol]:
    return [molecule for molecule in Chem.SDMolSupplier(str(path), removeHs=False) if molecule is not None]


def parse_scored_poses(
    job: dict[str, Any], input_path: Path, scored_path: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    input_molecules = read_sdf(input_path)
    scored_molecules = read_sdf(scored_path)
    if len(input_molecules) != len(scored_molecules):
        raise RuntimeError(
            f"GNINA output count mismatch for {job['source_run']}: "
            f"{len(input_molecules)} input versus {len(scored_molecules)} scored"
        )

    rows: list[dict[str, Any]] = []
    coordinate_mismatches = 0
    required = ("CNNscore", "CNNaffinity", "CNNaffinity_variance", "CNN_VS", "minimizedAffinity")
    for source_molecule, scored_molecule in zip(input_molecules, scored_molecules):
        missing = [name for name in required if not scored_molecule.HasProp(name)]
        if missing:
            raise RuntimeError(f"GNINA output lacks {missing} in {scored_path}")
        if heavy_atom_coordinate_signature(source_molecule) != heavy_atom_coordinate_signature(scored_molecule):
            coordinate_mismatches += 1
        rmsd = (
            float(scored_molecule.GetProp("vina_rmsd_to_crystal_angstrom"))
            if scored_molecule.HasProp("vina_rmsd_to_crystal_angstrom")
            else None
        )
        rows.append(
            {
                "pose_rank": int(scored_molecule.GetProp("vina_pose_rank")),
                "vina_affinity_kcal_mol": float(scored_molecule.GetProp("vina_affinity_kcal_mol")),
                "vina_rmsd_to_crystal_angstrom": rmsd,
                "gnina_score_only_vina_kcal_mol": float(scored_molecule.GetProp("minimizedAffinity")),
                "cnn_score": float(scored_molecule.GetProp("CNNscore")),
                "cnn_affinity": float(scored_molecule.GetProp("CNNaffinity")),
                "cnn_affinity_variance": float(scored_molecule.GetProp("CNNaffinity_variance")),
                "cnn_vs": float(scored_molecule.GetProp("CNN_VS")),
            }
        )

    rows.sort(key=lambda row: int(row["pose_rank"]))
    qc = {
        "input_pose_count": len(input_molecules),
        "scored_pose_count": len(scored_molecules),
        "heavy_atom_coordinate_mismatches_at_0.001_angstrom": coordinate_mismatches,
        "coordinates_preserved": coordinate_mismatches == 0,
    }
    if coordinate_mismatches:
        raise RuntimeError(f"GNINA changed coordinates for {coordinate_mismatches} poses in {job['source_run']}")
    return rows, qc


def container_path(host_path: Path) -> str:
    return "/work/" + str(host_path.resolve().relative_to(ROOT.resolve()))


def build_gnina_command(
    *, binary: Path, image: str, device: int, receptor: Path, output_dir: Path
) -> list[str]:
    uid = os.getuid()
    gid = os.getgid()
    return [
        "docker",
        "run",
        "--rm",
        "--gpus",
        f"device={device}",
        "--cpus",
        str(CPU_LIMIT),
        "--cpuset-cpus",
        CPUSET,
        "--user",
        f"{uid}:{gid}",
        "-e",
        f"OMP_NUM_THREADS={CPU_LIMIT}",
        "-e",
        "OPENBLAS_NUM_THREADS=1",
        "-e",
        "MKL_NUM_THREADS=1",
        "-v",
        f"{binary.resolve()}:/opt/gnina:ro",
        "-v",
        f"{ROOT.resolve()}:/work:ro",
        "-v",
        f"{output_dir.resolve()}:/out:rw",
        image,
        "/opt/gnina",
        "--receptor",
        container_path(receptor),
        "--ligand",
        "/out/input_poses.sdf",
        "--score_only",
        "--cnn_scoring",
        "rescore",
        "--cpu",
        str(CPU_LIMIT),
        "--device",
        "0",
        "--out",
        "/out/scored_poses.sdf",
    ]


def collect_runtime_metadata(binary: Path, image: str, device: int) -> dict[str, Any]:
    image_info = run_command(
        ["docker", "image", "inspect", image, "--format", "{{index .RepoDigests 0}}|{{.Id}}|{{.Size}}"]
    ).stdout.strip()
    image_digest, image_id, image_size = image_info.split("|")
    version_command = [
        "docker", "run", "--rm", "--gpus", f"device={device}", "--cpus", str(CPU_LIMIT),
        "--cpuset-cpus", CPUSET, "-v", f"{binary.resolve()}:/opt/gnina:ro", image,
        "/opt/gnina", "--version",
    ]
    version_result = run_command(version_command)
    version_text = version_result.stdout + version_result.stderr
    version_line = next((line for line in version_text.splitlines() if line.startswith("gnina v")), version_text.strip())
    gpu_result = run_command(
        [
            "nvidia-smi", f"--id={device}",
            "--query-gpu=index,name,driver_version,memory.total",
            "--format=csv,noheader,nounits",
        ]
    )
    git_commit = run_command(["git", "rev-parse", "HEAD"]).stdout.strip()
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": "Pose-preserving GNINA score_only with the built-in default CNN ensemble",
        "gnina_version": version_line,
        "gnina_binary": str(binary),
        "gnina_binary_sha256": sha256(binary),
        "expected_binary_sha256": EXPECTED_BINARY_SHA256,
        "binary_checksum_matches_release": sha256(binary) == EXPECTED_BINARY_SHA256,
        "container_image": image,
        "container_repo_digest": image_digest,
        "container_image_id": image_id,
        "container_image_size_bytes": int(image_size),
        "gpu": gpu_result.stdout.strip(),
        "requested_gpu_device": device,
        "container_visible_device": 0,
        "cpu_limit": CPU_LIMIT,
        "cpuset": CPUSET,
        "jobs_serial": True,
        "cnn_scoring": "rescore",
        "cnn_model": "built-in default ensemble (implicit; no --cnn override)",
        "coordinate_optimization": False,
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "meeko_version": importlib.metadata.version("meeko"),
        "rdkit_version": importlib.metadata.version("rdkit"),
        "repository_commit_at_launch": git_commit,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true", help="Run only BACE1 canonical seed 1103")
    parser.add_argument("--group", help="Run only jobs whose group label exactly matches this value")
    parser.add_argument("--overwrite", action="store_true", help="Recompute jobs with existing successful result.json")
    parser.add_argument("--device", type=int, default=0, help="Host GPU index (default: 0)")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--image", default=DEFAULT_IMAGE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.binary.is_file():
        raise FileNotFoundError(args.binary)
    binary_digest = sha256(args.binary)
    if binary_digest != EXPECTED_BINARY_SHA256:
        raise RuntimeError(f"GNINA binary checksum mismatch: {binary_digest}")

    allowed = sorted(os.sched_getaffinity(0))[:CPU_LIMIT]
    os.sched_setaffinity(0, allowed)
    os.environ.update(
        {
            "OMP_NUM_THREADS": str(CPU_LIMIT),
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    RDLogger.DisableLog("rdApp.warning")

    jobs = discover_jobs(args.pilot, args.group)
    if not jobs:
        raise RuntimeError("No GNINA jobs matched the requested filters")

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    runtime_metadata = collect_runtime_metadata(args.binary, args.image, args.device)
    runtime_metadata["job_count"] = len(jobs)
    runtime_metadata["pilot"] = args.pilot
    metadata_name = "runtime_metadata.json"
    if args.group:
        safe_group = re.sub(r"[^A-Za-z0-9]+", "_", args.group).strip("_").lower()
        metadata_name = f"runtime_metadata_{safe_group}.json"
    (ANALYSIS / metadata_name).write_text(json.dumps(runtime_metadata, indent=2) + "\n")

    print(
        f"GNINA jobs={len(jobs)}; host affinity={allowed}; Docker cpuset={CPUSET}; "
        f"cpu={CPU_LIMIT}; GPU={args.device}; jobs are serial",
        flush=True,
    )
    completed_count = 0
    skipped_count = 0
    for index, job in enumerate(jobs, start=1):
        output_dir = SYSTEMS / job["system"] / "gnina_rescoring" / job["source_run"].name
        result_path = output_dir / "result.json"
        if result_path.exists() and not args.overwrite:
            existing = json.loads(result_path.read_text())
            if existing.get("status") == "complete" and existing.get("qc", {}).get("coordinates_preserved"):
                print(f"SKIP {index:02d}/{len(jobs)} {job['group']} seed={job['seed']}", flush=True)
                skipped_count += 1
                continue

        output_dir.mkdir(parents=True, exist_ok=True)
        input_path = output_dir / "input_poses.sdf"
        scored_path = output_dir / "scored_poses.sdf"
        pose_count = write_pose_sdf(job, input_path)
        command = build_gnina_command(
            binary=args.binary,
            image=args.image,
            device=args.device,
            receptor=job["receptor"],
            output_dir=output_dir,
        )
        print(
            f"START {index:02d}/{len(jobs)} {job['group']} seed={job['seed']} poses={pose_count}",
            flush=True,
        )
        started = datetime.now(timezone.utc)
        completed = run_command(command, check=False)
        ended = datetime.now(timezone.utc)
        log_path = output_dir / "gnina.log"
        log_path.write_text(completed.stdout + completed.stderr)
        if completed.returncode != 0:
            failure = {
                "status": "failed",
                "returncode": completed.returncode,
                "group": job["group"],
                "system": job["system"],
                "seed": job["seed"],
                "source_result": str(job["source_result"].relative_to(ROOT)),
                "log": str(log_path.relative_to(ROOT)),
            }
            result_path.write_text(json.dumps(failure, indent=2) + "\n")
            raise RuntimeError(f"GNINA failed for {job['source_run']}; see {log_path}")

        rows, qc = parse_scored_poses(job, input_path, scored_path)
        result = {
            "status": "complete",
            "group": job["group"],
            "system": job["system"],
            "seed": job["seed"],
            "source_run": job["source_run"].name,
            "source_result": str(job["source_result"].relative_to(ROOT)),
            "source_poses": str(job["source_poses"].relative_to(ROOT)),
            "receptor": str(job["receptor"].relative_to(ROOT)),
            "input_poses_sdf": str(input_path.relative_to(ROOT)),
            "scored_poses_sdf": str(scored_path.relative_to(ROOT)),
            "log": str(log_path.relative_to(ROOT)),
            "started_utc": started.isoformat(),
            "ended_utc": ended.isoformat(),
            "elapsed_seconds": round((ended - started).total_seconds(), 3),
            "protocol": {
                "gnina_mode": "score_only",
                "cnn_scoring": "rescore",
                "cnn_model": "built-in default ensemble (implicit; no --cnn override)",
                "coordinate_optimization": False,
                "cpu": CPU_LIMIT,
                "cpuset": CPUSET,
                "host_gpu_device": args.device,
            },
            "qc": qc,
            "poses": rows,
        }
        result_path.write_text(json.dumps(result, indent=2) + "\n")
        completed_count += 1
        print(
            f"DONE  {index:02d}/{len(jobs)} {job['group']} seed={job['seed']} "
            f"coordinate_qc=PASS elapsed={result['elapsed_seconds']:.1f}s",
            flush=True,
        )

    print(f"Complete: new={completed_count}, skipped={skipped_count}, total={len(jobs)}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise

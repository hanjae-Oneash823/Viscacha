"""M2 — protein structures: canonical via AFDB/ColabFold, alt/isoform always
locally predicted with a materially different protocol.

Canonical sequence: AlphaFold DB REST lookup by uniprot_acc first (free,
instant, no GPU -- the canonical protein IS the plain UniProt/MANE sequence
by construction, so AFDB usually already has it). colabfold_batch with
default settings (templates on) only on an AFDB miss -- a template hit here
is expected and desirable, this is the well-studied fold.

Alt/isoform sequence: ALWAYS predicted locally, never via AFDB (it's novel
by construction). Uses a protocol that differs from a normal single-protein
fold specifically because this is an isoform-structure question:
  - templates OFF -- otherwise ColabFold tends to pull the canonical
    protein's own crystal structure as a template and bias the model toward
    the fold this analysis exists to question.
  - --num-recycle well above the AFDB default (elevated recycles help the
    model settle on a self-consistent structure rather than an early,
    template-biased guess), following the published no-template +
    elevated-recycle protocol used for AD-variant (Presenilin-1) structure
    prediction.
  - --num-seeds > 1 for a small ensemble -- the isoform-altered region is
    precisely where model-to-model variability signals genuine structural
    ambiguity, worth seeing, not averaging away.
  - a single ESMFold (MSA-free) prediction as a cross-check specifically for
    the splice-altered span: MMseqs2 pulls homologs of the CANONICAL
    sequence, so a truncated/frameshifted/novel alt region can have
    thin-to-zero aligned homology, starving AF2's co-evolutionary signal
    exactly where it matters most.

Everything is cached by sequence hash under STRUCTURE_CACHE_DIR, so identical
sequences (e.g. the same alt transcript hit in multiple cell types) are
computed once and reused.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor.config import (
    AFDB_API, AFDB_TIMEOUT, COLABFOLD_ALT_ARGS, COLABFOLD_BIN,
    COLABFOLD_CANONICAL_ARGS, COLABFOLD_TIMEOUT_S, CUDA_VISIBLE_DEVICES,
    ESMFOLD_CHUNK_SIZE, ESMFOLD_ENV_PYTHON, ESMFOLD_SERVER_FOLD_TIMEOUT_S,
    ESMFOLD_SERVER_HEALTHCHECK_TIMEOUT_S, ESMFOLD_SERVER_URL,
    FASTA_CACHE_DIR, STRUCTURE_CACHE_DIR,
)
from master_surveyor.utils.http import HTTPError, get_json

_SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
_ESMFOLD_SCRIPT = _SCRIPTS_DIR / "esmfold_predict.py"


def _log(msg: str) -> None:
    print(f"  [m2] {msg}", file=sys.stderr, flush=True)


def seq_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()[:16]


@dataclass
class CanonicalStructure:
    path: Path
    source: str   # "afdb" | "colabfold"


@dataclass
class AltStructureResult:
    seed_models: list[Path]         # top-1 model per seed, the ensemble
    scores_json: list[Path]         # matching per-seed scores (plddt + pae)
    esmfold_model: Path | None      # None if ESMFold prediction failed
    esmfold_error: str | None = None


# ---------------------------------------------------------------------------
# Canonical
# ---------------------------------------------------------------------------

def _fetch_afdb_pdb(uniprot_acc: str) -> str | None:
    try:
        entries = get_json(AFDB_API.format(acc=uniprot_acc), timeout=AFDB_TIMEOUT, retries=1)
    except HTTPError as e:
        if e.status == 404:
            return None
        _log(f"{uniprot_acc}: AFDB lookup failed ({e})")
        return None
    if not entries:
        return None
    import requests
    pdb_resp = requests.get(entries[0]["pdbUrl"], timeout=AFDB_TIMEOUT * 2)
    pdb_resp.raise_for_status()
    return pdb_resp.text


def get_canonical_structure(uniprot_acc: str, sequence: str) -> CanonicalStructure:
    h = seq_hash(sequence)
    out_dir = STRUCTURE_CACHE_DIR / h
    afdb_path = out_dir / "canonical_afdb.pdb"
    if afdb_path.exists():
        return CanonicalStructure(afdb_path, "afdb")

    colabfold_dir = out_dir / "canonical_colabfold"
    existing = sorted(colabfold_dir.glob("*_unrelaxed_rank_001_*.pdb")) if colabfold_dir.exists() else []
    if existing:
        return CanonicalStructure(existing[0], "colabfold")

    pdb_text = _fetch_afdb_pdb(uniprot_acc) if uniprot_acc else None
    if pdb_text:
        out_dir.mkdir(parents=True, exist_ok=True)
        afdb_path.write_text(pdb_text)
        _log(f"{uniprot_acc}: canonical structure from AFDB")
        return CanonicalStructure(afdb_path, "afdb")

    _log(f"{uniprot_acc or '(no accession)'}: AFDB miss, folding canonical via ColabFold")
    _run_colabfold(sequence, colabfold_dir, COLABFOLD_CANONICAL_ARGS, jobname="canonical")
    models = sorted(colabfold_dir.glob("*_unrelaxed_rank_001_*.pdb"))
    if not models:
        raise RuntimeError(f"ColabFold produced no rank_001 model for canonical structure in {colabfold_dir}")
    return CanonicalStructure(models[0], "colabfold")


# ---------------------------------------------------------------------------
# Alt / isoform
# ---------------------------------------------------------------------------

def get_alt_structure(sequence: str) -> AltStructureResult:
    """COLABFOLD_ALT_ARGS pins --num-models 1, so --num-seeds N produces
    exactly N predictions total, globally ranked rank_001..rank_00N -- one
    per seed, unambiguously (verified empirically: colabfold_batch ranks
    ALL model x seed combinations together in one global list, so without
    --num-models 1 there would be no clean "one file per seed" glob).
    """
    h = seq_hash(sequence)
    out_dir = STRUCTURE_CACHE_DIR / h / "alt_colabfold"

    seed_models = sorted(out_dir.glob("*_unrelaxed_rank_*.pdb")) if out_dir.exists() else []
    if not seed_models:
        _log(f"{h}: folding alt/isoform sequence (no templates, elevated recycles, multi-seed ensemble)")
        _run_colabfold(sequence, out_dir, COLABFOLD_ALT_ARGS, jobname="alt")
        seed_models = sorted(out_dir.glob("*_unrelaxed_rank_*.pdb"))
        if not seed_models:
            raise RuntimeError(f"ColabFold produced no rank_* models in {out_dir}")
    scores_json = [_scores_path_for(model) for model in seed_models]

    esmfold_path = STRUCTURE_CACHE_DIR / h / "esmfold.pdb"
    esmfold_error = None
    if not esmfold_path.exists():
        esmfold_error = _run_esmfold(sequence, esmfold_path)

    return AltStructureResult(
        seed_models=seed_models,
        scores_json=[p for p in scores_json if p is not None],
        esmfold_model=esmfold_path if esmfold_path.exists() else None,
        esmfold_error=esmfold_error,
    )


def _scores_path_for(model_pdb: Path) -> Path | None:
    """colabfold_batch's scores json shares the model's stem, swapping the
    '_unrelaxed_' segment for '_scores_' and the .pdb extension for .json.
    """
    candidate = model_pdb.parent / model_pdb.name.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json")
    return candidate if candidate.exists() else None


# ---------------------------------------------------------------------------
# Subprocess runners
# ---------------------------------------------------------------------------

def _run_colabfold(sequence: str, out_dir: Path, extra_args: list[str], jobname: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    FASTA_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    fasta_path = FASTA_CACHE_DIR / f"{jobname}_{seq_hash(sequence)}.fasta"
    fasta_path.write_text(f">{jobname}_{seq_hash(sequence)}\n{sequence}\n")

    env = {**os.environ, "CUDA_VISIBLE_DEVICES": CUDA_VISIBLE_DEVICES}
    cmd = [COLABFOLD_BIN, *extra_args, str(fasta_path), str(out_dir)]
    _log(f"running: {' '.join(cmd)}")
    result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=COLABFOLD_TIMEOUT_S)
    if result.returncode != 0:
        raise RuntimeError(
            f"colabfold_batch failed (exit {result.returncode}) for {jobname} {seq_hash(sequence)}:\n"
            f"{result.stderr[-4000:]}"
        )


def _esmfold_server_available() -> bool:
    import requests
    try:
        resp = requests.get(f"{ESMFOLD_SERVER_URL}/health", timeout=ESMFOLD_SERVER_HEALTHCHECK_TIMEOUT_S)
        return resp.status_code == 200 and resp.json().get("ok")
    except requests.RequestException:
        return False


def _run_esmfold_via_server(sequence: str, out_pdb: Path) -> str | None:
    import requests
    try:
        resp = requests.post(
            f"{ESMFOLD_SERVER_URL}/fold",
            json={"sequence": sequence, "out_pdb": str(out_pdb), "chunk_size": ESMFOLD_CHUNK_SIZE},
            timeout=ESMFOLD_SERVER_FOLD_TIMEOUT_S,
        )
        payload = resp.json()
    except requests.RequestException as exc:
        return f"ESMFold server request failed: {exc}"
    if not payload.get("ok"):
        return payload.get("error") or "ESMFold server returned an unspecified error"
    return None


def _run_esmfold_via_subprocess(sequence: str, out_pdb: Path) -> str | None:
    chunk_arg = str(ESMFOLD_CHUNK_SIZE) if ESMFOLD_CHUNK_SIZE else "none"
    cmd = [ESMFOLD_ENV_PYTHON, str(_ESMFOLD_SCRIPT), sequence, str(out_pdb), chunk_arg]
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": CUDA_VISIBLE_DEVICES}   # pin GPU0, same as _run_colabfold
    try:
        result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=COLABFOLD_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return "ESMFold prediction timed out"
    if result.returncode != 0:
        _log(f"ESMFold cross-check failed: {result.stderr[-2000:]}")
        return result.stderr[-2000:] or "ESMFold prediction failed (no stderr captured)"
    return None


def _run_esmfold(sequence: str, out_pdb: Path) -> str | None:
    """Returns None on success, an error string on failure -- ESMFold is a
    cross-check, not a hard dependency, so a failure here must not abort
    the alt-structure prediction as a whole (see AltStructureResult.esmfold_error).

    Tries the warm server (scripts/esmfold_server.py) first -- avoids
    reloading the ~2GB model per call (measured ~45s of pure overhead via
    the subprocess path) -- and transparently falls back to the
    subprocess-per-call script if the server isn't running, so starting
    the server is an optional speedup, not a hard dependency.
    """
    if _esmfold_server_available():
        return _run_esmfold_via_server(sequence, out_pdb)
    return _run_esmfold_via_subprocess(sequence, out_pdb)
    return None

#!/usr/bin/env python
"""Run ESMFold on one sequence -- MUST be invoked with the dedicated
`esmfold` conda env's interpreter (config.ESMFOLD_ENV_PYTHON), never with
oneash_dtu: ESMFold's torch/transformers stack lives in an isolated env so
it can't be imported directly from m2_structures.py's process, the same way
junior_surveyor shells out to a separate env's hmmscan/TransDecoder binaries
(see junior_surveyor/config.py's HMMSCAN_BIN/TRANSDECODER_* absolute paths).

Uses transformers' facebook/esmfold_v1 port, not fair-esm+openfold: the HF
port vendors the folding trunk in pure PyTorch, needing no separate openfold
pip package or CUDA-kernel compilation -- much more reliable to install on a
shared box. `model.infer_pdb()` returns a ready PDB string with per-residue
pLDDT already embedded as the B-factor column (same convention AlphaFold2/
ColabFold PDBs use), so m2b_structure_qc.py can read confidence uniformly
from any of these PDBs via Bio.PDB, regardless of which engine produced it.

Usage:
    esmfold_predict.py <sequence> <out_pdb_path> [chunk_size]

Exit code 0 + PDB written on success. Exit code 1 + message on stderr on
failure (including CUDA OOM) -- callers should catch this and skip the
ESMFold cross-check for that hit rather than fail the whole export.
"""

from __future__ import annotations

import sys
from pathlib import Path

MODEL_NAME = "facebook/esmfold_v1"


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: esmfold_predict.py <sequence> <out_pdb_path> [chunk_size]", file=sys.stderr)
        return 1

    sequence = sys.argv[1]
    out_path = Path(sys.argv[2])
    chunk_size = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "none" else None

    import torch
    from transformers import EsmForProteinFolding

    try:
        model = EsmForProteinFolding.from_pretrained(MODEL_NAME, low_cpu_mem_usage=True)
        model = model.cuda()
        # Half-precision LM backbone (HuggingFace's documented pattern) -- in fp32
        # the backbone alone used ~21GB of this box's 24GB, OOMing on anything but
        # short sequences. Trunk/structure module stays fp32 for stability.
        model.esm = model.esm.half()
        model.trunk.set_chunk_size(chunk_size)   # None restores full (fast, more memory)

        with torch.no_grad():
            pdb_str = model.infer_pdb(sequence)
    except torch.cuda.OutOfMemoryError as exc:
        print(f"CUDA OOM predicting sequence of length {len(sequence)}: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001 -- surfaced to caller's subprocess.run, not swallowed
        print(f"ESMFold prediction failed: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(pdb_str)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Retry the 7 sequences from run_colabfold_drug_repurposing.py that timed
out at the standard COLABFOLD_TIMEOUT_S=14400s (4h) -- log:
outputs/master_surveyor/colabfold_drug_repurposing.log, "canonical done:
165 ok ... 3 failed" / "alt done: 165 ok, 4 failed". All 7 are unusually
large proteins (canonical 4547-4967 aa, alt 2446-3641 aa); PLEC and
DYNC1H1 timed out in both the canonical and alt passes, so this is a size
problem, not a fluke. Retries with COLABFOLD_TIMEOUT_S monkeypatched up to
28800s (8h) per sequence -- same colabfold args/GPU otherwise.

Canonical sequences matched by uniprot_acc (unambiguous, get_canonical_structure's
own cache key). Alt sequences matched by seq_hash against the fixed hash
list below (from the FAILED lines in colabfold_drug_repurposing.log).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from master_surveyor import m2_structures
from master_surveyor.config import HITS_CSV

RETRY_TIMEOUT_S = 28800  # 8h, up from the default 14400s (4h)
m2_structures.COLABFOLD_TIMEOUT_S = RETRY_TIMEOUT_S

# PLEC (4547 aa) timed out at 28800s partway through model 2/5 -- colabfold's
# default 5-model canonical ensemble runs ~70min/recycle for proteins this
# large (model 1 alone took 16753s for 3 recycles), so 5 models would need
# ~17-20h, far past any reasonable timeout. Pin to 1 model, matching the
# alt/isoform protocol's own --num-models 1, so a single model fits inside
# RETRY_TIMEOUT_S instead of guaranteeing another timeout.
m2_structures.COLABFOLD_CANONICAL_ARGS = ["--num-models", "1"]

FAILED_CANONICAL_ACCS = {
    "Q15149": "PLEC (4547 aa, timeout)",
    "Q92736": "RYR2 (4967 aa, timeout)",
    "Q14204": "DYNC1H1 (4646 aa, timeout)",
}

FAILED_ALT_HASHES = {
    "d123d8079590d039": "RIF1 rank0 (2446 aa, timeout)",
    "c5a2b56bc1dcec1a": "PLEC rank0 (3447 aa, timeout)",
    "02e767cbec1f8e57": "BPTF rank0 (2457 aa, timeout)",
    "f6cef3489741a3b7": "DYNC1H1 rank0 (3641 aa, timeout)",
}


def _log(msg: str) -> None:
    print(f"[run_colabfold_giant_retry] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    df = pd.read_csv(HITS_CSV)
    dr = df[df["master_group"] == "drug_repurposing_candidate"]

    _log(f"COLABFOLD_TIMEOUT_S raised to {RETRY_TIMEOUT_S}s ({RETRY_TIMEOUT_S / 3600:.0f}h) for this retry")

    canon = dr[["gene_name", "uniprot_acc", "canonical_protein_seq"]].drop_duplicates(subset=["canonical_protein_seq"])
    canon = canon[canon["uniprot_acc"].isin(FAILED_CANONICAL_ACCS)]
    _log(f"retrying {len(canon)}/{len(FAILED_CANONICAL_ACCS)} previously failed canonical sequences")

    canon_done, canon_failed = 0, 0
    for i, (_, row) in enumerate(canon.iterrows(), 1):
        gene, acc = row["gene_name"], row["uniprot_acc"]
        t0 = time.monotonic()
        try:
            result = m2_structures.get_canonical_structure(acc, row["canonical_protein_seq"])
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(canon)}] {gene} ({acc}) ok in {elapsed:.0f}s -- source={result.source}")
            canon_done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(canon)}] {gene} ({acc}) FAILED after {elapsed:.0f}s: {exc}")
            canon_failed += 1
    _log(f"canonical done: {canon_done:,} ok, {canon_failed:,} failed")

    alts = dr[["gene_name", "alt_rank", "alt_protein_seq"]].drop_duplicates(subset=["alt_protein_seq"])
    alts = alts[alts["alt_protein_seq"].apply(lambda s: m2_structures.seq_hash(s) in FAILED_ALT_HASHES)]

    found_hashes = {m2_structures.seq_hash(s) for s in alts["alt_protein_seq"]}
    missing = FAILED_ALT_HASHES.keys() - found_hashes
    if missing:
        _log(f"WARNING: {len(missing)} expected failed alt hash(es) not found in HITS_CSV: {sorted(missing)}")
    _log(f"retrying {len(alts)}/{len(FAILED_ALT_HASHES)} previously failed alt sequences")

    alt_done, alt_failed = 0, 0
    for i, (_, row) in enumerate(alts.iterrows(), 1):
        seq = row["alt_protein_seq"]
        gene = row["gene_name"]
        rank = row["alt_rank"]
        h = m2_structures.seq_hash(seq)
        t0 = time.monotonic()
        try:
            result = m2_structures.get_alt_structure(seq)
            elapsed = time.monotonic() - t0
            _log(
                f"  [{i}/{len(alts)}] {gene} rank{rank} ({h}, {len(seq)} aa) ok in {elapsed:.0f}s -- "
                f"{len(result.seed_models)} seed models, esmfold_error={result.esmfold_error}"
            )
            alt_done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(alts)}] {gene} rank{rank} ({h}, {len(seq)} aa) FAILED after {elapsed:.0f}s: {exc}")
            alt_failed += 1
    _log(f"alt done: {alt_done:,} ok, {alt_failed:,} failed")


if __name__ == "__main__":
    main()

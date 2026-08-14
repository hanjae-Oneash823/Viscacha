"""Fold the trial_failure_candidate structures that have never been attempted
(not failures -- simply never run): 14 canonical sequences (AFDB miss +
never folded locally) and 46 alt/isoform sequences, identified by checking
STRUCTURE_CACHE_DIR directly against get_canonical_structure/get_alt_structure's
own cache-hit conditions (canonical_afdb.pdb, canonical_colabfold's
rank_001 glob, alt_colabfold's rank_* glob).

Giant-sequence handling (same lesson as run_colabfold_giant_retry.py):
BIRC6 (canonical 4857aa; alt ranks 1/3/4 at 4857/4837/3842aa) is larger
than PLEC (4547aa, which needed ~16753s/~4.7h for a single model at the
default 3 recycles). ACACA alt ranks 1/3 (2383/2288aa) already failed once
under the standard alt protocol even with an extended 4h timeout (see
run_colabfold_failed_retry.py's FAILED_HASHES). DOCK10 alt ranks 1/2
(~2180aa) are in the same range.

The alt protocol (COLABFOLD_ALT_ARGS: --num-recycle 20 --num-seeds 5) is
105 forward passes per sequence -- at PLEC's ~70min/recycle for proteins
this large, that is days, not hours. For these three genes' alt
sequences, downgrade to --num-recycle 3 --num-seeds 1 --num-models 1
(matching the giant-canonical protocol's cost) and raise the timeout to
28800s (8h), same as run_colabfold_giant_retry.py. BIRC6's canonical is
also pinned to --num-models 1 with the 8h timeout for the same reason
PLEC needed it. Everything else runs under standard config.py defaults.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from master_surveyor import config, m2_structures
from master_surveyor.config import HITS_CSV, STRUCTURE_CACHE_DIR

GIANT_TIMEOUT_S = 28800  # 8h, up from the default 14400s (4h)
GIANT_CANONICAL_ARGS = ["--num-models", "1"]
GIANT_ALT_ARGS = ["--num-recycle", "3", "--num-seeds", "1", "--num-models", "1"]

GIANT_CANONICAL_GENES = {"BIRC6"}
GIANT_ALT_GENES = {"BIRC6", "ACACA", "DOCK10"}


def _log(msg: str) -> None:
    print(f"[run_colabfold_trial_failure_missing] {msg}", file=sys.stderr, flush=True)


def _canonical_folded(seq: str) -> bool:
    h = m2_structures.seq_hash(seq)
    out_dir = STRUCTURE_CACHE_DIR / h
    if (out_dir / "canonical_afdb.pdb").exists():
        return True
    cf_dir = out_dir / "canonical_colabfold"
    return bool(cf_dir.exists() and list(cf_dir.glob("*_unrelaxed_rank_001_*.pdb")))


def _alt_folded(seq: str) -> bool:
    h = m2_structures.seq_hash(seq)
    out_dir = STRUCTURE_CACHE_DIR / h / "alt_colabfold"
    return bool(out_dir.exists() and list(out_dir.glob("*_unrelaxed_rank_*.pdb")))


def main() -> None:
    df = pd.read_csv(HITS_CSV)
    tf = df[df["master_group"] == "trial_failure_candidate"]

    canon = tf[["gene_name", "uniprot_acc", "canonical_protein_seq"]].dropna(subset=["canonical_protein_seq"])
    canon = canon.drop_duplicates(subset=["canonical_protein_seq"])
    canon = canon[~canon["canonical_protein_seq"].apply(_canonical_folded)]
    canon = canon.assign(is_giant=canon["gene_name"].isin(GIANT_CANONICAL_GENES)).sort_values("is_giant")

    alts = tf[["gene_name", "alt_rank", "alt_protein_seq"]].dropna(subset=["alt_protein_seq"])
    alts = alts.drop_duplicates(subset=["alt_protein_seq"])
    alts = alts[~alts["alt_protein_seq"].apply(_alt_folded)]
    alts = alts.assign(is_giant=alts["gene_name"].isin(GIANT_ALT_GENES)).sort_values("is_giant")

    _log(f"{len(canon)} missing canonical sequences ({canon['is_giant'].sum()} giant), "
         f"{len(alts)} missing alt sequences ({alts['is_giant'].sum()} giant)")

    canon_done, canon_failed = 0, 0
    for i, (_, row) in enumerate(canon.iterrows(), 1):
        gene, acc, is_giant = row["gene_name"], row["uniprot_acc"], row["is_giant"]
        if is_giant:
            m2_structures.COLABFOLD_CANONICAL_ARGS = GIANT_CANONICAL_ARGS
            m2_structures.COLABFOLD_TIMEOUT_S = GIANT_TIMEOUT_S
        else:
            m2_structures.COLABFOLD_CANONICAL_ARGS = config.COLABFOLD_CANONICAL_ARGS
            m2_structures.COLABFOLD_TIMEOUT_S = config.COLABFOLD_TIMEOUT_S
        t0 = time.monotonic()
        try:
            result = m2_structures.get_canonical_structure(acc, row["canonical_protein_seq"])
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(canon)}] {gene} ({acc}) ok in {elapsed:.0f}s -- source={result.source}"
                 f"{' [giant]' if is_giant else ''}")
            canon_done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(canon)}] {gene} ({acc}) FAILED after {elapsed:.0f}s: {exc}")
            canon_failed += 1
    _log(f"canonical done: {canon_done:,} ok, {canon_failed:,} failed")

    alt_done, alt_failed = 0, 0
    for i, (_, row) in enumerate(alts.iterrows(), 1):
        seq, gene, rank, is_giant = row["alt_protein_seq"], row["gene_name"], row["alt_rank"], row["is_giant"]
        if is_giant:
            m2_structures.COLABFOLD_ALT_ARGS = GIANT_ALT_ARGS
            m2_structures.COLABFOLD_TIMEOUT_S = GIANT_TIMEOUT_S
        else:
            m2_structures.COLABFOLD_ALT_ARGS = config.COLABFOLD_ALT_ARGS
            m2_structures.COLABFOLD_TIMEOUT_S = config.COLABFOLD_TIMEOUT_S
        h = m2_structures.seq_hash(seq)
        t0 = time.monotonic()
        try:
            result = m2_structures.get_alt_structure(seq)
            elapsed = time.monotonic() - t0
            _log(
                f"  [{i}/{len(alts)}] {gene} rank{rank} ({h}, {len(seq)} aa) ok in {elapsed:.0f}s -- "
                f"{len(result.seed_models)} seed models, esmfold_error={result.esmfold_error}"
                f"{' [giant]' if is_giant else ''}"
            )
            alt_done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(alts)}] {gene} rank{rank} ({h}, {len(seq)} aa) FAILED after {elapsed:.0f}s: {exc}")
            alt_failed += 1
    _log(f"alt done: {alt_done:,} ok, {alt_failed:,} failed")


if __name__ == "__main__":
    main()

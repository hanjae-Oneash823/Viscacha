"""ColabFold pass over ALL ranked alt-isoform sequences for the 43
trial_failure hits (114 raw alt_rank rows, 104 unique sequences after
dedup) -- not just the single representative alt per hit that
run_colabfold_alts.py covers. Mirrors run_esmfold_all_alts.py's scope on
the ESMFold side.

No length-based exclusions this time (unlike run_colabfold_alts.py's
LARGE_PROTEIN_EXCLUSIONS): COLABFOLD_TIMEOUT_S (3600s, in config.py) already
caps any single sequence's wall-clock cost via _run_colabfold's
subprocess.run(..., timeout=...), which kills a stalled colabfold_batch
call and lets the batch move on rather than hanging indefinitely -- this is
what let ACACA's earlier 39+ minute stall go unbounded before (it was run
without this loop's per-sequence timeout wrapping in practice), so it's now
safe to include everything.

Uses get_alt_structure() directly (not a raw _run_colabfold call), so
caching, the ESMFold cross-check, and file layout stay identical to what
run_master_surveyor.run_one() would produce -- and any sequence already
folded by run_colabfold_alts.py's earlier partial run is skipped for free.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from master_surveyor import m2_structures
from master_surveyor.config import HITS_CSV


def _log(msg: str) -> None:
    print(f"[run_colabfold_all_alts] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    df = pd.read_csv(HITS_CSV)
    tf_raw = df[df["master_group"] == "trial_failure_candidate"]
    n_hits = tf_raw.groupby(["gene_name", "cell_type", "hit_ENST_ID"]).ngroups
    _log(f"{len(tf_raw):,} raw alt_rank rows across {n_hits} trial_failure hits")

    seqs = tf_raw[["gene_name", "alt_rank", "alt_protein_seq"]].drop_duplicates(subset=["alt_protein_seq"])
    _log(f"{len(seqs):,} unique alt sequences to fold (deduped from {len(tf_raw):,} rows)")

    done, failed = 0, 0
    for i, (_, row) in enumerate(seqs.iterrows(), 1):
        seq = row["alt_protein_seq"]
        gene = row["gene_name"]
        rank = row["alt_rank"]
        t0 = time.monotonic()
        try:
            result = m2_structures.get_alt_structure(seq)
            elapsed = time.monotonic() - t0
            _log(
                f"  [{i}/{len(seqs)}] {gene} rank{rank} ({len(seq)} aa) ok in {elapsed:.0f}s -- "
                f"{len(result.seed_models)} seed models, esmfold_error={result.esmfold_error}"
            )
            done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(seqs)}] {gene} rank{rank} ({len(seq)} aa) FAILED after {elapsed:.0f}s: {exc}")
            failed += 1

    _log(f"done: {done:,} ok, {failed:,} failed")


if __name__ == "__main__":
    main()

"""ESMFold pass over ALL ranked alt-isoform sequences for the 43 trial_failure
hits, not just the single gate-driving alt m0_select.representative_row()
picks per hit. The raw hits_deep.csv is long-format for trial_failure_candidate
(one row per alt_rank) -- 114 rows across the 43 hits, 2-6 ranked alts each.

Reads master_group-scoped rows directly (pre-collapse), skipping m0_select's
representative_row() collapse entirely, so every ranked alt gets folded, not
just the top one. Same seq_hash cache path as run_esmfold_pilot.py, so any alt
sequence already folded (e.g. the 35 gate-driver alts from the prior pilot) is
skipped, and results here are reused by the full pipeline later too.

GPU0 only: relies on config.CUDA_VISIBLE_DEVICES="0", used by both
_run_esmfold_via_server (the already-running warm server, itself launched with
CUDA_VISIBLE_DEVICES=0) and _run_esmfold_via_subprocess's fallback path.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from master_surveyor import m2_structures
from master_surveyor.config import HITS_CSV, STRUCTURE_CACHE_DIR


def _log(msg: str) -> None:
    print(f"[run_esmfold_all_alts] {msg}", file=sys.stderr, flush=True)


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
        h = m2_structures.seq_hash(seq)
        out_pdb = STRUCTURE_CACHE_DIR / h / "esmfold.pdb"

        if out_pdb.exists():
            _log(f"  [{i}/{len(seqs)}] {gene} rank{rank} ({len(seq)} aa) already cached, skipping")
            done += 1
            continue

        t0 = time.monotonic()
        err = m2_structures._run_esmfold(seq, out_pdb)
        elapsed = time.monotonic() - t0

        if err:
            _log(f"  [{i}/{len(seqs)}] {gene} rank{rank} ({len(seq)} aa) FAILED after {elapsed:.0f}s: {err[:300]}")
            failed += 1
        else:
            _log(f"  [{i}/{len(seqs)}] {gene} rank{rank} ({len(seq)} aa) ok in {elapsed:.0f}s")
            done += 1

    _log(f"done: {done:,} ok, {failed:,} failed")


if __name__ == "__main__":
    main()

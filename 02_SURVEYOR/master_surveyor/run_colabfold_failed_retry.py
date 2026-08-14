"""Retry just the 14 alt/isoform sequences that failed in the full
run_colabfold_all_alts.py pass (log: outputs/master_surveyor/colabfold_all_alts.log,
"done: 90 ok, 14 failed"). Reruns on GPU1 (config.CUDA_VISIBLE_DEVICES="1")
with COLABFOLD_TIMEOUT_S raised to 14400s (from 3600s) -- 13 of the 14
failures were hitting the old timeout on long sequences (up to ACACA's
2383 aa); one (CHD1 rank2, hash a5ed87f7cd4a3089) crashed with exit -11
(segfault) after 1877s rather than timing out, so the longer timeout won't
help it directly, but GPU1 currently has more free VRAM which may.

Sequences are looked up by seq_hash(alt_protein_seq) against the fixed
hash list below (from the FAILED lines in colabfold_all_alts.log), not by
(gene, alt_rank), so this is exact regardless of row ordering.

Reuses get_alt_structure() (same caching/file layout as the full run), so
if a hash somehow already has cached rank_* models this is a no-op for it.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from master_surveyor import m2_structures
from master_surveyor.config import HITS_CSV

FAILED_HASHES = {
    "93cfe8189194eba3": "ARID4A rank2 (1257 aa, timeout)",
    "c86266cbe482b3e1": "TUT7 rank2 (1495 aa, timeout)",
    "c04c2215c32a67ad": "CHD1 rank3 (1710 aa, timeout)",
    "b7628d6de8c66aae": "ACACA rank1 (2383 aa, timeout)",
    "3897523e1e43471e": "GAPVD1 rank4 (1460 aa, timeout)",
    "69d3c72ad0146285": "TNIK rank1 (1360 aa, timeout)",
    "2c1a8c5f1398bb51": "TUT7 rank3 (1259 aa, timeout)",
    "a5ed87f7cd4a3089": "CHD1 rank2 (1798 aa, exit -11 segfault)",
    "6ba40f0df65a354c": "ACACA rank3 (2288 aa, timeout)",
    "081b3d3250ae903f": "GAPVD1 rank2 (1433 aa, timeout)",
    "3f79c1b66d961479": "GAPVD1 rank3 (1487 aa, timeout)",
    "ead96aa172d1b132": "GAPVD1 rank5 (1439 aa, timeout)",
    "576a786dfce5983e": "SACM1L rank4 (609 aa, no rank_* models produced)",
    "89b28998ae31b2c8": "TNIK rank2 (1352 aa, timeout)",
}


def _log(msg: str) -> None:
    print(f"[run_colabfold_failed_retry] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    df = pd.read_csv(HITS_CSV)
    tf_raw = df[df["master_group"] == "trial_failure_candidate"]
    seqs = tf_raw[["gene_name", "alt_rank", "alt_protein_seq"]].drop_duplicates(subset=["alt_protein_seq"])
    seqs = seqs[seqs["alt_protein_seq"].apply(lambda s: m2_structures.seq_hash(s) in FAILED_HASHES)]

    found_hashes = {m2_structures.seq_hash(s) for s in seqs["alt_protein_seq"]}
    missing = FAILED_HASHES.keys() - found_hashes
    if missing:
        _log(f"WARNING: {len(missing)} expected failed hash(es) not found in HITS_CSV: {sorted(missing)}")
    _log(f"retrying {len(seqs)}/{len(FAILED_HASHES)} previously failed sequences")

    done, failed = 0, 0
    for i, (_, row) in enumerate(seqs.iterrows(), 1):
        seq = row["alt_protein_seq"]
        gene = row["gene_name"]
        rank = row["alt_rank"]
        h = m2_structures.seq_hash(seq)
        t0 = time.monotonic()
        try:
            result = m2_structures.get_alt_structure(seq)
            elapsed = time.monotonic() - t0
            _log(
                f"  [{i}/{len(seqs)}] {gene} rank{rank} ({h}, {len(seq)} aa) ok in {elapsed:.0f}s -- "
                f"{len(result.seed_models)} seed models, esmfold_error={result.esmfold_error}"
            )
            done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(seqs)}] {gene} rank{rank} ({h}, {len(seq)} aa) FAILED after {elapsed:.0f}s: {exc}")
            failed += 1

    _log(f"done: {done:,} ok, {failed:,} failed")


if __name__ == "__main__":
    main()

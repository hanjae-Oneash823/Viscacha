"""ESMFold pass over the 43 trial_failure hits' CANONICAL sequences --
independent of get_canonical_structure()'s normal AFDB-first/ColabFold-fallback
path. Gives the viewer an ESMFold option for canonical structures without
waiting on AFDB lookups or a ColabFold run, same rationale as
run_esmfold_all_alts.py for the alt side: ESMFold is MSA-free and
non-iterative, so it's a fast way to get *something* viewable while the
slower/more authoritative ColabFold pipeline hasn't run yet for most hits.

Writes to the same STRUCTURE_CACHE_DIR/<seq_hash>/esmfold.pdb path keyed by
sequence content, so this is never wasted work and doesn't collide with the
alt-isoform ESMFold cache (canonical and alt sequences hash differently).

GPU0 only: relies on config.CUDA_VISIBLE_DEVICES="0", used by both
_run_esmfold_via_server (the warm server, itself launched with
CUDA_VISIBLE_DEVICES=0) and _run_esmfold_via_subprocess's fallback path.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor import m0_select, m2_structures
from master_surveyor.config import STRUCTURE_CACHE_DIR


def _log(msg: str) -> None:
    print(f"[run_esmfold_canonical] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    shortlist = m0_select.run()
    tf = shortlist[shortlist["master_group"] == "trial_failure_candidate"]
    _log(f"{len(tf):,} trial_failure hits")

    seqs = tf[["gene_name", "canonical_protein_seq"]].drop_duplicates(subset=["canonical_protein_seq"])
    _log(f"{len(seqs):,} unique canonical sequences to fold (deduped from {len(tf):,} hits)")

    done, failed = 0, 0
    for i, (_, row) in enumerate(seqs.iterrows(), 1):
        seq = row["canonical_protein_seq"]
        gene = row["gene_name"]
        h = m2_structures.seq_hash(seq)
        out_pdb = STRUCTURE_CACHE_DIR / h / "esmfold.pdb"

        if out_pdb.exists():
            _log(f"  [{i}/{len(seqs)}] {gene} ({len(seq)} aa) already cached, skipping")
            done += 1
            continue

        t0 = time.monotonic()
        err = m2_structures._run_esmfold(seq, out_pdb)
        elapsed = time.monotonic() - t0

        if err:
            _log(f"  [{i}/{len(seqs)}] {gene} ({len(seq)} aa) FAILED after {elapsed:.0f}s: {err[:300]}")
            failed += 1
        else:
            _log(f"  [{i}/{len(seqs)}] {gene} ({len(seq)} aa) ok in {elapsed:.0f}s")
            done += 1

    _log(f"done: {done:,} ok, {failed:,} failed")


if __name__ == "__main__":
    main()

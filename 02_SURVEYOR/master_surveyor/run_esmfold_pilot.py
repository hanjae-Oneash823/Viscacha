"""ESMFold-only pass over the trial_failure pilot batch's alt sequences --
decoupled from ColabFold specifically because ESMFold is MSA-free and
non-iterative (one forward pass, no recycling), so it shouldn't be nearly as
sensitive to sequence length as ColabFold's recycled search (which stalled
39+ minutes on ACACA's 2288-residue alt sequence without completing).

Writes results to the exact cache path get_alt_structure() already checks
before invoking ESMFold, so this is never wasted work: the full pipeline
will find these cached and skip re-folding them.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor import m0_select, m2_structures
from master_surveyor.config import STRUCTURE_CACHE_DIR


def _log(msg: str) -> None:
    print(f"[run_esmfold_pilot] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    shortlist = m0_select.run()
    tf = shortlist[shortlist["master_group"] == "trial_failure_candidate"]
    _log(f"{len(tf):,} trial_failure hits")

    seqs = tf[["gene_name", "cell_type", "alt_protein_seq"]].drop_duplicates(subset=["alt_protein_seq"])
    _log(f"{len(seqs):,} unique alt sequences to fold (deduped from {len(tf):,} hits)")

    done, failed = 0, 0
    for i, (_, row) in enumerate(seqs.iterrows(), 1):
        seq = row["alt_protein_seq"]
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

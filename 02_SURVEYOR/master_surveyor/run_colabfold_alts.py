"""ColabFold pass over the 43 trial_failure hits' representative alt/isoform
sequences -- the template-free, elevated-recycle, multi-seed protocol
documented in m2_structures.py's module docstring (COLABFOLD_ALT_ARGS: 20
recycles x 5 seeds). Complements the earlier ESMFold-only pre-pass already
cached for these same sequences: this is the slower, more authoritative fold
get_alt_structure() needs for a real export, and it reuses that function
directly so caching, the ESMFold cross-check, and file layout all stay
identical to what run_master_surveyor.run_one() would produce.

Same LARGE_PROTEIN_EXCLUSIONS as run_pilot_tf_filtered.py: ACACA's alt
sequence alone stalled ColabFold 39+ minutes without completing earlier this
session; CHD1/GAPVD1/TNIK are excluded on the same signal (all also failed
ESMFold's much cheaper single-pass fold at similar lengths, >1250 aa). These
4 need a dedicated strategy (windowing, most likely) and shouldn't block
getting real ColabFold results for the other 35 first.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor import m0_select, m2_structures

LARGE_PROTEIN_EXCLUSIONS = {"ACACA", "CHD1", "GAPVD1", "TNIK"}


def _log(msg: str) -> None:
    print(f"[run_colabfold_alts] {msg}", file=sys.stderr, flush=True)


def main() -> None:
    shortlist = m0_select.run()
    tf = shortlist[shortlist["master_group"] == "trial_failure_candidate"]
    tf = tf[~tf["gene_name"].isin(LARGE_PROTEIN_EXCLUSIONS)]

    seqs = tf[["gene_name", "alt_protein_seq"]].drop_duplicates(subset=["alt_protein_seq"])
    _log(f"{len(seqs):,} unique alt sequences to fold (excluding {sorted(LARGE_PROTEIN_EXCLUSIONS)})")

    done, failed = 0, 0
    for i, (_, row) in enumerate(seqs.iterrows(), 1):
        seq = row["alt_protein_seq"]
        gene = row["gene_name"]
        t0 = time.monotonic()
        try:
            result = m2_structures.get_alt_structure(seq)
            elapsed = time.monotonic() - t0
            _log(
                f"  [{i}/{len(seqs)}] {gene} ({len(seq)} aa) ok in {elapsed:.0f}s -- "
                f"{len(result.seed_models)} seed models, esmfold_error={result.esmfold_error}"
            )
            done += 1
        except Exception as exc:
            elapsed = time.monotonic() - t0
            _log(f"  [{i}/{len(seqs)}] {gene} ({len(seq)} aa) FAILED after {elapsed:.0f}s: {exc}")
            failed += 1

    _log(f"done: {done:,} ok, {failed:,} failed")


if __name__ == "__main__":
    main()

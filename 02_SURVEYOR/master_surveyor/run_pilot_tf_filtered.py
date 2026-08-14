"""Pilot batch: trial_failure_candidate hits, EXCLUDING the 4 large-protein
outliers already known to exceed this GPU's memory (ACACA 2288aa, CHD1
1798aa, GAPVD1 1439aa, TNIK 1352aa -- confirmed via the ESMFold pre-pass;
ACACA separately stalled ColabFold for 39+ min without completing). These 4
need a dedicated strategy (windowing around the altered region, most likely)
and shouldn't block getting clean results for the other 39 hits first.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor import m0_select
from master_surveyor.run_master_surveyor import _log, run_one

LARGE_PROTEIN_EXCLUSIONS = {"ACACA", "CHD1", "GAPVD1", "TNIK"}


def main() -> None:
    shortlist = m0_select.run()
    tf = shortlist[shortlist["master_group"] == "trial_failure_candidate"]
    tf = tf[~tf["gene_name"].isin(LARGE_PROTEIN_EXCLUSIONS)]
    _log(f"pilot batch: {len(tf):,} trial_failure hits (excluding {sorted(LARGE_PROTEIN_EXCLUSIONS)})")

    exported, failed = 0, 0
    for i, (_, row) in enumerate(tf.iterrows(), 1):
        t0 = time.monotonic()
        result = run_one(row.to_dict(), thresholds_used={})
        elapsed = time.monotonic() - t0
        if result is not None:
            exported += 1
        else:
            failed += 1
        _log(f"  [{i}/{len(tf)}] {row['gene_name']} {elapsed:.0f}s elapsed  (running total: {exported} ok, {failed} failed)")

    _log(f"pilot done: {exported:,} exported, {failed:,} failed")


if __name__ == "__main__":
    main()

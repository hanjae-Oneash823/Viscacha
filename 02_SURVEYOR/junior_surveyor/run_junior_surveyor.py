"""JUNIOR_SURVEYOR — main entry point.

Splits junior_pass hits by candidate_group into two branch-specific analyses,
then runs shared J2 -> J3 -> J4 on the merged long-format result:

  new_target_candidate    -> j1_sequences   (alt_rank=0: DTU transcript vs
                              MANE canonical -- does DTU reveal a novel target?)
  trial_failure_candidate -> j1b_isoform_ranking (alt_rank=1..N: every isoform
                              of the gene ranked by AD-donor usage, all diffed
                              against MANE canonical -- did a poor-drug-partner
                              alternate take over when the canonical form's
                              usage collapsed in AD?)

Writes outputs/junior_surveyor/hits_deep.csv (long format: one row per
alt_rank; new_target hits are exactly one row, trial_failure hits are one row
per ranked isoform of that gene, no usage cutoff applied).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Ensure repo root (for classify_hit_scenarios_mane, imported by j1_canonical)
# and 02_SURVEYOR (this package's parent, for the self-referential
# `junior_surveyor.x` imports below) are on path regardless of invocation cwd.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from junior_surveyor.config import CACHE_DIR, HITS_CSV, OUT_DIR
from junior_surveyor import (
    j1_sequences, j1b_isoform_ranking, j1c_alt_biotype,
    j2_protein_diff, j3_drug_targets, j4_gate,
)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load ────────────────────────────────────────────────────────────────
    print("Loading hits_enriched.csv …")
    all_hits = pd.read_csv(HITS_CSV)
    hits = all_hits[all_hits["junior_pass"]].copy()
    print(f"  {len(hits):,} junior_pass hits  |  "
          f"{hits['gene_name'].nunique():,} genes\n")

    new_target_hits    = hits[hits["candidate_group"] == "new_target_candidate"]
    trial_failure_hits = hits[hits["candidate_group"] == "trial_failure_candidate"]
    print(f"  new_target_candidate:    {len(new_target_hits):,} hits\n"
          f"  trial_failure_candidate: {len(trial_failure_hits):,} hits\n")

    # ── J1 — branch-specific sequence resolution ──────────────────────────────
    print("J1 — new_target: DTU transcript vs MANE canonical …")
    nt_df = j1_sequences.run(new_target_hits)
    print()

    print("J1b — trial_failure: ranking every isoform by AD-donor usage …")
    tf_df = j1b_isoform_ranking.run(trial_failure_hits)
    print()

    df = pd.concat([nt_df, tf_df], ignore_index=True, sort=False)
    print(f"  merged long-format frame: {len(df):,} rows "
          f"({df['gene_name'].nunique():,} genes)\n")

    # ── J1c — per-row alt biotype classification + quality filter ──────────────
    print("J1c — classifying alt_biotype_class + filtering low-quality rows …")
    df = j1c_alt_biotype.run(df)
    print()

    # ── J2: protein diff ────────────────────────────────────────────────────
    print("J2 — aligning canonical vs alternative proteins …")
    df = j2_protein_diff.run(df)
    print()

    # ── J3: drug targets ────────────────────────────────────────────────────
    print("J3 — querying ChEMBL + DGIdb + Open Targets + Pharos …")
    df = j3_drug_targets.run(df)
    print()

    # ── J4: unified gate ─────────────────────────────────────────────────────
    print("J4 — gating hits for master_surveyor …")
    df = j4_gate.run(df)
    print()

    # ── Save ────────────────────────────────────────────────────────────────
    out_csv = OUT_DIR / "hits_deep.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved → {out_csv}")
    print(f"  {len(df):,} rows  |  "
          f"{df['gene_name'].nunique():,} genes  |  "
          f"{df.shape[1]} columns")


if __name__ == "__main__":
    main()

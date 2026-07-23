#!/usr/bin/env python3
"""
INITIAL_FILTER — runs before ASSISTANT_SURVEYOR.

Filters the 2,599 permutation-significant DTU hits down to the two
highest-confidence candidate groups, combining empirical isoform dominance
(Control-PSI rank per gene x cell type, from classify_hit_scenarios.py) with
MANE Select transcript role (classify_hit_scenarios_mane.py):

  trial_failure_candidate = Dominant isoform, MANE Canonical, CT_enriched
                            (the reference/major transcript loses usage in AD
                            -- may explain why a drug targeting it failed)
  new_target_candidate    = Minor isoform, MANE Alternate, AD_enriched
                            (a non-reference transcript gains usage in AD
                            -- a candidate not previously drugged)

The other two role x MANE combinations (Dominant x Alternate, Minor x
Canonical, regardless of direction) are labeled "other"; anything with no
MANE coverage is labeled "no_MANE_coverage". Both are dead ends -- they are
NOT passed to ASSISTANT_SURVEYOR -- but are kept (with full 2,599-row
coverage) in ALL_GROUPS_OUT_CSV purely so the junior_surveyor Sankey can show
where the excluded hits went.

~1,200 of the 2,599 hits pass. ASSISTANT_SURVEYOR reads OUT_CSV, not the
full DIU_significant_hits_combined.csv -- see assistant_surveyor/config.py:HITS_CSV.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python initial_filter.py
"""

from pathlib import Path

import pandas as pd

from classify_hit_scenarios import load_dominance
from classify_hit_scenarios_mane import load_mane_select

REPO_ROOT         = Path(__file__).resolve().parent
HITS_CSV          = REPO_ROOT / "outputs/layer1_1/DIU_significant_hits_combined.csv"
TX_ID_MAP         = REPO_ROOT / "outputs/layer1_1/annotation/tx_id_map.csv"
OUT_CSV           = REPO_ROOT / "outputs/layer1_1/DIU_significant_hits_initial_filter.csv"
ALL_GROUPS_OUT_CSV = REPO_ROOT / "outputs/layer1_1/DIU_significant_hits_all_candidate_groups.csv"


def _tx_role(row) -> str:
    if pd.isna(row["ENST_ID"]) or pd.isna(row["mane_ENST_ID"]):
        return "no_MANE_coverage"
    return "Canonical" if row["ENST_ID"] == row["mane_ENST_ID"] else "Alternate"


def main() -> None:
    hits = pd.read_csv(HITS_CSV)
    n_hits = len(hits)

    dom = load_dominance()
    hits = hits.merge(dom, on=["transcript_name", "gene_name", "cell_type"], how="left")

    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID", "ENSG_ID", "gene_name"]]
    hits = hits.merge(id_map, on=["transcript_name", "gene_name"], how="left")

    mane = load_mane_select()
    hits = hits.merge(mane, on="ENSG_ID", how="left")
    hits["tx_role_mane"] = hits.apply(_tx_role, axis=1)
    assert len(hits) == n_hits, "row count changed during merge — check for duplicate keys"

    is_trial_failure = (
        (hits["isoform_role"] == "Dominant")
        & (hits["tx_role_mane"] == "Canonical")
        & (hits["usage_direction"] == "CT_enriched")
    )
    is_new_target = (
        (hits["isoform_role"] == "Minor")
        & (hits["tx_role_mane"] == "Alternate")
        & (hits["usage_direction"] == "AD_enriched")
    )

    is_no_mane = hits["tx_role_mane"] == "no_MANE_coverage"

    hits["candidate_group"] = "other"
    hits.loc[is_no_mane, "candidate_group"] = "no_MANE_coverage"
    hits.loc[is_trial_failure, "candidate_group"] = "trial_failure_candidate"
    hits.loc[is_new_target, "candidate_group"] = "new_target_candidate"

    # Drop the join scaffolding used only to compute tx_role_mane -- these
    # collide (ENST_ID/ENSG_ID) with columns assistant_surveyor's own L1
    # stage merges in from the same tx_id_map.csv further down the pipeline.
    hits_clean = hits.drop(columns=["ENST_ID", "ENSG_ID", "mane_ENST_ID", "symbol"])

    filtered = hits_clean[hits_clean["candidate_group"].isin(
        ["trial_failure_candidate", "new_target_candidate"]
    )].copy()

    print(f"Input hits: {n_hits}")
    print(hits_clean["candidate_group"].value_counts().to_string())
    print(f"Total passing initial_filter: {len(filtered)} ({len(filtered) / n_hits:.1%} of input)")
    print(f"Unique genes (passing): {filtered['gene_name'].nunique()}")

    filtered.to_csv(OUT_CSV, index=False)
    print(f"Saved -> {OUT_CSV}")

    hits_clean.to_csv(ALL_GROUPS_OUT_CSV, index=False)
    print(f"Saved -> {ALL_GROUPS_OUT_CSV} (all {n_hits} hits, 4-way candidate_group, for Sankey use)")


if __name__ == "__main__":
    main()

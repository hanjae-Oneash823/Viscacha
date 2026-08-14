"""M0 — final filter + collapse to one row per hit.

Restricts junior_surveyor's hits_deep.csv to master_surveyor's own scope
(trial_failure_candidate + drug_repurposing_candidate, per
docs/MASTER_SURVEYOR_plan.md), collapses trial_failure's long-format
multi-alt_rank rows to one representative row per hit, then applies
group-specific filter predicates.

Unlike the plan doc's original sketch, thresholds are a **parameter**, not
fixed module constants -- dossier_server passes live UI values straight
through with no code change. Defaults (config.DEFAULT_THRESHOLDS) are the
same deliberately permissive placeholders: nothing is dropped until the
user tightens a cutoff.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor.config import DEFAULT_THRESHOLDS, HITS_CSV, MASTER_SURVEYOR_GROUPS


def _log(msg: str) -> None:
    print(f"  [m0] {msg}", file=sys.stderr, flush=True)


def representative_row(g: pd.DataFrame) -> pd.Series:
    """One row per hit.

    trial_failure_candidate hits are long-format (one row per ranked
    alternate isoform) -- collapse to the gate-driving alt with the largest
    AD usage gain. drug_repurposing_candidate (and any other group) is
    already one row per hit, first row is representative.

    This is the single implementation of logic that used to be duplicated
    three ways: master_surveyor/plot_results.py's _tf_representative_row,
    dossier/generate_index.py's _representative_row, and (implicitly) this
    module's own m0. Both existing callers now import this function instead.
    """
    if g["master_group"].iloc[0] == "trial_failure_candidate":
        drivers = g[g["alt_is_gate_driver"]]
        pool = drivers if not drivers.empty else g
        return pool.loc[pool["alt_usage_delta"].idxmax()]
    return g.iloc[0]


def _collapse(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per (gene_name, cell_type, hit_ENST_ID).

    hit_ENST_ID (not alt_ENST_ID) is the stable hit identity: for
    trial_failure's long-format rows it stays constant across every ranked
    alt_rank row of the same hit (so grouping by it correctly collapses the
    whole ranked series down to one representative row), while still
    keeping genuinely distinct hits on the same (gene_name, cell_type) pair
    (e.g. two different DTU-significant transcripts of the same gene/cell
    type) separate, since those carry two different hit_ENST_ID values.
    alt_ENST_ID varies per alt_rank row WITHIN a single trial_failure hit,
    so grouping by it instead would silently defeat the collapse (verified
    empirically: it produced 99 "trial_failure hits" instead of a real
    collapse). Matches dossier/generate_index.py's build_manifest(), which
    already groups this same way.
    """
    rows = [
        representative_row(g)
        for _, g in df.groupby(["gene_name", "cell_type", "hit_ENST_ID"], dropna=False)
    ]
    return pd.DataFrame(rows).reset_index(drop=True)


def _apply_thresholds(df: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
    kept = df
    n0 = len(kept)

    tf_mask = kept["master_group"] == "trial_failure_candidate"
    dr_mask = kept["master_group"] == "drug_repurposing_candidate"

    min_delta = thresholds.get("tf_min_abs_delta_usage", 0.0)
    if min_delta and min_delta > 0:
        drop = tf_mask & (kept["alt_usage_delta"].abs() < min_delta)
        kept = kept[~drop]
        _log(f"kept {len(kept):,}/{n0:,} after tf_min_abs_delta_usage>={min_delta}")

    if thresholds.get("tf_require_domain_overlap"):
        n1 = len(kept)
        tf_mask = kept["master_group"] == "trial_failure_candidate"
        no_domain = kept["affected_domain"].isna() | (kept["affected_domain"] == "none") \
            | (kept["affected_domain"].astype(str).str.strip() == "")
        drop = tf_mask & no_domain
        kept = kept[~drop]
        _log(f"kept {len(kept):,}/{n1:,} after tf_require_domain_overlap")

    min_phase = thresholds.get("dr_min_chembl_or_ot_phase", 0)
    if min_phase and min_phase > 0:
        n1 = len(kept)
        dr_mask = kept["master_group"] == "drug_repurposing_candidate"
        best_phase = kept[["chembl_max_phase", "ot_max_phase"]].fillna(0).max(axis=1)
        drop = dr_mask & (best_phase < min_phase)
        kept = kept[~drop]
        _log(f"kept {len(kept):,}/{n1:,} after dr_min_chembl_or_ot_phase>={min_phase}")

    if thresholds.get("dr_require_structural_change"):
        n1 = len(kept)
        dr_mask = kept["master_group"] == "drug_repurposing_candidate"
        drop = dr_mask & (kept["protein_change_type"] == "substitution")
        kept = kept[~drop]
        _log(f"kept {len(kept):,}/{n1:,} after dr_require_structural_change")

    return kept


def run(df: pd.DataFrame | None = None, thresholds: dict | None = None) -> pd.DataFrame:
    """df defaults to a fresh read of HITS_CSV; thresholds defaults to the
    permissive DEFAULT_THRESHOLDS. Pass a partial dict -- missing keys fall
    back to the default for that key.
    """
    if df is None:
        df = pd.read_csv(HITS_CSV)
    merged_thresholds = {**DEFAULT_THRESHOLDS, **(thresholds or {})}

    scoped = df[df["master_group"].isin(MASTER_SURVEYOR_GROUPS)].copy()
    _log(f"{len(scoped):,}/{len(df):,} rows in scope ({MASTER_SURVEYOR_GROUPS})")

    collapsed = _collapse(scoped)
    n_tf = int((collapsed["master_group"] == "trial_failure_candidate").sum())
    n_dr = int((collapsed["master_group"] == "drug_repurposing_candidate").sum())
    _log(f"collapsed to {len(collapsed):,} hits ({n_tf:,} trial_failure, {n_dr:,} drug_repurposing)")

    return _apply_thresholds(collapsed, merged_thresholds)


if __name__ == "__main__":
    result = run()
    print(f"\n{len(result):,} hits after m0_select (default thresholds)")
    dup_check = result.groupby(["gene_name", "cell_type"])["hit_ENST_ID"].nunique()
    multi = dup_check[dup_check > 1]
    if not multi.empty:
        print(f"{len(multi)} (gene, cell_type) pairs carry >1 distinct hit:")
        print(multi)

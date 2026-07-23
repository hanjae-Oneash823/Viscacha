"""J4 — Selection for the next surveyor stage.

A hit is selected_for_next_stage when BOTH hold:
  - protein_change_type is a real protein-level change (not identical/no_sequence)
  - existing drug/target evidence (chembl_max_phase >= 1 OR dgidb_interactions > 0)
"""

from __future__ import annotations

import pandas as pd

from junior_surveyor.config import NULL_PROTEIN_CHANGE_TYPES


def _has_protein_change(row: pd.Series) -> bool:
    return row.get("protein_change_type", "no_sequence") not in NULL_PROTEIN_CHANGE_TYPES


def _has_drug_evidence(row: pd.Series) -> bool:
    phase = int(row.get("chembl_max_phase", 0) or 0)
    dgidb = int(row.get("dgidb_interactions", 0) or 0)
    return phase >= 1 or dgidb > 0


def _select(row: pd.Series) -> bool:
    return _has_protein_change(row) and _has_drug_evidence(row)


def run(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["selected_for_next_stage"] = df.apply(_select, axis=1)

    n_sel = int(df["selected_for_next_stage"].sum())
    n_genes = df.loc[df["selected_for_next_stage"], "gene_name"].nunique()
    print(f"  [j4] selected_for_next_stage: {n_sel:,} / {len(df):,} hits "
          f"({n_genes:,} genes)")
    return df

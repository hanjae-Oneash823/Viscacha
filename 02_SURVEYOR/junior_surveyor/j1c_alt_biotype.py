"""J1C — per-row biotype classification + quality filtering for ranked alternates.

Every J1/J1b row carries an alt_transcript_name/alt_ENST_ID -- the isoform
actually being compared to canonical in that row -- but the biotype_class
column inherited from assistant_surveyor's L1 describes the ORIGINAL hit
transcript, not this row's alt. For new_target rows alt IS the hit, so that
label is already correct; for trial_failure's ranked-alternate rows (all but
the one row where the alt happens to be canonical itself) it isn't -- most
never went through L1 at all, since J1b introduces isoforms L1 never saw.

Two groups of rows are passed through untouched, without classification or
filtering, because the answer is already known or the filter can't change
the outcome:

  - new_target rows: alt_biotype_class is just copied from the existing
    biotype_class (re-deriving it would reproduce the same value, since that
    biotype was already computed against this same alt transcript and
    already passed this same 5-category filter in assistant_surveyor's
    junior_pass gate). The 5% AD-usage filter never applied to new_target
    either (single alt_rank=0 row, already a significant DTU hit by
    construction).

  - trial_failure's is_canonical rows (exactly one per hit -- the gene's own
    canonical transcript reappearing in its own J1b ranking): always
    alt_biotype_class = PC_CDS trivially (alt IS canonical, same transcript),
    and always protein_change_type = "identical" downstream in J2 (canonical
    vs itself), so it can never be a gate driver regardless of its usage.
    Applying the 5% usage filter to it would only ever discard informational
    rows for no gate-outcome benefit -- worse, if canonical's own AD usage
    happens to have collapsed below 5% (its usage did drop, that's the whole
    premise of trial_failure), the filter would delete the one row showing
    that collapse, and could leave a hit with zero rows at all even though
    canonical clearly still exists in the data.

Only trial_failure's non-canonical ranked alternates do real classification
+ filtering work here, reusing assistant_surveyor's exact biotype_class logic
(same annotation TSV, same MANE-based CDS diff, same _BIOTYPE_CLASS_MAP),
applied per-row to alt_transcript_name instead of the hit's own transcript_name,
using the canonical each row already resolved via MANE (canonical_enst) -- no
need to re-resolve canonical.

run() adds alt_biotype_class and returns the row-filtered df:
  - new_target: passed through unchanged
  - trial_failure is_canonical rows: passed through unchanged (alt_biotype_class
    set to PC_CDS directly)
  - trial_failure's other ranked alternates: keeps PC_CDS / PC_UTR / PC_CDS_ND
    / novel / TEC, drops RI / NMD / other; then drops rows below MIN_AD_USAGE_PCT
  - warns (does not drop) if any PC_CDS/PC_UTR row has no resolved
    alt_protein_seq -- GENCODE says these biotypes should always have a
    defined CDS, so a miss there is a real gap worth a manual look, unlike
    PC_CDS_ND/novel/TEC where an unresolved sequence is expected.

A trial_failure hit can end up with only its is_canonical row surviving (every
other ranked isoform got dropped by the biotype or usage filter) -- that hit
can never pass J4's gate (the is_canonical row is always protein-identical to
itself, never a driver), but it's kept rather than dropped so the pipeline
stays queryable about *why* a hit failed, not just *whether* it did. run() adds
has_viable_alt (True/False, broadcast to every row of the hit) so this can be
filtered on directly instead of re-deriving it with a groupby. Always True for
new_target (single alt_rank=0 row, no separate "viable alt" question there).
"""

from __future__ import annotations

import sys

import pandas as pd

from assistant_surveyor.config import ANNOT_TSV, JUNIOR_PASS_BIOTYPES
from assistant_surveyor.l1_biotype import _ANN_COLS, _BIOTYPE_CLASS_MAP, _cds_diff_bp, _cds_intervals
from junior_surveyor.config import MIN_AD_USAGE_PCT, TX_ID_MAP


def _log(msg: str) -> None:
    print(f"  [j1c] {msg}", file=sys.stderr, flush=True)


def run(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    nt_mask       = df["candidate_group"] == "new_target_candidate"
    is_canon_bool = df["is_canonical"].fillna(False).astype(bool)
    canon_mask    = (~nt_mask) & is_canon_bool
    alt_mask      = (~nt_mask) & (~is_canon_bool)

    nt    = df[nt_mask].copy()
    canon = df[canon_mask].copy()
    tf    = df[alt_mask].copy()

    nt["alt_biotype_class"] = nt["biotype_class"]
    _log(f"new_target: {len(nt):,} rows passed through unchanged "
         f"(alt_biotype_class copied from biotype_class, no filter applies)")

    canon["alt_biotype_class"] = "PC_CDS"
    _log(f"trial_failure is_canonical rows: {len(canon):,} passed through "
         f"unchanged (alt_biotype_class = PC_CDS trivially, no filter applies)")

    if tf.empty:
        out = pd.concat([nt, canon], ignore_index=True, sort=False)
        _log(f"{len(out):,} rows remain ({len(nt):,} new_target + "
             f"{len(canon):,} trial_failure is_canonical + 0 trial_failure alt)")
        return out

    gene_set = set(tf["gene_name"].unique())
    ann = pd.read_csv(ANNOT_TSV, sep="\t", usecols=_ANN_COLS, low_memory=False)
    ann = ann[
        (ann["feature"].isin(["transcript", "CDS", "UTR"])) & (ann["gene_name"].isin(gene_set))
    ].copy()
    tx_type = ann[ann["feature"] == "transcript"].set_index("transcript_name")["transcript_type"].to_dict()

    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID"]]
    enst_to_name = dict(zip(
        id_map["ENST_ID"].astype(str).str.split(".").str[0], id_map["transcript_name"]
    ))

    def _classify(row: pd.Series) -> str:
        tx_type_ = tx_type.get(row["alt_transcript_name"], "unknown")
        if tx_type_ in _BIOTYPE_CLASS_MAP:
            return _BIOTYPE_CLASS_MAP[tx_type_]
        if tx_type_ == "protein_coding":
            alt_cds = _cds_intervals(ann, row["alt_transcript_name"])
            canon_name = enst_to_name.get(str(row["canonical_enst"]))
            canon_cds = _cds_intervals(ann, canon_name) if canon_name else []
            return "PC_CDS" if _cds_diff_bp(alt_cds, canon_cds) > 0 else "PC_UTR"
        if pd.isna(tx_type_):
            return "novel"
        return "other"

    _log(f"classifying alt_biotype_class for {len(tf):,} trial_failure "
         f"(non-canonical) rows …")
    tf["alt_biotype_class"] = tf.apply(_classify, axis=1)
    _log("alt_biotype_class distribution (trial_failure, non-canonical):")
    for cls, cnt in tf["alt_biotype_class"].value_counts().items():
        _log(f"  {cls}: {cnt}")

    keep = tf["alt_biotype_class"].isin(JUNIOR_PASS_BIOTYPES)
    _log(f"dropping {(~keep).sum():,} trial_failure rows with alt_biotype_class in "
         f"{{RI, NMD, other}} (kept: {keep.sum():,})")
    tf = tf[keep].copy()

    under_min = tf["alt_usage_pct_AD"] < MIN_AD_USAGE_PCT
    _log(f"dropping {under_min.sum():,} trial_failure rows with "
         f"alt_usage_pct_AD < {MIN_AD_USAGE_PCT:.0%}")
    tf = tf[~under_min].copy()

    out = pd.concat([nt, canon, tf], ignore_index=True, sort=False)

    out["has_viable_alt"] = True
    tf_out_mask = out["candidate_group"] == "trial_failure_candidate"
    hits_with_alt = set(map(tuple, tf[["gene_name", "cell_type"]].drop_duplicates().values))
    out.loc[tf_out_mask, "has_viable_alt"] = [
        (g, c) in hits_with_alt
        for g, c in zip(out.loc[tf_out_mask, "gene_name"], out.loc[tf_out_mask, "cell_type"])
    ]
    n_no_alt = (~out.loc[tf_out_mask, "has_viable_alt"]).sum()
    n_hits_no_alt = (
        out.loc[tf_out_mask & ~out["has_viable_alt"], ["gene_name", "cell_type"]]
        .drop_duplicates().shape[0]
    )
    _log(f"has_viable_alt = False: {n_no_alt:,} rows ({n_hits_no_alt:,} trial_failure "
         f"hits with no surviving non-canonical alt -- can never pass J4)")

    pc_rows = out[out["alt_biotype_class"].isin({"PC_CDS", "PC_UTR"})]
    no_seq  = pc_rows[pc_rows["alt_protein_seq"].fillna("") == ""]
    if len(no_seq):
        _log(f"WARNING: {len(no_seq)} PC_CDS/PC_UTR rows have no resolved "
             f"alt_protein_seq -- unexpected, needs manual review:")
        for _, r in no_seq.iterrows():
            _log(f"    {r['gene_name']} / {r['cell_type']} / {r['alt_transcript_name']}")

    _log(f"{len(out):,} rows remain ({len(nt):,} new_target + "
         f"{len(canon):,} trial_failure is_canonical + {len(tf):,} trial_failure alt)")
    return out

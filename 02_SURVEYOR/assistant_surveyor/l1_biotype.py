"""L1 — Biotype & CDS Classification.

Reads the annotation TSV once (filtered to transcript/CDS/UTR rows for genes
in the hit set), determines the canonical transcript per gene, computes CDS
interval diffs between hit and canonical, and assigns biotype_class.

No external I/O — all annotation is local.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from assistant_surveyor.config import ANNOT_TSV, HITS_CSV, JUNIOR_PASS_BIOTYPES, TX_ID_MAP

# Columns we need from the 30-column annotation TSV.
_ANN_COLS = [
    "feature", "start", "end", "strand",
    "gene_name", "transcript_name", "transcript_type",
    "tag", "ccdsid", "transcript_support_level",
]

_BIOTYPE_CLASS_MAP: dict[str, str] = {
    "retained_intron":                "RI",
    "nonsense_mediated_decay":        "NMD",
    "protein_coding_CDS_not_defined": "PC_CDS_ND",
    "TEC":                            "TEC",
}


# ---------------------------------------------------------------------------
# Canonical transcript selection (same hierarchy as plot_tx_structure.py)
# ---------------------------------------------------------------------------
"""
Takes all transcript rows for one gene and returns the name of the most authoritative isoform.
Three-level fallback:

1. If any transcript has Ensembl_canonical in its tag field → that's the one. Done.
2. Otherwise, filter to transcripts with basic tag and a CCDS ID (i.e., confirmed coding sequence),
then pick the one with the lowest Transcript Support Level number (TSL 1 = most evidence, TSL 5 = least).
3. If no CCDS, fall back to any basic-tagged transcript with the lowest TSL.
"""

def _get_canonical(tx_rows: pd.DataFrame) -> str | None:
    tag = tx_rows["tag"].fillna("")
    has_canonical = tag.str.contains("Ensembl_canonical", regex=False)
    if has_canonical.any():
        return tx_rows[has_canonical].iloc[0]["transcript_name"]

    has_basic = tag.str.contains("basic", regex=False)
    has_ccds  = tx_rows["ccdsid"].notna()

    cands = tx_rows[has_basic & has_ccds].copy()
    if not cands.empty:
        cands = cands.sort_values("transcript_support_level", na_position="last")
        return cands.iloc[0]["transcript_name"]

    basic_only = tx_rows[has_basic].copy()
    if not basic_only.empty:
        basic_only = basic_only.sort_values("transcript_support_level", na_position="last")
        return basic_only.iloc[0]["transcript_name"]

    return None


# ---------------------------------------------------------------------------
# CDS interval helpers
# ---------------------------------------------------------------------------

# Pulls all CDS rows for a given transcript from the annotation and returns them
# as a list of (start, end) coordinate tuples. A transcript with no CDS rows
# (e.g. a retained intron or UTR-only transcript) returns an empty list.

def _cds_intervals(ann: pd.DataFrame, tx_name: str) -> list[tuple[int, int]]:
    rows = ann[(ann["transcript_name"] == tx_name) & (ann["feature"] == "CDS")]
    return [(int(r["start"]), int(r["end"])) for _, r in rows.iterrows()]


def _interval_bp(intervals: list[tuple[int, int]]) -> set[int]:
    """Expand a list of (start, end) intervals into a set of bp positions."""
    positions: set[int] = set()
    for s, e in intervals:
        positions.update(range(s, e + 1))
    return positions


def _cds_diff_bp(hit_intervals: list[tuple[int, int]],
                 ref_intervals: list[tuple[int, int]]) -> int:
    """Total bp in the symmetric difference of two CDS interval sets."""
    hit_bp = _interval_bp(hit_intervals)
    ref_bp = _interval_bp(ref_intervals)
    return len(hit_bp.symmetric_difference(ref_bp))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(hits: pd.DataFrame) -> pd.DataFrame:
    """Enrich *hits* with L1 columns. Returns a new DataFrame."""
    print("[L1] Loading tx_id_map ...", flush=True)
    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID", "ENSG_ID"]]

    hit_genes = set(hits["gene_name"].unique())
    hit_txs   = set(hits["transcript_name"].unique())
    print(f"[L1] {len(hit_genes)} unique genes, {len(hit_txs)} unique transcripts",
          flush=True)

    print("[L1] Loading annotation TSV (transcript/CDS/UTR rows) ...", flush=True)
    ann = pd.read_csv(
        ANNOT_TSV, sep="\t", usecols=_ANN_COLS, low_memory=False,
    )
    ann = ann[
        (ann["feature"].isin(["transcript", "CDS", "UTR"]))
        & (ann["gene_name"].isin(hit_genes))
    ].copy()
    print(f"[L1] Annotation filtered to {len(ann):,} rows "
          f"({ann['feature'].value_counts().to_dict()})", flush=True)

    tx_rows = ann[ann["feature"] == "transcript"]

    # -- canonical per gene --------------------------------------------------
    print("[L1] Identifying canonical transcripts per gene ...", flush=True)
    canonical: dict[str, str | None] = {}
    for gene, grp in tx_rows.groupby("gene_name"):
        canonical[gene] = _get_canonical(grp)

    # -- transcript_type per tx ---------------------------------------------
    tx_type: dict[str, str] = (
        tx_rows.set_index("transcript_name")["transcript_type"]
        .to_dict()
    )

    # -- ct_count per transcript (precomputed for O(1) lookup in the loop) --
    ct_count_map: dict[str, int] = (
        hits.groupby("transcript_name")["cell_type"].nunique().to_dict()
    )

    # -- per-hit enrichment --------------------------------------------------
    print("[L1] Computing CDS diffs and assigning biotype_class ...", flush=True)

    rows_out = []
    for _, hit in hits.iterrows():
        tx_name  = hit["transcript_name"]
        gene     = hit["gene_name"]
        canon    = canonical.get(gene)
        tx_type_ = tx_type.get(tx_name, "unknown")

        # CDS intervals for this transcript and its canonical counterpart
        hit_cds   = _cds_intervals(ann, tx_name)
        canon_cds = _cds_intervals(ann, canon) if canon else []
        has_cds   = len(hit_cds) > 0
        diff_bp   = _cds_diff_bp(hit_cds, canon_cds)

        # biotype_class assignment
        # NaN transcript_type = no GENCODE curation at all -> IsoQuant-sourced
        # novel transcript. Any other unmapped GENCODE biotype (lncRNA,
        # non_stop_decay, IG/TR genes, ...) is a real reference transcript
        # and must not be conflated with "novel" -> falls to "other".
        if tx_type_ in _BIOTYPE_CLASS_MAP:
            bc = _BIOTYPE_CLASS_MAP[tx_type_]
        elif tx_type_ == "protein_coding":
            bc = "PC_CDS" if diff_bp > 0 else "PC_UTR"
        elif pd.isna(tx_type_):
            bc = "novel"
        else:
            bc = "other"

        n_ct = ct_count_map.get(tx_name, 1)
        rows_out.append({
            "transcript_type": tx_type_,
            "biotype_class":   bc,
            "junior_pass":     bc in JUNIOR_PASS_BIOTYPES,
            "has_CDS":         has_cds,
            "canonical_tx":    canon,
            "cds_diff_bp":     diff_bp,
            "ct_count":        n_ct,
            "multi_ct":        n_ct > 1,
        })

    l1 = pd.DataFrame(rows_out, index=hits.index)

    # merge ENSG_ID / ENST_ID
    hits_with_ids = hits.merge(
        id_map, on="transcript_name", how="left"
    )

    result = pd.concat([hits_with_ids.reset_index(drop=True),
                        l1.reset_index(drop=True)], axis=1)

    # summary
    print("[L1] biotype_class distribution:", flush=True)
    for cls, cnt in result["biotype_class"].value_counts().items():
        print(f"       {cls}: {cnt}", flush=True)
    print(f"[L1] junior_pass: {result['junior_pass'].sum()} / {len(result)}", flush=True)
    print(f"[L1] multi_ct transcripts: {result['multi_ct'].sum()}", flush=True)

    return result

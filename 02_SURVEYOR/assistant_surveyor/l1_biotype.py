"""L1 — Biotype & CDS Classification.

Reads the annotation TSV once (filtered to transcript/CDS/UTR rows for genes
in the hit set), resolves each gene's canonical transcript via MANE Select
(the same source used everywhere else in the pipeline -- classify_hit_scenarios_mane.py,
junior_surveyor's j1_canonical.py), computes CDS interval diffs between hit
and canonical, and assigns biotype_class.

Every hit reaching this file already has MANE coverage guaranteed --
initial_filter.py only keeps trial_failure_candidate / new_target_candidate,
both of which require tx_role_mane in {Canonical, Alternate}, never
no_MANE_coverage -- so the MANE lookup below never needs a fallback.

No external I/O — all annotation is local.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from classify_hit_scenarios_mane import load_mane_select

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
    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID", "ENSG_ID", "gene_name"]]

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

    # -- canonical per gene, via MANE Select ---------------------------------
    print("[L1] Resolving canonical transcripts via MANE Select ...", flush=True)
    mane = load_mane_select()
    ensg_to_mane_enst = dict(zip(mane["ENSG_ID"], mane["mane_ENST_ID"]))
    enst_to_name = dict(zip(
        id_map["ENST_ID"].astype(str).str.split(".").str[0], id_map["transcript_name"]
    ))
    gene_to_ensg = dict(zip(id_map["gene_name"], id_map["ENSG_ID"]))

    canonical: dict[str, str | None] = {}
    for gene in hit_genes:
        mane_enst = ensg_to_mane_enst.get(gene_to_ensg.get(gene))
        canonical[gene] = enst_to_name.get(mane_enst)

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
        tx_type_ = tx_type.get(tx_name, "unknown")

        # tx_role_mane == "Canonical" means this hit transcript IS the MANE
        # canonical itself (trial_failure_candidate hits are Canonical x
        # CT_enriched by construction) -- diffing it against canonical would
        # just trivially return 0 (comparing a transcript to itself), so skip
        # the interval diff and assign directly. PC_CDS is used here (rather
        # than falling through to the general diff_bp==0 -> PC_UTR rule)
        # since "same CDS as itself" isn't a meaningful UTR-only-variant
        # claim -- it's the reference transcript, not a variant of it.
        if hit.get("tx_role_mane") == "Canonical":
            canon     = tx_name
            has_cds   = True
            diff_bp   = 0
            bc        = "PC_CDS"
        else:
            canon = canonical.get(gene)

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
        id_map, on=["transcript_name", "gene_name"], how="left"
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

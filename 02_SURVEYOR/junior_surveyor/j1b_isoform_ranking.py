"""J1B — ranked alternate isoforms for the trial_failure branch.

For each trial_failure_candidate hit (a gene x cell_type pair where the MANE
canonical isoform is dominant but losing usage in AD), ranks EVERY isoform of
that gene present in the primary pseudobulk by mean AD-donor PSI (no usage
cutoff -- all isoforms are kept; downstream filtering is a manual/gate decision,
not a ranking-stage one). Each ranked isoform becomes a row in the long-format
output, sequence-fetched and later diffed against the gene's canonical protein
by J2 -- this is the "list of alternate protein sequences expressed for that
cell type for that transcript, ranked by usage percent" the trial_failure
analysis needs.

PSI per donor = count(isoform) / sum(counts of all the gene's isoforms), from
outputs/00_PreAggregation_QC/pseudobulk/counts_{cell_type}.csv -- the same primary
pseudobulk 01_ViscachaDTU_Analysis's satuRn DTU test was run on (mirrors the _donor_psi pattern in
layer2/m07_expression.py, generalized here to rank all isoforms of a gene
instead of profiling one pre-selected transcript).
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from junior_surveyor.config import PSEUDOBULK_DIR, TX_ID_MAP
from junior_surveyor import j1_canonical, j1_sequences

_AD_LABEL = "AD"
_CT_LABEL = "Control"


def _log(msg: str) -> None:
    print(f"  [j1b] {msg}", file=sys.stderr, flush=True)


def _load_pb(cell_type: str) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    cpath = PSEUDOBULK_DIR / f"counts_{cell_type}.csv"
    mpath = PSEUDOBULK_DIR / f"metadata_{cell_type}.csv"
    if not cpath.exists() or not mpath.exists():
        return None, None
    counts = pd.read_csv(cpath, index_col=0)
    meta = pd.read_csv(mpath, index_col=0)
    return counts, meta


def _rank_isoforms(counts: pd.DataFrame, meta: pd.DataFrame,
                    gene_cols: list[str]) -> list[dict]:
    """Rank every isoform column of a gene by mean AD-donor PSI, descending."""
    gene_total = counts[gene_cols].sum(axis=1)
    ad_donors = [d for d in meta.index[meta["condition"] == _AD_LABEL] if d in counts.index]
    ct_donors = [d for d in meta.index[meta["condition"] == _CT_LABEL] if d in counts.index]

    results = []
    for col in gene_cols:
        psi = counts[col] / gene_total.replace(0, np.nan)
        ad_vals = psi.loc[ad_donors].dropna()
        ct_vals = psi.loc[ct_donors].dropna()
        results.append({
            "transcript_name":   col,
            "usage_pct_AD":      float(ad_vals.mean()) if len(ad_vals) else float("nan"),
            "usage_pct_control": float(ct_vals.mean()) if len(ct_vals) else float("nan"),
            "n_donors_AD":       len(ad_vals),
        })

    results.sort(key=lambda r: (-1.0 if r["usage_pct_AD"] == r["usage_pct_AD"] else 1.0,
                                 -(r["usage_pct_AD"] if r["usage_pct_AD"] == r["usage_pct_AD"] else 0.0)))
    return results


def run(trial_failure_hits: pd.DataFrame) -> pd.DataFrame:
    _log(f"{len(trial_failure_hits):,} trial_failure hits  "
         f"({trial_failure_hits['gene_name'].nunique():,} unique genes)")

    enst_to_cds, ensg_to_enst_list = j1_sequences.load_fasta()

    gene_df = trial_failure_hits[["ENSG_ID", "canonical_tx"]].drop_duplicates("ENSG_ID")
    canon_map = j1_canonical.resolve_canonical(gene_df, ensg_to_enst_list, enst_to_cds)

    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID"]]
    name_to_enst = dict(zip(id_map["transcript_name"], id_map["ENST_ID"]))

    pb_cache: dict[str, tuple] = {}
    all_rows: list[dict] = []
    n_no_pb, n_no_gene_cols = 0, 0

    for i, (_, hit) in enumerate(trial_failure_hits.iterrows()):
        gene, ct, ensg = hit["gene_name"], hit["cell_type"], hit["ENSG_ID"]

        if ct not in pb_cache:
            pb_cache[ct] = _load_pb(ct)
        counts, meta = pb_cache[ct]
        if counts is None:
            n_no_pb += 1
            continue

        gene_cols = [c for c in counts.columns if c.startswith(gene + "-")]
        if not gene_cols:
            n_no_gene_cols += 1
            continue

        ranked = _rank_isoforms(counts, meta, gene_cols)

        canon_enst, canon_source = canon_map.get(ensg, ("", ""))
        canon_cds     = enst_to_cds.get(canon_enst, "")
        canon_protein = j1_sequences.translate(canon_cds)

        for rank, iso in enumerate(ranked, start=1):
            row = hit.to_dict()
            alt_enst = str(name_to_enst.get(iso["transcript_name"], "")).split(".")[0]
            ad_pct, ct_pct = iso["usage_pct_AD"], iso["usage_pct_control"]
            usage_delta = (ad_pct - ct_pct
                            if ad_pct == ad_pct and ct_pct == ct_pct  # both non-NaN
                            else float("nan"))
            row.update({
                "alt_rank":              rank,
                "hit_ENST_ID":           hit["ENST_ID"],
                "hit_transcript_name":   hit["transcript_name"],
                "alt_transcript_name":   iso["transcript_name"],
                "alt_ENST_ID":           alt_enst,
                "ENST_ID":               alt_enst,   # downstream J2/J3 compat
                "alt_usage_pct_AD":      ad_pct,
                "alt_usage_pct_control": ct_pct,
                "alt_usage_delta":       usage_delta,
                "is_canonical":          bool(canon_enst) and alt_enst == canon_enst,
                "n_donors_AD":           iso["n_donors_AD"],
                "canonical_enst":        canon_enst,
                "canonical_source":      canon_source,
                "canonical_cds_seq":     canon_cds,
                "canonical_protein_seq": canon_protein,
            })
            all_rows.append(row)

        if (i + 1) % 50 == 0:
            _log(f"  ranked {i+1:,}/{len(trial_failure_hits):,} hits "
                 f"({len(all_rows):,} alt rows so far)")

    if n_no_pb:
        _log(f"  {n_no_pb:,} hits skipped — no pseudobulk for cell type")
    if n_no_gene_cols:
        _log(f"  {n_no_gene_cols:,} hits skipped — gene not in pseudobulk")

    out = pd.DataFrame(all_rows)
    _log(f"ranked-isoform rows: {len(out):,}  "
         f"({out['gene_name'].nunique():,} genes, "
         f"avg {len(out) / max(trial_failure_hits['gene_name'].nunique(), 1):.1f} rows/gene)")

    _log("fetching + translating ranked alternate sequences …")
    out = j1_sequences.fetch_alt_sequences(out, enst_to_cds)

    return out

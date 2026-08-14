"""J1_CANONICAL — unified canonical-transcript resolver, shared by both branches.

Priority per gene (ENSG_ID):
  1. MANE Select (authoritative; loaded via classify_hit_scenarios_mane.load_mane_select)
  2. tx_id_map's canonical_tx name -> ENST_ID (Ensembl_canonical tag / TSL heuristic,
     computed upstream by assistant_surveyor's L1 stage)
  3. Longest CDS among the gene's transcripts in the Ensembl FASTA

Genes with no MANE Select entry (non-coding biotypes) fall through to 2/3, matching
the "no_MANE_coverage" concept already used upstream in initial_filter.py.
"""

from __future__ import annotations

import sys

import pandas as pd

from classify_hit_scenarios_mane import load_mane_select
from junior_surveyor.config import TX_ID_MAP


def _log(msg: str) -> None:
    print(f"  [j1_canonical] {msg}", file=sys.stderr, flush=True)


def resolve_canonical(
    gene_df: pd.DataFrame,
    ensg_to_enst_list: dict[str, list[str]],
    enst_to_cds: dict[str, str],
) -> dict[str, tuple[str, str]]:
    """Return {ENSG_ID: (canonical_enst_noversion, source)}.

    source is one of "MANE", "heuristic_tx_id_map", "heuristic_longest_cds".
    gene_df must have columns ENSG_ID, canonical_tx (one row per ENSG_ID).
    """
    mane = load_mane_select()
    mane_map = dict(zip(mane["ENSG_ID"], mane["mane_ENST_ID"]))

    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID"]]
    name_to_enst = dict(zip(id_map["transcript_name"], id_map["ENST_ID"]))

    result: dict[str, tuple[str, str]] = {}
    source_counts = {"MANE": 0, "heuristic_tx_id_map": 0, "heuristic_longest_cds": 0}

    for _, row in gene_df.iterrows():
        ensg = row["ENSG_ID"]

        mane_enst = mane_map.get(ensg)
        if mane_enst and mane_enst in enst_to_cds:
            result[ensg] = (mane_enst, "MANE")
            source_counts["MANE"] += 1
            continue

        canon_name = row.get("canonical_tx", "")
        if canon_name and canon_name in name_to_enst:
            enst = str(name_to_enst[canon_name]).split(".")[0]
            if enst in enst_to_cds:
                result[ensg] = (enst, "heuristic_tx_id_map")
                source_counts["heuristic_tx_id_map"] += 1
                continue

        candidates = ensg_to_enst_list.get(ensg, [])
        if candidates:
            best = max(candidates, key=lambda e: len(enst_to_cds.get(e, "")))
            if enst_to_cds.get(best):
                result[ensg] = (best, "heuristic_longest_cds")
                source_counts["heuristic_longest_cds"] += 1

    _log(f"canonical resolved: {len(result):,}/{len(gene_df):,} genes  "
         f"(MANE={source_counts['MANE']:,}, "
         f"tx_id_map={source_counts['heuristic_tx_id_map']:,}, "
         f"longest_cds={source_counts['heuristic_longest_cds']:,})")

    return result

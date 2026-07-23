"""J1 — CDS + protein sequences from local Ensembl FASTA.

Parses Homo_sapiens.GRCh38.cds.all.fa.gz (23MB, 123k sequences) into a
local dict keyed by ENST_ID.  No REST API calls needed.

FASTA header format:
  >ENST00000641515.2 cds chromosome:GRCh38:... gene:ENSG00000186092.7 ...

Canonical transcript per gene: resolved via local tx_id_map first (covers
~80%), then by selecting the longest CDS among the gene's transcripts in
the FASTA for the remainder.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from junior_surveyor.config import CACHE_DIR, JUNIOR_PASS_BIOTYPES
from junior_surveyor import j0_novel_orfs

_FASTA_GZ = Path(__file__).resolve().parent.parent / \
    "outputs/junior_surveyor/ref/Homo_sapiens.GRCh38.cds.all.fa.gz"

_TX_ID_MAP = Path(__file__).resolve().parent.parent / \
    "outputs/layer1_1/annotation/tx_id_map.csv"


def _log(msg: str) -> None:
    print(f"  [j1] {msg}", file=sys.stderr, flush=True)


def _load_fasta() -> tuple[dict[str, str], dict[str, list[str]]]:
    """Parse CDS FASTA → (enst_to_cds, ensg_to_enst_list).

    enst_to_cds      : {ENST_ID_noversion: cds_sequence}
    ensg_to_enst_list: {ENSG_ID_noversion: [ENST_ID, ...]}
    """
    cache = CACHE_DIR / "j1_fasta_index.json"
    if cache.exists():
        _log("FASTA index: loading from cache")
        data = json.loads(cache.read_text())
        return data["enst_to_cds"], data["ensg_to_enst_list"]

    _log(f"parsing {_FASTA_GZ.name} …")
    enst_to_cds:       dict[str, str]        = {}
    ensg_to_enst_list: dict[str, list[str]]  = {}

    cur_enst = cur_ensg = ""
    seq_parts: list[str] = []

    def _commit() -> None:
        if cur_enst and seq_parts:
            enst_to_cds[cur_enst] = "".join(seq_parts)
            ensg_to_enst_list.setdefault(cur_ensg, []).append(cur_enst)

    with gzip.open(_FASTA_GZ, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                _commit()
                seq_parts = []
                parts = line[1:].split()
                # Strip version: ENST00000641515.2 → ENST00000641515
                cur_enst = parts[0].split(".")[0]
                cur_ensg = ""
                for p in parts:
                    if p.startswith("gene:ENSG"):
                        cur_ensg = p.split(":")[1].split(".")[0]
                        break
            else:
                seq_parts.append(line)
    _commit()

    _log(f"  {len(enst_to_cds):,} ENST sequences, {len(ensg_to_enst_list):,} genes")
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"enst_to_cds": enst_to_cds,
                                  "ensg_to_enst_list": ensg_to_enst_list}))
    return enst_to_cds, ensg_to_enst_list


def _resolve_canonical(
    gene_df: pd.DataFrame,
    ensg_to_enst_list: dict[str, list[str]],
    enst_to_cds: dict[str, str],
) -> dict[str, str]:
    """Return {ensg_id: canonical_enst_id}.

    Priority:
      1. tx_id_map canonical_tx name → ENST_ID (local, fast)
      2. Longest CDS among the gene's transcripts in the FASTA
    """
    id_map = pd.read_csv(_TX_ID_MAP)[["transcript_name", "ENST_ID"]]
    name_to_enst = dict(zip(id_map["transcript_name"], id_map["ENST_ID"]))

    result: dict[str, str] = {}
    for _, row in gene_df.iterrows():
        ensg       = row["ENSG_ID"]
        canon_name = row.get("canonical_tx", "")

        # Priority 1: tx_id_map lookup
        if canon_name and canon_name in name_to_enst:
            enst = name_to_enst[canon_name].split(".")[0]
            if enst in enst_to_cds:
                result[ensg] = enst
                continue

        # Priority 2: longest CDS in FASTA
        candidates = ensg_to_enst_list.get(ensg, [])
        if candidates:
            best = max(candidates, key=lambda e: len(enst_to_cds.get(e, "")))
            if enst_to_cds.get(best):
                result[ensg] = best

    return result


def _translate(cds: str) -> str:
    if not cds or len(cds) < 3:
        return ""
    cds = cds[: (len(cds) // 3) * 3]
    return str(Seq(cds).translate(to_stop=False)).rstrip("*")


def run(df: pd.DataFrame) -> pd.DataFrame:
    sub = df[df["biotype_class"].isin(JUNIOR_PASS_BIOTYPES)].copy()
    _log(f"{len(sub):,} protein-affecting hits  "
         f"({sub['ENST_ID'].nunique():,} unique DIU ENST, "
         f"{sub['ENSG_ID'].nunique():,} unique genes)")

    # ── Load FASTA index ─────────────────────────────────────────────────────
    enst_to_cds, ensg_to_enst_list = _load_fasta()

    # ── Alt CDS (by ENST_ID) — Ensembl FASTA first ──────────────────────────
    # Strip version suffix from real Ensembl IDs (e.g. ENST00000641515.2 → .515)
    # Novel IDs (transcript39631.chr6.nic) don't have version suffixes so this
    # is harmless — they simply won't match anything in the Ensembl dict.
    sub["_enst_novers"] = sub["ENST_ID"].str.split(".").str[0]
    sub["alt_cds_seq"]  = sub["_enst_novers"].map(enst_to_cds).fillna("")
    sub.drop(columns=["_enst_novers"], inplace=True)

    n_ensembl = (sub["alt_cds_seq"] != "").sum()
    _log(f"alt CDS from Ensembl FASTA: {n_ensembl:,}/{len(sub):,}")

    # ── Alt CDS fallback — TransDecoder ORF prediction for novel transcripts ─
    missing_mask = sub["alt_cds_seq"] == ""
    novel_ids    = set(sub.loc[missing_mask, "ENST_ID"].tolist())
    novel_orfs   = j0_novel_orfs.run(novel_ids) if novel_ids else {}

    if novel_orfs:
        # Fill gaps where TransDecoder found an ORF
        fill_mask = missing_mask & sub["ENST_ID"].isin(novel_orfs)
        sub.loc[fill_mask, "alt_cds_seq"] = (
            sub.loc[fill_mask, "ENST_ID"].map(novel_orfs)
        )
        n_td = fill_mask.sum()
        _log(f"alt CDS from TransDecoder: {n_td:,}/{missing_mask.sum():,} novel transcripts")

    # ── Source flag ─────────────────────────────────────────────────────────
    # Records where each alt CDS came from so the final CSV is traceable.
    def _cds_source(row: pd.Series) -> str:
        if row["alt_cds_seq"] == "":
            return "none"
        if row["ENST_ID"] in novel_orfs:
            return "transdecoder"
        return "ensembl"

    sub["alt_cds_source"] = sub.apply(_cds_source, axis=1)
    _log(f"alt_cds_source counts: {sub['alt_cds_source'].value_counts().to_dict()}")

    # ── Canonical CDS (per gene, always from Ensembl FASTA) ─────────────────
    gene_df = sub[["ENSG_ID", "canonical_tx"]].drop_duplicates("ENSG_ID")
    ensg_to_canon = _resolve_canonical(gene_df, ensg_to_enst_list, enst_to_cds)
    _log(f"canonical ENST resolved: {len(ensg_to_canon):,}/{len(gene_df):,} genes")

    sub["canonical_enst"]    = sub["ENSG_ID"].map(ensg_to_canon)
    sub["canonical_cds_seq"] = sub["canonical_enst"].map(enst_to_cds).fillna("")

    n_can = (sub["canonical_cds_seq"] != "").sum()
    _log(f"canonical CDS matched: {n_can:,}/{len(sub):,}")

    # ── Translate ────────────────────────────────────────────────────────────
    _log("translating sequences …")
    sub["alt_protein_seq"]       = sub["alt_cds_seq"].apply(_translate)
    sub["canonical_protein_seq"] = sub["canonical_cds_seq"].apply(_translate)

    _log(f"alt protein OK:       {(sub['alt_protein_seq'] != '').sum():,}/{len(sub):,}")
    _log(f"canonical protein OK: {(sub['canonical_protein_seq'] != '').sum():,}/{len(sub):,}")

    return sub

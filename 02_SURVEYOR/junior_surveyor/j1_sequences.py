"""J1 — CDS + protein sequences for the new_target branch (alt_rank=0 rows).

Parses Homo_sapiens.GRCh38.cds.all.fa.gz (23MB, 123k sequences) into a
local dict keyed by ENST_ID.  No REST API calls needed.

For each new_target_candidate hit, resolves the gene's canonical transcript via
j1_canonical (MANE-first) and fetches/translates the DTU-tested transcript's own
CDS as the "alt" -- the single pairwise comparison the new_target branch needs
(candidate DTU transcript vs canonical).

FASTA loading (`load_fasta`), translation (`translate`), and alt-sequence
fetching (`fetch_alt_sequences`, incl. TransDecoder fallback for novel
transcripts) are reusable and imported by junior_surveyor/j1b_isoform_ranking.py
for the trial_failure branch, which applies the same logic to many ranked
alternates per gene instead of one.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from junior_surveyor.config import CACHE_DIR
from junior_surveyor import j0_novel_orfs, j0b_sqanti_orfs, j1_canonical

_FASTA_GZ = Path(__file__).resolve().parent.parent.parent / \
    "outputs/junior_surveyor/ref/Homo_sapiens.GRCh38.cds.all.fa.gz"


def _log(msg: str) -> None:
    print(f"  [j1] {msg}", file=sys.stderr, flush=True)


def load_fasta() -> tuple[dict[str, str], dict[str, list[str]]]:
    """Parse CDS FASTA -> (enst_to_cds, ensg_to_enst_list).

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


def translate(cds: str) -> str:
    if not cds or len(cds) < 3:
        return ""
    cds = cds[: (len(cds) // 3) * 3]
    return str(Seq(cds).translate(to_stop=False)).rstrip("*")


def _substring_compatible(a: str, b: str) -> bool:
    """True if the shorter sequence is a substring of the longer (or equal).

    Used to distinguish an interpretable ORF-boundary difference (e.g. one
    caller picked an upstream start codon) from a genuinely conflicting call
    (different reading frame entirely).
    """
    if not a or not b:
        return False
    return a in b or b in a


def fetch_alt_sequences(df: pd.DataFrame, enst_to_cds: dict[str, str]) -> pd.DataFrame:
    """Fill alt_cds_seq / alt_cds_source / alt_protein_seq from df["ENST_ID"].

    Three tiers, in priority order:
      1. Ensembl CDS FASTA (strip version suffix so novel IDs like
         transcript39631.chr6.nic simply miss, which is fine).
      2. SQANTI3's own ORF call (j0b_sqanti_orfs) -- already run on the full
         long-read catalog, resolves ~70% of what tier 3 misses (see
         j0b_sqanti_orfs.py docstring for the reliability comparison behind
         this ordering).
      3. TransDecoder ab-initio ORF prediction (j0_novel_orfs), last resort,
         only for transcripts SQANTI3 also has nothing for.

    Also flags alt_seq_source_conflict = True where a cached TransDecoder
    prediction exists for a SQANTI3-sourced row (from a prior run, or this
    one) and the two calls are NOT substring-compatible -- i.e. a genuine
    reading-frame disagreement worth a manual look, not just a different
    start-codon choice. Doesn't force fresh TransDecoder runs to check this;
    only compares against whatever's already cached.

    Mutates and returns df.
    """
    df["_enst_novers"] = df["ENST_ID"].str.split(".").str[0]
    df["alt_cds_seq"]  = df["_enst_novers"].map(enst_to_cds).fillna("")
    df.drop(columns=["_enst_novers"], inplace=True)

    n_ensembl = (df["alt_cds_seq"] != "").sum()
    _log(f"alt CDS from Ensembl FASTA: {n_ensembl:,}/{len(df):,}")

    missing_mask = df["alt_cds_seq"] == ""

    # ── Tier 2: SQANTI3 ─────────────────────────────────────────────────────
    sqanti_orfs = j0b_sqanti_orfs.load()
    sqanti_hit_mask = missing_mask & df["ENST_ID"].isin(sqanti_orfs)
    n_sqanti = sqanti_hit_mask.sum()
    _log(f"alt protein from SQANTI3: {n_sqanti:,}/{missing_mask.sum():,} missing transcripts")

    # ── Tier 3: TransDecoder, only for what's still missing ─────────────────
    still_missing = missing_mask & ~sqanti_hit_mask
    novel_ids  = set(df.loc[still_missing, "ENST_ID"].tolist())
    novel_orfs = j0_novel_orfs.run(novel_ids) if novel_ids else {}

    if novel_orfs:
        fill_mask = still_missing & df["ENST_ID"].isin(novel_orfs)
        df.loc[fill_mask, "alt_cds_seq"] = df.loc[fill_mask, "ENST_ID"].map(novel_orfs)
        n_td = fill_mask.sum()
        _log(f"alt CDS from TransDecoder: {n_td:,}/{still_missing.sum():,} remaining transcripts")

    source = pd.Series("none", index=df.index)
    source[~missing_mask] = "ensembl"
    source[sqanti_hit_mask] = "sqanti3"
    td_fill_mask = still_missing & df["ENST_ID"].isin(novel_orfs)
    source[td_fill_mask] = "transdecoder"
    df["alt_cds_source"] = source

    # alt_protein_seq: translate for ensembl/transdecoder (CDS nucleotide),
    # use directly for sqanti3 (already a protein).
    df["alt_protein_seq"] = df["alt_cds_seq"].apply(translate)
    sqanti_rows = df["alt_cds_source"] == "sqanti3"
    df.loc[sqanti_rows, "alt_protein_seq"] = df.loc[sqanti_rows, "ENST_ID"].map(sqanti_orfs)

    # ── Conflict flag: SQANTI3 vs. any already-cached TransDecoder call ─────
    td_cache = j0_novel_orfs.load_cache()
    def _conflict(row: pd.Series) -> bool:
        if row["alt_cds_source"] != "sqanti3":
            return False
        td_cds = td_cache.get(row["ENST_ID"])
        if not td_cds:
            return False
        td_protein = translate(td_cds)
        sq_protein = row["alt_protein_seq"]
        if not td_protein or not sq_protein:
            return False
        return not _substring_compatible(sq_protein, td_protein)

    df["alt_seq_source_conflict"] = df.apply(_conflict, axis=1)
    n_conflict = df["alt_seq_source_conflict"].sum()
    if n_conflict:
        _log(f"alt_seq_source_conflict: {n_conflict:,} rows where SQANTI3 and a "
             f"cached TransDecoder call disagree on the reading frame")

    n_alt = (df["alt_protein_seq"].fillna("") != "").sum()
    _log(f"alt protein OK: {n_alt:,}/{len(df):,}")
    return df


def run(new_target_hits: pd.DataFrame) -> pd.DataFrame:
    """New-target branch: one row per hit, alt_rank=0 (DTU transcript vs canonical)."""
    sub = new_target_hits.copy()
    _log(f"{len(sub):,} new_target hits  "
         f"({sub['ENST_ID'].nunique():,} unique DIU ENST, "
         f"{sub['ENSG_ID'].nunique():,} unique genes)")

    enst_to_cds, ensg_to_enst_list = load_fasta()

    sub = fetch_alt_sequences(sub, enst_to_cds)

    # ── Canonical CDS (per gene, MANE-first via j1_canonical) ───────────────
    gene_df = sub[["ENSG_ID", "canonical_tx"]].drop_duplicates("ENSG_ID")
    canon_map = j1_canonical.resolve_canonical(gene_df, ensg_to_enst_list, enst_to_cds)

    sub["canonical_enst"]    = sub["ENSG_ID"].map(lambda g: canon_map.get(g, ("", ""))[0])
    sub["canonical_source"]  = sub["ENSG_ID"].map(lambda g: canon_map.get(g, ("", ""))[1])
    sub["canonical_cds_seq"] = sub["canonical_enst"].map(enst_to_cds).fillna("")

    n_can = (sub["canonical_cds_seq"] != "").sum()
    _log(f"canonical CDS matched: {n_can:,}/{len(sub):,}")

    _log("translating canonical sequences …")
    sub["canonical_protein_seq"] = sub["canonical_cds_seq"].apply(translate)
    _log(f"canonical protein OK: {(sub['canonical_protein_seq'] != '').sum():,}/{len(sub):,}")

    # ── Unified long-format identity columns ─────────────────────────────────
    sub["alt_rank"]              = 0
    sub["hit_ENST_ID"]           = sub["ENST_ID"]
    sub["hit_transcript_name"]   = sub["transcript_name"]
    sub["alt_ENST_ID"]           = sub["ENST_ID"]
    sub["alt_transcript_name"]   = sub["transcript_name"]
    sub["alt_usage_pct_AD"]      = sub["AD"]
    sub["alt_usage_pct_control"] = sub["Control"]

    return sub

"""J0B — SQANTI3 ORF sequences, primary fallback ahead of TransDecoder.

SQANTI3 already classified every isoform in the long-read catalog (known +
novel transcripts alike), including its own ORF/CDS prediction per transcript
(the ORF_seq column of isoforms_classification_with_tx_name_and_gene_name.csv).

Cross-checked against this pipeline's own TransDecoder ab-initio predictions
for the 129 transcripts where both existed: 74% exact protein match, 15%
where SQANTI3's call is a shorter substring of TransDecoder's (TransDecoder's
longest-ORF search is known to favor upstream start codons -- not necessarily
wrong, just a different convention), and 9% genuinely unrelated (different
reading frame -- real ambiguity). SQANTI3 is used as the primary fallback
(ahead of TransDecoder, only after the Ensembl CDS FASTA lookup misses)
because it also independently resolves ~70% of the transcripts TransDecoder's
longest-ORF heuristic missed entirely (mostly PC_CDS_ND biotype, where
Ensembl itself couldn't confidently call a CDS).

Cache: CACHE_DIR / "sqanti_orfs.json" -- {isoform: protein_seq}, built once
from the classification CSV's ORF_seq column (non-null rows only). Keyed by
"isoform" (unversioned ENST accession or novel long-read ID), the same
identifier space as this pipeline's ENST_ID/alt_ENST_ID columns -- NOT the
CSV's "transcript_name" column, which is the gene-symbol-style name used
elsewhere in this pipeline for known transcripts only.
"""

from __future__ import annotations

import json
import sys

import pandas as pd

from junior_surveyor.config import CACHE_DIR, SQANTI_CLASSIFICATION_CSV

_CACHE = CACHE_DIR / "sqanti_orfs.json"


def _log(msg: str) -> None:
    print(f"  [j0b] {msg}", file=sys.stderr, flush=True)


def load() -> dict[str, str]:
    """Return {isoform: protein_seq} for every SQANTI3-called ORF."""
    if _CACHE.exists():
        _log("SQANTI3 ORFs: loading from cache")
        return json.loads(_CACHE.read_text())

    _log(f"parsing {SQANTI_CLASSIFICATION_CSV.name} …")
    df = pd.read_csv(SQANTI_CLASSIFICATION_CSV, usecols=["isoform", "ORF_seq"])
    df = df.dropna(subset=["ORF_seq", "isoform"])
    orfs = {
        iso: seq.rstrip("*")
        for iso, seq in zip(df["isoform"], df["ORF_seq"])
    }
    _log(f"  {len(orfs):,} transcripts with a SQANTI3 ORF")

    _CACHE.parent.mkdir(parents=True, exist_ok=True)
    _CACHE.write_text(json.dumps(orfs))
    return orfs

"""J2 helper — Pfam domain scanning via pyhmmer (in-process HMMER3).

Replaces subprocess hmmscan with pyhmmer, which runs HMMER3 natively in Python:
  - No temporary FASTA or domtblout files needed for the scan itself
  - No subprocess overhead or parsing of text output
  - Equivalent results: same GA thresholds, same coordinate system

Cache strategy (checked in order):
  1. j2_pfam_hits.json  — fast JSON cache of parsed domain hits (written by pyhmmer path)
  2. j2_pfam.domtblout  — legacy cache from previous hmmscan runs; parsed once, then
                          promoted to JSON so future runs skip the domtblout parser

This means the first run after switching to pyhmmer will either:
  - Read the existing domtblout and write a JSON cache (if domtblout exists), or
  - Run a fresh pyhmmer scan and write JSON cache (if neither cache exists).
Subsequent runs always hit the JSON cache.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pyhmmer

from junior_surveyor.config import CACHE_DIR, PFAM_CPU, PFAM_HMM


def _log(msg: str) -> None:
    print(f"  [j2_pfam] {msg}", file=sys.stderr, flush=True)


_CACHE_JSON = CACHE_DIR / "j2_pfam_hits.json"
_CACHE_DOMTBLOUT = CACHE_DIR / "j2_pfam.domtblout"


def _truncate_at_stop(seq: str) -> str:
    idx = seq.find("*")
    return seq[:idx] if idx != -1 else seq


def _collect_sequences(df: pd.DataFrame) -> dict[str, str]:
    """Return {seq_id: sequence} for all unique canonical + alt proteins in df."""
    seqs: dict[str, str] = {}
    for _, row in df.iterrows():
        canon_enst = row.get("canonical_enst", "")
        canon_seq  = row.get("canonical_protein_seq", "")
        if canon_enst and canon_seq:
            seqs.setdefault(f"{canon_enst}__canonical", _truncate_at_stop(canon_seq))
        alt_enst = row.get("ENST_ID", "")
        alt_seq  = row.get("alt_protein_seq", "")
        if alt_enst and alt_seq:
            seqs.setdefault(f"{alt_enst}__alt", _truncate_at_stop(alt_seq))
    return seqs


def _scan_pyhmmer(seqs: dict[str, str]) -> dict[str, list[dict]]:
    """Run pyhmmer hmmscan with GA cutoffs. Returns {seq_id: [{acc,name,start,end,evalue}]}."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital_seqs = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
        pyhmmer.easel.TextSequence(name=seq_id, sequence=seq).digitize(alphabet)
        for seq_id, seq in seqs.items() if seq
    ])

    hits: dict[str, list[dict]] = {}
    with pyhmmer.plan7.HMMFile(str(PFAM_HMM)) as hmm_file:
        for top_hits in pyhmmer.hmmscan(
            digital_seqs, hmm_file,
            cpus=PFAM_CPU,
            bit_cutoffs="gathering",
        ):
            seq_id  = top_hits.query.name
            entries = []
            for hit in top_hits:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    if not domain.included:
                        continue
                    aln = domain.alignment
                    # target_from/to: 1-based positions in the sequence (verified
                    # against existing domtblout — matches ali-from / ali-to columns).
                    entries.append({
                        "acc":    hit.accession,
                        "name":   hit.name,
                        "start":  aln.target_from,
                        "end":    aln.target_to,
                        "evalue": domain.i_evalue,
                    })
            if entries:
                hits[seq_id] = entries

    return hits


def _parse_domtblout(domtab_path: Path) -> dict[str, list[dict]]:
    """Parse legacy hmmscan --domtblout (Biopython SearchIO)."""
    from Bio import SearchIO

    hits: dict[str, list[dict]] = {}
    for query in SearchIO.parse(str(domtab_path), "hmmscan3-domtab"):
        seq_id  = query.id
        entries = []
        for hit in query:
            for hsp in hit:
                entries.append({
                    "acc":    hit.id,
                    "name":   hit.description,
                    "start":  int(hsp.query_start) + 1,   # SearchIO 0-based → 1-based
                    "end":    int(hsp.query_end),
                    "evalue": float(hsp.evalue),
                })
        if entries:
            hits[seq_id] = entries

    n = sum(len(v) for v in hits.values())
    _log(f"parsed {n:,} domain hits across {len(hits):,} sequences")
    return hits


def run(df: pd.DataFrame) -> dict[str, list[dict]]:
    """Full Pfam pipeline. Returns {seq_id: [{acc, name, start, end, evalue}]}."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. JSON cache (fastest path) ──────────────────────────────────────────
    if _CACHE_JSON.exists():
        with _CACHE_JSON.open() as f:
            hits = json.load(f)
        n = sum(len(v) for v in hits.values())
        _log(f"loaded {n:,} domain hits from JSON cache ({len(hits):,} sequences)")
        return hits

    # ── 2. Legacy domtblout cache ──────────────────────────────────────────────
    if _CACHE_DOMTBLOUT.exists():
        _log("promoting legacy domtblout to JSON cache …")
        hits = _parse_domtblout(_CACHE_DOMTBLOUT)
        with _CACHE_JSON.open("w") as f:
            json.dump(hits, f)
        _log(f"JSON cache written → {_CACHE_JSON.name}")
        return hits

    # ── 3. Fresh pyhmmer scan ──────────────────────────────────────────────────
    seqs = _collect_sequences(df)
    _log(f"scanning {len(seqs):,} sequences via pyhmmer "
         f"(GA cutoffs, {PFAM_CPU} CPUs) …")
    hits = _scan_pyhmmer(seqs)

    n = sum(len(v) for v in hits.values())
    _log(f"found {n:,} domain hits across {len(hits):,} sequences")

    with _CACHE_JSON.open("w") as f:
        json.dump(hits, f)
    _log(f"JSON cache written → {_CACHE_JSON.name}")

    return hits

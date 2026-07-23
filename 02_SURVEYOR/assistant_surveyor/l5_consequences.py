"""L5 — Lightweight Consequence Classification.

Assigns proxy_type as an annotation-level approximation of the Surveyor's
Type A/B/C/D/N system. No sequence alignment — uses biotype_class (L1) and
has_structural_feat (L4).

proxy_type values:
  N   — no protein-level difference expected (UTR-only or retained intron)
  NMD — NMD isoform usage shift (regulatory interest, no docking)
  D   — CDS-altering, gene has no known structural domains
  C   — CDS-altering, gene has structural domain(s) that could be affected

Note: frameshift is not determined here. Coordinate-based net CDS length
is an unreliable proxy for true reading frame disruption. Frameshift
detection is deferred to the sequence analysis step.
"""

from __future__ import annotations

import pandas as pd


_PROXY_ORDER = {"C": 0, "D": 1, "frame": 2, "NMD": 3, "N": 4}


def _classify_row(row: pd.Series) -> str:
    bc   = row["biotype_class"]
    feat = bool(row.get("has_structural_feat", False))

    if bc in ("RI", "PC_UTR"):
        return "N"
    if bc == "NMD":
        return "NMD"
    if bc in ("PC_CDS", "PC_CDS_ND", "novel"):
        return "C" if feat else "D"
    return "N"  # unknown/other → conservative


def run(hits: pd.DataFrame) -> pd.DataFrame:
    """Add proxy_type and splice_in_cds columns. Returns new DataFrame."""
    print("[L5] Classifying splice consequences ...", flush=True)

    result = hits.copy()
    result["proxy_type"]   = result.apply(_classify_row, axis=1)
    result["splice_in_cds"] = result["cds_diff_bp"] > 0

    print("[L5] proxy_type distribution:", flush=True)
    for pt, cnt in result["proxy_type"].value_counts().items():
        print(f"       {pt}: {cnt}", flush=True)

    return result

"""JUNIOR_SURVEYOR — main entry point.

Runs j1 → j2 → j3 → j4 on the protein-affecting hits from assistant_surveyor.
Writes outputs/junior_surveyor/hits_deep.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Ensure repo root is on path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from junior_surveyor.config import CACHE_DIR, HITS_CSV, OUT_DIR
from junior_surveyor import j1_sequences, j2_protein_diff, j3_drug_targets, j4_score


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load ────────────────────────────────────────────────────────────────
    print("Loading hits_enriched.csv …")
    all_hits = pd.read_csv(HITS_CSV)
    print(f"  {len(all_hits):,} total hits  |  "
          f"{all_hits['gene_name'].nunique():,} genes\n")

    # ── J1: sequences ───────────────────────────────────────────────────────
    print("J1 — fetching transcript + protein sequences …")
    df = j1_sequences.run(all_hits)
    print(f"  → {len(df):,} protein-affecting hits\n")

    # ── J2: protein diff ────────────────────────────────────────────────────
    print("J2 — aligning canonical vs alternative proteins …")
    df = j2_protein_diff.run(df)
    print()

    # ── J3: drug targets ────────────────────────────────────────────────────
    print("J3 — querying ChEMBL + DGIdb …")
    df = j3_drug_targets.run(df)
    print()

    # ── J4: select hits for next stage ──────────────────────────────────────
    print("J4 — selecting hits for next surveyor stage …")
    df = j4_score.run(df)
    print()

    # ── Save ────────────────────────────────────────────────────────────────
    out_csv = OUT_DIR / "hits_deep.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved → {out_csv}")
    print(f"  {len(df):,} hits  |  "
          f"{df['gene_name'].nunique():,} genes  |  "
          f"{df.shape[1]} columns")


if __name__ == "__main__":
    main()

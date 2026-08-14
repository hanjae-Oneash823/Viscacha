#!/usr/bin/env python3
"""
4-scenario classification of the 2,599 permutation-significant DTU hits using
MANE Select (canonical) vs everything else (alternate) in place of the
empirical Control-PSI dominance rank used in classify_hit_scenarios.py.

  usage_direction (AD_enriched / CT_enriched) x tx_role (Canonical / Alternate)

  Canonical x AD_enriched   -> canonical isoform gaining further usage in AD
  Canonical x CT_enriched   -> canonical isoform losing usage in AD (trial-failure-like)
  Alternate x AD_enriched   -> alternate isoform gaining usage in AD (new-target-like)
  Alternate x CT_enriched   -> alternate isoform losing further usage in AD

MANE Select only covers protein-coding genes and is keyed by ENST, so hits
whose gene has no MANE Select entry (lncRNA/other biotypes) or whose
transcript has no ENST at all (long-read novel calls) get tx_role =
"no_MANE_coverage" instead of being forced into Canonical/Alternate.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python classify_hit_scenarios_mane.py
"""

import gzip
from pathlib import Path

import pandas as pd

REPO_ROOT   = Path(__file__).resolve().parent
HITS_CSV    = REPO_ROOT / "outputs/DIU_significant_hits/DIU_significant_hits_combined.csv"
TX_ID_MAP   = REPO_ROOT / "outputs/annotation/tx_id_map.csv"
MANE_FILE   = REPO_ROOT / "outputs/reference/mane/MANE.GRCh38.v1.5.summary.txt.gz"
OUT_DIR     = REPO_ROOT / "outputs/scenario_analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def strip_version(ensembl_id: str) -> str:
    return ensembl_id.split(".")[0]


def load_mane_select() -> pd.DataFrame:
    """One row per protein-coding gene: gene_name/ENSG -> MANE Select ENST (unversioned)."""
    with gzip.open(MANE_FILE, "rt") as fh:
        mane = pd.read_csv(fh, sep="\t")
    mane = mane[mane["MANE_status"] == "MANE Select"].copy()
    mane["ENSG_ID"] = mane["Ensembl_Gene"].map(strip_version)
    mane["mane_ENST_ID"] = mane["Ensembl_nuc"].map(strip_version)
    dupes = mane["ENSG_ID"].duplicated().sum()
    if dupes:
        print(f"WARNING: {dupes} genes have >1 MANE Select row (unexpected)")
    return mane[["ENSG_ID", "symbol", "mane_ENST_ID"]]


def main() -> None:
    hits = pd.read_csv(HITS_CSV)
    n_hits = len(hits)

    id_map = pd.read_csv(TX_ID_MAP)[["transcript_name", "ENST_ID", "ENSG_ID", "gene_name"]]
    hits = hits.merge(id_map, on=["transcript_name", "gene_name"], how="left", suffixes=("", "_map"))
    n_no_enst = hits["ENST_ID"].isna().sum()
    print(f"{n_no_enst} / {n_hits} hits have no ENST_ID in tx_id_map (novel/long-read transcripts)")

    mane = load_mane_select()
    print(f"Loaded {len(mane)} MANE Select genes")

    hits = hits.merge(mane, on="ENSG_ID", how="left")
    n_no_mane_gene = hits["mane_ENST_ID"].isna().sum()
    print(f"{n_no_mane_gene} / {n_hits} hits belong to a gene with no MANE Select entry (non-coding biotype)")

    def tx_role(row):
        if pd.isna(row["ENST_ID"]) or pd.isna(row["mane_ENST_ID"]):
            return "no_MANE_coverage"
        return "Canonical" if row["ENST_ID"] == row["mane_ENST_ID"] else "Alternate"

    hits["tx_role_mane"] = hits.apply(tx_role, axis=1)
    assert len(hits) == n_hits, "row count changed during merge — check for duplicate keys"

    out_csv = OUT_DIR / "hit_scenarios_mane.csv"
    hits.to_csv(out_csv, index=False)
    print(f"\nSaved {len(hits)} rows -> {out_csv}")

    print("\n=== tx_role_mane distribution ===")
    print(hits["tx_role_mane"].value_counts(dropna=False).to_string())

    print("\n=== 4-scenario crosstab: tx_role_mane x usage_direction ===")
    ct = pd.crosstab(hits["tx_role_mane"], hits["usage_direction"], margins=True)
    print(ct.to_string())

    covered = hits[hits["tx_role_mane"] != "no_MANE_coverage"]
    print(f"\n=== Same crosstab, MANE-covered hits only (n={len(covered)}) ===")
    ct2 = pd.crosstab(covered["tx_role_mane"], covered["usage_direction"], margins=True)
    print(ct2.to_string())


if __name__ == "__main__":
    main()

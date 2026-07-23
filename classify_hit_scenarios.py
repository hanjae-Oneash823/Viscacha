#!/usr/bin/env python3
"""
Classify the 2,599 permutation-significant DTU hits by isoform role
(Dominant/Minor, from empirical Control-PSI rank) and usage_direction
(AD_enriched/CT_enriched, already in the hits CSV), and flag the two
priority groups:

  new_target_candidate      = Minor isoform, AD_enriched (any DGE)      -> scenarios MIEx/MIC/MIR
  trial_failure_explanation = Dominant isoform, CT_enriched (any DGE)   -> scenarios DIDi/DIDe/DIEr

These two flags depend only on isoform_role x usage_direction, NOT on
gene-level DGE direction — verified empirically: priority group counts
are IDENTICAL whether DGE is computed from long-read (LR) or short-read
(SR) pseudobulk (1093 / 1019 / 487 either way), even though the two DGE
sources barely correlate with each other (see compare_lr_sr_dge section
below: genome-wide Pearson r=0.094, zero overlapping significant genes
in every cell type). So priority flagging is used as the primary,
robust output; the full 12-scenario sub-label (which DOES depend on
DGE direction, only 38% LR/SR agreement) is retained as descriptive/
exploratory only, reported for both sources side by side rather than
picking one as ground truth.

Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/python classify_hit_scenarios.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

REPO_ROOT   = Path(__file__).resolve().parent
HITS_CSV    = REPO_ROOT / "outputs/layer1_1/DIU_significant_hits_combined.csv"
DIU_DIR     = Path("/node212data/welcome3/Grad_proj_2026/DATA/DIU_result_with_permutation")
GENE_DE_LR  = REPO_ROOT / "outputs/layer1/gene_level_de_results.csv"
GENE_DE_SR  = REPO_ROOT / "outputs/layer1/gene_level_de_results_SR.csv"
OUT_DIR     = REPO_ROOT / "outputs/scenario_analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CELL_TYPES = [
    "Astrocyte", "Excitatory_neuron", "Inhibitory_neuron", "Lymphocyte",
    "Microglia", "Oligodendrocyte", "OPC", "Vascular_cell",
]

DGE_LOG2FC_FLOOR = 0.585  # ~1.5-fold; magnitude-only, no padj gate (descriptive sub-label only)

SCENARIO_TABLE = {
    ("AD_enriched", "Dominant", "up"):   ("DIEx", "Dominant Isoform Expansion"),
    ("AD_enriched", "Minor",    "up"):   ("MIEx", "Minor Isoform Expansion"),
    ("AD_enriched", "Dominant", "down"): ("DIC",  "Dominant Isoform Concentration"),
    ("AD_enriched", "Minor",    "down"): ("MIC",  "Minor Isoform Concentration"),
    ("AD_enriched", "Dominant", "flat"): ("DIR",  "Dominant Isoform Redistribution"),
    ("AD_enriched", "Minor",    "flat"): ("MIR",  "Minor Isoform Redistribution"),
    ("CT_enriched", "Dominant", "up"):   ("DIDi", "Dominant Isoform Dilution"),
    ("CT_enriched", "Minor",    "up"):   ("MIDi", "Minor Isoform Dilution"),
    ("CT_enriched", "Dominant", "down"): ("DIDe", "Dominant Isoform Depletion"),
    ("CT_enriched", "Minor",    "down"): ("MIDe", "Minor Isoform Depletion"),
    ("CT_enriched", "Dominant", "flat"): ("DIEr", "Dominant Isoform Erosion"),
    ("CT_enriched", "Minor",    "flat"): ("MIEr", "Minor Isoform Erosion"),
}

# Priority is defined directly on (isoform_role, usage_direction) -- DGE-independent.
PRIORITY_RULES = {
    ("Minor",    "AD_enriched"): "new_target_candidate",
    ("Dominant", "CT_enriched"): "trial_failure_explanation",
}


def load_dominance() -> pd.DataFrame:
    """Rank every isoform of every gene by Control PSI, per cell type.
    Rank 1 = Dominant, everything else = Minor. Uses the full per-cell-type
    DIU files (all isoforms tested, not just significant hits), so dominance
    is empirical (Control PSI), not the (unreliable) Ensembl_canonical tag.
    """
    frames = []
    for ct in CELL_TYPES:
        fpath = DIU_DIR / f"DIU_result_with_Permutation_10000_donor_umi_cutoff_{ct}.csv"
        df = pd.read_csv(fpath, usecols=["transcript_name", "gene_name", "Control"])
        df["cell_type"] = ct
        df["dominance_rank"] = (
            df.groupby("gene_name")["Control"].rank(method="first", ascending=False)
        )
        df["isoform_role"] = np.where(df["dominance_rank"] == 1, "Dominant", "Minor")
        frames.append(df[["transcript_name", "gene_name", "cell_type", "isoform_role"]])
    return pd.concat(frames, ignore_index=True)


def classify_dge_magnitude(log2fc: float) -> str:
    if pd.isna(log2fc) or abs(log2fc) < DGE_LOG2FC_FLOOR:
        return "flat"
    return "up" if log2fc > 0 else "down"


def attach_scenario_label(hits: pd.DataFrame, gene_de_path: Path, suffix: str) -> pd.DataFrame:
    """Descriptive-only 12-scenario sub-label for one DGE source. Does not affect priority."""
    gene_de = pd.read_csv(gene_de_path).rename(
        columns={"gene_id": "gene_name", "log2FC": "gene_log2FC", "padj": "gene_padj"}
    )
    merged = hits.merge(
        gene_de[["gene_name", "cell_type", "gene_log2FC"]],
        on=["gene_name", "cell_type"], how="left",
    )
    dge_class = merged["gene_log2FC"].apply(classify_dge_magnitude)
    keys = list(zip(merged["usage_direction"], merged["isoform_role"], dge_class))
    codes, _names = zip(*(SCENARIO_TABLE.get(k, (np.nan, np.nan)) for k in keys))
    hits[f"scenario_code_{suffix}"] = codes
    return hits


def compare_de_sources(hits_keys: pd.DataFrame) -> None:
    lr = pd.read_csv(GENE_DE_LR).rename(columns={"gene_id": "gene_name", "log2FC": "log2FC_LR", "padj": "padj_LR"})
    sr = pd.read_csv(GENE_DE_SR).rename(columns={"gene_id": "gene_name", "log2FC": "log2FC_SR", "padj": "padj_SR"})
    merged = lr.merge(sr, on=["gene_name", "cell_type"], how="inner")
    print(f"\n=== LR vs SR gene-level DE comparison (context only) ===")
    print(f"Matched {len(merged)} gene x cell_type pairs tested in both")

    valid = merged.dropna(subset=["log2FC_LR", "log2FC_SR"])
    r_p, _ = pearsonr(valid["log2FC_LR"], valid["log2FC_SR"])
    print(f"All matched genes (n={len(valid)}): Pearson r={r_p:.3f}")

    hit_merged = hits_keys.merge(merged, on=["gene_name", "cell_type"], how="inner")
    valid_hits = hit_merged.dropna(subset=["log2FC_LR", "log2FC_SR"])
    if len(valid_hits) > 2:
        r_p2, _ = pearsonr(valid_hits["log2FC_LR"], valid_hits["log2FC_SR"])
        print(f"Hit gene x cell_type combos only (n={len(valid_hits)}): Pearson r={r_p2:.3f}")


def main() -> None:
    hits = pd.read_csv(HITS_CSV)
    n_hits = len(hits)

    dom = load_dominance()
    hits = hits.merge(dom, on=["transcript_name", "gene_name", "cell_type"], how="left")
    n_missing_dom = hits["isoform_role"].isna().sum()
    if n_missing_dom:
        print(f"WARNING: {n_missing_dom} hits had no dominance match")

    # Priority: DGE-independent, defined purely on isoform_role x usage_direction
    hits["priority"] = [
        PRIORITY_RULES.get((role, direction))
        for role, direction in zip(hits["isoform_role"], hits["usage_direction"])
    ]

    # Descriptive-only 12-scenario sub-labels, one column per DGE source
    hits = attach_scenario_label(hits, GENE_DE_LR, "LR")
    hits = attach_scenario_label(hits, GENE_DE_SR, "SR")
    assert len(hits) == n_hits, "row count changed during merge — check for duplicate keys"

    out_csv = OUT_DIR / "hit_scenarios.csv"
    hits.to_csv(out_csv, index=False)
    print(f"Saved {len(hits)} rows -> {out_csv}")

    print("\n=== Priority group counts (DGE-independent) ===")
    print(hits["priority"].value_counts(dropna=False).to_string())

    for label in ("new_target_candidate", "trial_failure_explanation"):
        sub = hits[hits["priority"] == label]
        print(f"\n=== {label}: top 15 by |delta_usage| ===")
        print(sub[["gene_name", "transcript_name", "cell_type", "isoform_role",
                    "usage_direction", "delta_usage", "scenario_code_LR", "scenario_code_SR"]]
              .sort_values("delta_usage", key=abs, ascending=False).head(15).to_string(index=False))

    compare_de_sources(hits[["gene_name", "cell_type"]].drop_duplicates())


if __name__ == "__main__":
    main()

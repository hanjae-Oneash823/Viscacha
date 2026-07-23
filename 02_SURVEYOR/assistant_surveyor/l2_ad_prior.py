"""L2 — AD Gene Prior.

Assigns an AD genetic prior category to each hit based on three hardcoded gene
lists. No external I/O — pure lookup.

Priority order: causal > gwas > pathway > none
"""

from __future__ import annotations

import pandas as pd

from assistant_surveyor.config import AD_CAUSAL, AD_GWAS, AD_PATHWAY


def _classify(gene: str) -> str:
    if gene in AD_CAUSAL:
        return "causal"
    if gene in AD_GWAS:
        return "gwas"
    if gene in AD_PATHWAY:
        return "pathway"
    return "none"


def run(hits: pd.DataFrame) -> pd.DataFrame:
    """Add ad_prior_category and ad_prior_flag columns. Returns new DataFrame."""
    print("[L2] Assigning AD gene priors ...", flush=True)

    cat = hits["gene_name"].map(_classify)
    result = hits.copy()
    result["ad_prior_category"] = cat
    result["ad_prior_flag"]     = cat != "none"

    prior_genes = hits.loc[result["ad_prior_flag"], "gene_name"].unique()
    print(f"[L2] {result['ad_prior_flag'].sum()} hits from {len(prior_genes)} "
          f"AD-prior genes", flush=True)
    if len(prior_genes) > 0:
        print(f"       genes: {', '.join(sorted(prior_genes))}", flush=True)

    return result

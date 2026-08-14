"""DOSSIER — data loading for a single (gene, cell_type[, transcript]) candidate.

All reads are cached in-process (functools.lru_cache) since a single CLI
invocation renders one dossier and re-reading the same wide pseudobulk CSVs
per donor row would be wasteful.
"""

from __future__ import annotations

import functools
import json

import pandas as pd

from dossier.config import HITS_CSV, MASTER_SURVEYOR_GROUPS, PFAM_CACHE, PSEUDOBULK_DIR, TX_ID_MAP_CSV


@functools.lru_cache(maxsize=1)
def load_hits() -> pd.DataFrame:
    return pd.read_csv(HITS_CSV)


@functools.lru_cache(maxsize=1)
def load_tx_id_map() -> pd.DataFrame:
    return pd.read_csv(TX_ID_MAP_CSV)


@functools.lru_cache(maxsize=1)
def load_pfam_domains() -> dict[str, list[dict]]:
    with PFAM_CACHE.open() as f:
        return json.load(f)


def domains_for(enst_id: str, kind: str) -> list[dict]:
    """kind is 'canonical' or 'alt', matching j2_pfam_hits.json's key suffix."""
    return load_pfam_domains().get(f"{enst_id}__{kind}", [])


def get_hit_rows(gene: str, cell_type: str, hit_enst: str | None = None) -> pd.DataFrame:
    """All hits_deep.csv rows for this (gene, cell_type) -- one row for
    drug_repurposing/novel_target hits, one row per ranked alt_rank for
    trial_failure hits.

    (gene_name, cell_type) is NOT always a unique hit key -- a gene/cell_type
    can carry a trial_failure hit (canonical collapses) and a separate
    drug_repurposing/novel_target hit (a different transcript rises) at the
    same time, each with its own hit_ENST_ID. Pass hit_enst to pick one; if
    omitted and the pair is ambiguous, this raises rather than silently
    mixing both hits' rows into one dossier.
    """
    df = load_hits()
    rows = df[(df["gene_name"] == gene) & (df["cell_type"] == cell_type)]
    if hit_enst is not None:
        rows = rows[rows["hit_ENST_ID"] == hit_enst]
    if rows.empty:
        raise ValueError(f"no hit found for gene={gene!r} cell_type={cell_type!r} hit_enst={hit_enst!r}")
    if hit_enst is None and rows["hit_ENST_ID"].nunique() > 1:
        options = sorted(rows["hit_ENST_ID"].unique())
        raise ValueError(
            f"gene={gene!r} cell_type={cell_type!r} has {len(options)} distinct hits "
            f"({options}) -- pass hit_enst to disambiguate"
        )
    return rows.sort_values("alt_rank") if "alt_rank" in rows.columns else rows


def resolve_hit_enst(gene: str, cell_type: str, prefer: str | None = None) -> str:
    """One concrete hit_ENST_ID for (gene, cell_type) -- `prefer` if given
    and valid, else the first (sorted) hit found. Used by the multi-cell-type
    overview sections (stacked bars, donor dots), which show one
    representative hit per cell type the gene appears in; a rare
    two-hits-per-pair collision in some OTHER cell type than the one this
    dossier is actually about must never crash the render, so this always
    resolves to something instead of raising like get_hit_rows does.
    """
    df = load_hits()
    rows = df[(df["gene_name"] == gene) & (df["cell_type"] == cell_type)]
    if rows.empty:
        raise ValueError(f"no hit found for gene={gene!r} cell_type={cell_type!r}")
    options = sorted(rows["hit_ENST_ID"].unique())
    return prefer if prefer in options else options[0]


def dossier_filename(gene: str, cell_type: str, hit_enst: str) -> str:
    """The exact filename generate_all.py writes for this hit: plain
    <gene>_<cell_type>.html, except the pairs that carry >1 hit within
    master_surveyor's scope, which get _<hit_ENST_ID> appended. Shared here
    so the index page links to the real files without re-deriving the
    ambiguous-pair set itself.
    """
    df = load_hits()
    rows = df[
        (df["gene_name"] == gene) & (df["cell_type"] == cell_type)
        & df["master_group"].isin(MASTER_SURVEYOR_GROUPS)
    ]
    suffix = f"_{hit_enst}" if rows["hit_ENST_ID"].nunique() > 1 else ""
    return f"{gene}_{cell_type}{suffix}.html"


def get_gene_cell_types(gene: str) -> list[str]:
    """Every cell type this gene appears as a hit in (any master_group),
    for the per-cell-type stacked usage bars.
    """
    df = load_hits()
    return sorted(df.loc[df["gene_name"] == gene, "cell_type"].dropna().unique().tolist())


@functools.lru_cache(maxsize=8)
def _pseudobulk_counts(cell_type: str) -> pd.DataFrame:
    path = PSEUDOBULK_DIR / f"counts_{cell_type}.csv"
    return pd.read_csv(path, index_col=0)


@functools.lru_cache(maxsize=8)
def _pseudobulk_metadata(cell_type: str) -> pd.DataFrame:
    path = PSEUDOBULK_DIR / f"metadata_{cell_type}.csv"
    return pd.read_csv(path, index_col=0)


def get_gene_isoform_usage(gene: str, cell_type: str) -> pd.DataFrame:
    """Pooled usage (%) for EVERY isoform of the gene present in this cell
    type's pseudobulk -- not just the ones J1c ranked/named in
    hits_deep.csv. Pooled = sum(raw reads across donors in condition) /
    sum(gene-total reads across donors in condition), matching the
    convention the original DTU test itself uses (verified against
    DIU_result_with_Permutation_*.csv's own AD/Control columns, which
    chi_padj and delta_usage in hits_deep.csv are computed from) -- NOT
    hits_deep.csv's alt_usage_pct_AD/control, which is a separate,
    unweighted per-donor mean computed later just for isoform ranking and
    can diverge sharply from the pooled figure when one donor's depth
    dominates. Sums to ~1.0 per condition across all isoforms of the gene,
    same as the mean-of-donors convention did -- the basis for drawing the
    un-ranked remainder of the bar in grey.

    Returns columns: transcript_name, usage_pct_AD, usage_pct_control.
    """
    counts = _pseudobulk_counts(cell_type)
    meta = _pseudobulk_metadata(cell_type)

    tx_map = load_tx_id_map()
    gene_tx_names = tx_map.loc[tx_map["gene_name"] == gene, "transcript_name"].tolist()
    gene_tx_names = [t for t in gene_tx_names if t in counts.columns]
    if not gene_tx_names:
        return pd.DataFrame(columns=["transcript_name", "usage_pct_AD", "usage_pct_control"])

    gene_total = counts[gene_tx_names].sum(axis=1)
    condition = meta["condition"]
    denom = {c: gene_total[condition == c].sum() for c in ("AD", "Control")}

    rows = []
    for tx in gene_tx_names:
        numer = {c: counts.loc[condition == c, tx].sum() for c in ("AD", "Control")}
        rows.append({
            "transcript_name": tx,
            "usage_pct_AD": (numer["AD"] / denom["AD"]) if denom["AD"] else 0.0,
            "usage_pct_control": (numer["Control"] / denom["Control"]) if denom["Control"] else 0.0,
        })
    return pd.DataFrame(rows)


def get_donor_usage(gene: str, transcript_name: str, cell_type: str) -> pd.DataFrame:
    """Per-donor raw usage for one transcript in one cell type.

    Returns columns: donor, condition, raw_count, gene_total_count, psi, n_cells.
    """
    counts = _pseudobulk_counts(cell_type)
    meta   = _pseudobulk_metadata(cell_type)

    tx_map = load_tx_id_map()
    gene_tx_names = tx_map.loc[tx_map["gene_name"] == gene, "transcript_name"].tolist()
    gene_tx_names = [t for t in gene_tx_names if t in counts.columns]
    if transcript_name not in counts.columns:
        raise ValueError(
            f"transcript {transcript_name!r} not found in counts_{cell_type}.csv"
        )

    gene_total = counts[gene_tx_names].sum(axis=1)
    raw = counts[transcript_name]
    psi = (raw / gene_total).where(gene_total > 0)

    out = pd.DataFrame({
        "donor": counts.index,
        "raw_count": raw.values,
        "gene_total_count": gene_total.values,
        "psi": psi.values,
    }).set_index("donor")
    out["condition"] = meta["condition"]
    out["n_cells"] = meta["n_cells"]
    return out.reset_index().dropna(subset=["condition"])

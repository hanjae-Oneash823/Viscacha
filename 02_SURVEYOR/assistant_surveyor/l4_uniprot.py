"""L4 — UniProt Batch Domain Presence.

Checks whether each gene's reviewed human UniProt entry has structural features
(Domain, Binding site, Active site, Repeat, Zinc finger). No sequence alignment
— gene-level presence flag only. Used by L5 to distinguish proxy Type C from D.

Results cached per batch to outputs/assistant_surveyor/cache/.
"""

from __future__ import annotations

import io
import math
import sys
from pathlib import Path

import pandas as pd

from assistant_surveyor.config import (
    CACHE_DIR, UNIPROT_BACKOFF, UNIPROT_BATCH_SIZE, UNIPROT_ENDPOINT,
    UNIPROT_RETRIES, UNIPROT_TIMEOUT,
)
from assistant_surveyor.utils.http import HTTPError, get_text

_FIELDS = (
    "accession,gene_names,protein_name,length,"
    "ft_domain,ft_binding,ft_act_site,ft_repeat,ft_zn_fing"
)

_FEAT_COLS = ["ft_domain", "ft_binding", "ft_act_site", "ft_repeat", "ft_zn_fing"]


def _query_uniprot(gene_names: list[str]) -> pd.DataFrame | None:
    """Fetch TSV for a batch of gene names. Returns DataFrame or None on error."""
    or_clause = " OR ".join(f"gene:{g}" for g in gene_names)
    query = f"({or_clause}) AND (reviewed:true) AND (organism_id:9606)"
    params = {"query": query, "fields": _FIELDS, "format": "tsv", "size": "500"}
    try:
        text = get_text(
            UNIPROT_ENDPOINT, params=params,
            timeout=UNIPROT_TIMEOUT, retries=UNIPROT_RETRIES, backoff=UNIPROT_BACKOFF,
        )
        return pd.read_csv(io.StringIO(text), sep="\t")
    except (HTTPError, Exception) as exc:
        print(f"    [L4] UniProt request failed: {exc}", flush=True)
        return None


def _extract_domain_info(row: pd.Series) -> dict:
    """Parse feature columns from a single UniProt TSV row."""
    has_domain   = pd.notna(row.get("Domain [FT]")) and str(row.get("Domain [FT]", "")).strip() != ""
    has_binding  = pd.notna(row.get("Binding site [FT]")) and str(row.get("Binding site [FT]", "")).strip() != ""
    has_act_site = pd.notna(row.get("Active site [FT]")) and str(row.get("Active site [FT]", "")).strip() != ""
    has_repeat   = pd.notna(row.get("Repeat [FT]")) and str(row.get("Repeat [FT]", "")).strip() != ""
    has_zn       = pd.notna(row.get("Zinc finger [FT]")) and str(row.get("Zinc finger [FT]", "")).strip() != ""

    has_structural = has_domain or has_binding or has_act_site or has_repeat or has_zn

    # count domains and extract names (format: "DOMAIN start..end; /note="Name";...")
    domain_str = str(row.get("Domain [FT]", "") or "")
    domain_entries = [e.strip() for e in domain_str.split(";") if "/note=" in e]
    domain_names_list = []
    for entry in domain_entries:
        # extract value between /note=" and "
        if '/note="' in entry:
            name = entry.split('/note="')[1].rstrip('"').strip()
            if name and name not in domain_names_list:
                domain_names_list.append(name)
    n_domains    = len(domain_names_list)
    domain_names = ", ".join(domain_names_list[:3])

    return {
        "has_structural_feat": has_structural,
        "n_domains":           n_domains,
        "domain_names":        domain_names,
    }


def _match_gene(df: pd.DataFrame, gene: str) -> pd.Series | None:
    """Find the best reviewed-human entry for *gene* in a UniProt TSV DataFrame."""
    if df is None or df.empty:
        return None
    gene_col = "Gene Names" if "Gene Names" in df.columns else df.columns[1]
    acc_col  = "Entry" if "Entry" in df.columns else df.columns[0]
    for _, row in df.iterrows():
        names = str(row.get(gene_col, "") or "").upper().split()
        if gene.upper() in names:
            return row
    return None


def run(hits: pd.DataFrame) -> pd.DataFrame:
    """Add UniProt domain columns. Returns new DataFrame."""
    unique_genes = sorted(hits["gene_name"].dropna().unique().tolist())
    print(f"[L4] Querying UniProt for {len(unique_genes)} unique genes in batches "
          f"of {UNIPROT_BATCH_SIZE} ...", flush=True)

    n_batches = math.ceil(len(unique_genes) / UNIPROT_BATCH_SIZE)
    gene_info: dict[str, dict] = {}

    for i in range(n_batches):
        batch = unique_genes[i * UNIPROT_BATCH_SIZE : (i + 1) * UNIPROT_BATCH_SIZE]
        print(f"  [L4] batch {i+1}/{n_batches} ({len(batch)} genes) ...",
              end=" ", flush=True)

        cache_path = CACHE_DIR / f"uniprot_batch_{i:03d}.tsv"
        if cache_path.exists():
            try:
                df_batch = pd.read_csv(cache_path, sep="\t")
                print("cached", flush=True)
            except Exception:
                df_batch = None
                print("cache read failed, re-fetching", flush=True)
        else:
            df_batch = _query_uniprot(batch)
            if df_batch is not None and not df_batch.empty:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                df_batch.to_csv(cache_path, sep="\t", index=False)
            print(f"ok ({0 if df_batch is None else len(df_batch)} entries returned)",
                  flush=True)

        for gene in batch:
            row = _match_gene(df_batch, gene)
            if row is None:
                gene_info[gene] = {
                    "uniprot_acc":       None,
                    "has_structural_feat": False,
                    "n_domains":           0,
                    "domain_names":        "",
                }
            else:
                acc_col = "Entry" if "Entry" in row.index else row.index[0]
                info = _extract_domain_info(row)
                info["uniprot_acc"] = str(row.get(acc_col, ""))
                gene_info[gene] = info

    # -- attach to hits ------------------------------------------------------
    rows = []
    for _, hit in hits.iterrows():
        info = gene_info.get(hit["gene_name"], {
            "uniprot_acc": None, "has_structural_feat": False,
            "n_domains": 0, "domain_names": "",
        })
        rows.append(info)

    l4 = pd.DataFrame(rows, index=hits.index)
    result = pd.concat([hits.reset_index(drop=True), l4.reset_index(drop=True)], axis=1)

    n_with_feat = result["has_structural_feat"].sum()
    n_with_acc  = result["uniprot_acc"].notna().sum()
    print(f"[L4] UniProt resolved: {n_with_acc} hits have accession, "
          f"{n_with_feat} have structural features", flush=True)

    return result

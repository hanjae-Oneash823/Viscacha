"""J3 — Drug target enrichment via ChEMBL + DGIdb.

For each unique gene (by UniProt accession and gene symbol):
  - ChEMBL REST: find target entry, count drug molecules, get max clinical phase
  - DGIdb GraphQL: count drug-gene interactions

New columns added:
    chembl_target_id    — ChEMBL target ID (e.g. CHEMBL2363065) or ""
    chembl_max_phase    — highest clinical phase of any approved molecule (0–4)
    n_drugs             — number of ChEMBL drug molecules with binding data
    drug_names          — pipe-separated list of top-5 drug names by max_phase
    dgidb_interactions  — count of DGIdb drug-gene interactions
    is_druggable        — bool: has a ChEMBL target entry
"""

from __future__ import annotations

import json
import sys
import time
import math
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from junior_surveyor.config import (
    CACHE_DIR, CHEMBL_BASE, CHEMBL_BACKOFF, CHEMBL_RETRIES, CHEMBL_TIMEOUT,
    DGIDB_BATCH, DGIDB_BACKOFF, DGIDB_ENDPOINT, DGIDB_RETRIES, DGIDB_TIMEOUT,
)

_CHEMBL_WORKERS = 15   # concurrent ChEMBL request threads
_progress_lock  = threading.Lock()


def _log(msg: str) -> None:
    print(f"  [j3] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# ChEMBL
# ---------------------------------------------------------------------------

def _chembl_target_by_uniprot(uniprot_acc: str) -> dict:
    """Return target info dict for a UniProt accession."""
    import requests

    url = (f"{CHEMBL_BASE}/target"
           f"?target_components__accession={uniprot_acc}"
           f"&target_type=SINGLE+PROTEIN&format=json")
    try:
        resp = requests.get(url, timeout=CHEMBL_TIMEOUT)
        if resp.status_code != 200:
            return {}
        targets = resp.json().get("targets", [])
        if not targets:
            return {}
        return targets[0]   # take first SINGLE_PROTEIN hit
    except Exception:
        return {}


def _chembl_drugs_for_target(target_id: str) -> tuple[int, int, list[str]]:
    """Return (n_drugs, max_phase, top_drug_names) for a ChEMBL target.

    Fetches mechanism IDs, then retrieves all molecule info in one batch GET
    using molecule_chembl_id__in= rather than one call per molecule.
    max_phase comes back as a float string ("4.0") so we cast via float→int.
    """
    import requests

    mech_url = (f"{CHEMBL_BASE}/mechanism"
                f"?target_chembl_id={target_id}"
                f"&format=json&limit=100")
    try:
        resp = requests.get(mech_url, timeout=CHEMBL_TIMEOUT)
        if resp.status_code != 200:
            return 0, 0, []
        mechs = resp.json().get("mechanisms", [])
        if not mechs:
            return 0, 0, []

        molecule_ids = list({m.get("molecule_chembl_id", "") for m in mechs
                             if m.get("molecule_chembl_id")})[:20]

        # Batch molecule fetch — one request for all IDs
        ids_str  = ",".join(molecule_ids)
        mol_url  = (f"{CHEMBL_BASE}/molecule"
                    f"?molecule_chembl_id__in={ids_str}"
                    f"&format=json"
                    f"&fields=molecule_chembl_id,pref_name,max_phase"
                    f"&limit={len(molecule_ids)}")
        mr = requests.get(mol_url, timeout=CHEMBL_TIMEOUT)
        if mr.status_code != 200:
            return len(molecule_ids), 0, []

        molecules  = mr.json().get("molecules", [])
        drug_info: list[tuple[int, str]] = []
        for mol in molecules:
            phase = int(float(mol.get("max_phase") or 0))
            name  = mol.get("pref_name") or mol.get("molecule_chembl_id", "")
            drug_info.append((phase, name))

        if not drug_info:
            return len(molecule_ids), 0, []

        drug_info.sort(key=lambda x: -x[0])
        return len(drug_info), drug_info[0][0], [n for _, n in drug_info[:5]]

    except Exception:
        return 0, 0, []


def _fetch_one_gene(row: dict) -> tuple[str, dict]:
    """Fetch ChEMBL target + drug info for a single {gene_name, uniprot_acc} dict.

    Returns (gene_name, result_dict).  Designed for concurrent use via
    ThreadPoolExecutor — each call is fully independent.
    """
    gene  = row["gene_name"]
    unipr = row.get("uniprot_acc", "")
    empty = {"chembl_target_id": "", "chembl_max_phase": 0,
             "n_drugs": 0, "drug_names": "", "is_druggable": False}

    if not unipr or unipr != unipr:
        return gene, empty

    target = _chembl_target_by_uniprot(unipr)
    if not target:
        return gene, empty

    target_id = target.get("target_chembl_id", "")
    n_drugs, max_phase, names = _chembl_drugs_for_target(target_id)
    return gene, {
        "chembl_target_id": target_id,
        "chembl_max_phase": max_phase,
        "n_drugs":          n_drugs,
        "drug_names":       "|".join(names),
        "is_druggable":     True,
    }


# ---------------------------------------------------------------------------
# DGIdb
# ---------------------------------------------------------------------------

def _dgidb_query(gene_names: list[str]) -> dict[str, int]:
    """Return {gene_name: interaction_count} for a batch of gene names."""
    import requests

    names_gql = json.dumps(gene_names)
    query = f"""
    {{
      genes(names: {names_gql}) {{
        nodes {{
          name
          interactions {{ drug {{ name }} }}
        }}
      }}
    }}
    """
    for attempt in range(1, DGIDB_RETRIES + 1):
        try:
            resp = requests.post(
                DGIDB_ENDPOINT,
                json={"query": query},
                timeout=DGIDB_TIMEOUT,
            )
            if resp.status_code == 200:
                nodes = resp.json()["data"]["genes"]["nodes"]
                return {n["name"]: len(n["interactions"]) for n in nodes}
            time.sleep(DGIDB_BACKOFF * attempt)
        except Exception as exc:
            _log(f"DGIdb error (attempt {attempt}): {exc}")
            time.sleep(DGIDB_BACKOFF * attempt)
    return {}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(df: pd.DataFrame) -> pd.DataFrame:
    """Enrich df with ChEMBL + DGIdb drug-target columns."""

    # Deduplicate by gene for API efficiency
    gene_df = (df[["gene_name", "uniprot_acc"]]
               .drop_duplicates("gene_name")
               .to_dict("records"))
    unique_genes = [r["gene_name"] for r in gene_df]
    _log(f"{len(unique_genes):,} unique genes to query")

    # ── ChEMBL ──────────────────────────────────────────────────────────────
    chembl_cache = CACHE_DIR / "j3_chembl.json"
    if chembl_cache.exists():
        _log("ChEMBL: loading from cache")
        chembl_map: dict[str, dict] = json.loads(chembl_cache.read_text())
    else:
        _log(f"querying ChEMBL ({_CHEMBL_WORKERS} workers) …")
        chembl_map = {}
        done = 0
        total = len(gene_df)

        with ThreadPoolExecutor(max_workers=_CHEMBL_WORKERS) as pool:
            futures = {pool.submit(_fetch_one_gene, row): row["gene_name"]
                       for row in gene_df}
            for fut in as_completed(futures):
                gene, result = fut.result()
                chembl_map[gene] = result
                done += 1
                if done % 50 == 0 or done == total:
                    _log(f"  ChEMBL {done:,}/{total:,} genes")

        chembl_cache.parent.mkdir(parents=True, exist_ok=True)
        chembl_cache.write_text(json.dumps(chembl_map))

    n_druggable = sum(1 for v in chembl_map.values() if v.get("is_druggable"))
    _log(f"ChEMBL: {n_druggable:,}/{len(chembl_map):,} genes have a target entry")

    # ── DGIdb ───────────────────────────────────────────────────────────────
    dgidb_cache = CACHE_DIR / "j3_dgidb.json"
    if dgidb_cache.exists():
        _log("DGIdb: loading from cache")
        dgidb_map: dict[str, int] = json.loads(dgidb_cache.read_text())
    else:
        _log("querying DGIdb …")
        dgidb_map = {}
        n_batches = math.ceil(len(unique_genes) / DGIDB_BATCH)
        for i in range(n_batches):
            batch = unique_genes[i * DGIDB_BATCH:(i + 1) * DGIDB_BATCH]
            dgidb_map.update(_dgidb_query(batch))
            _log(f"  DGIdb batch {i+1}/{n_batches}")
            time.sleep(0.5)
        dgidb_cache.write_text(json.dumps(dgidb_map))

    n_dgidb = sum(1 for v in dgidb_map.values() if v > 0)
    _log(f"DGIdb: {n_dgidb:,} genes have interaction records")

    # ── Attach to dataframe ──────────────────────────────────────────────────
    for col in ["chembl_target_id", "chembl_max_phase", "n_drugs",
                "drug_names", "is_druggable"]:
        df[col] = df["gene_name"].map(
            lambda g, _col=col: chembl_map.get(g, {}).get(_col,
                False if _col == "is_druggable" else (0 if _col in {"chembl_max_phase","n_drugs"} else ""))
        )

    df["dgidb_interactions"] = df["gene_name"].map(
        lambda g: dgidb_map.get(g, 0)
    )

    return df

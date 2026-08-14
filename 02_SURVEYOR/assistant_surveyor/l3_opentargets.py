"""L3 — OpenTargets Batch Disease Query.

Fetches AD (+ PD, FTD, ALS) association scores for all unique genes in the
hit table via the disease-side GraphQL query with a gene-ID batch filter.

Key lesson from layer2/m05_disease.py: use the DISEASE-side query with the
`Bs` filter param. The target-side query with efoIds filter silently returns
0 rows and must NOT be used.

Responses are cached per-batch to outputs/assistant_surveyor/cache/.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import pandas as pd

from assistant_surveyor.config import (
    CACHE_DIR, OT_BACKOFF, OT_BATCH_SIZE, OT_EFO, OT_ENDPOINT,
    OT_LABEL_EMERGING, OT_LABEL_SUPPORTED, OT_RETRIES, OT_TIMEOUT,
)
from assistant_surveyor.utils.http import HTTPError, post_json


def _ot_label(score: float) -> str:
    if score >= OT_LABEL_SUPPORTED:
        return "supported"
    if score >= OT_LABEL_EMERGING:
        return "emerging"
    return "novel"


def _build_query(ensg_ids: list[str]) -> str:
    """Aliased 4-disease query filtered to a batch of ENSG IDs."""
    ids_json = json.dumps(ensg_ids)

    def block(alias: str, efo_id: str, include_dt: bool) -> str:
        dt_fragment = " datatypeScores { id score }" if include_dt else ""
        return (
            f'{alias}: disease(efoId: "{efo_id}") {{'
            f"  associatedTargets(Bs: {ids_json}, page: {{index: 0, size: {len(ensg_ids)}}}) {{"
            f"    rows {{ target {{ id approvedSymbol }} score{dt_fragment} }}"
            f"  }}"
            f"}}"
        )

    parts = [
        block("ad",  OT_EFO["AD"],  True),
        block("pd",  OT_EFO["PD"],  False),
        block("ftd", OT_EFO["FTD"], False),
        block("als", OT_EFO["ALS"], False),
    ]
    return "{ " + " ".join(parts) + " }"


def _parse_batch(data: dict) -> dict[str, dict]:
    """Parse one batch response into {ensg_id -> {ad, pd, ftd, als, datatypes}}."""
    results: dict[str, dict] = {}

    def extract(block: dict | None) -> dict[str, float]:
        rows = ((block or {}).get("associatedTargets") or {}).get("rows") or []
        return {
            r["target"]["id"]: float(r.get("score") or 0.0)
            for r in rows
            if r.get("target", {}).get("id")
        }

    def extract_dt(block: dict | None) -> dict[str, dict[str, float]]:
        rows = ((block or {}).get("associatedTargets") or {}).get("rows") or []
        out: dict[str, dict[str, float]] = {}
        for r in rows:
            ensg = (r.get("target") or {}).get("id")
            if ensg:
                out[ensg] = {
                    dt["id"]: float(dt["score"])
                    for dt in (r.get("datatypeScores") or [])
                }
        return out

    d = data.get("data") or {}
    ad_scores  = extract(d.get("ad"))
    pd_scores  = extract(d.get("pd"))
    ftd_scores = extract(d.get("ftd"))
    als_scores = extract(d.get("als"))
    dt_map     = extract_dt(d.get("ad"))

    for ensg, ad_score in ad_scores.items():
        results[ensg] = {
            "ad":  ad_score,
            "pd":  pd_scores.get(ensg, 0.0),
            "ftd": ftd_scores.get(ensg, 0.0),
            "als": als_scores.get(ensg, 0.0),
            "datatypes": dt_map.get(ensg, {}),
        }
    # genes present in other diseases but not AD
    for disease, smap in [("pd", pd_scores), ("ftd", ftd_scores), ("als", als_scores)]:
        for ensg, score in smap.items():
            if ensg not in results:
                results[ensg] = {
                    "ad": 0.0, "pd": 0.0, "ftd": 0.0, "als": 0.0, "datatypes": {},
                }
            results[ensg][disease] = max(results[ensg][disease], score)

    return results


def run(hits: pd.DataFrame) -> pd.DataFrame:
    """Add OT columns to hits. Returns new DataFrame."""
    print("[L3] Collecting unique ENSG IDs ...", flush=True)
    ensg_ids = [
        eid for eid in hits["ENSG_ID"].dropna().unique()
        if isinstance(eid, str) and eid.startswith("ENSG")
    ]
    print(f"[L3] Querying OpenTargets for {len(ensg_ids)} genes in batches of "
          f"{OT_BATCH_SIZE} ...", flush=True)

    n_batches = math.ceil(len(ensg_ids) / OT_BATCH_SIZE)
    all_scores: dict[str, dict] = {}

    for i in range(n_batches):
        batch = ensg_ids[i * OT_BATCH_SIZE : (i + 1) * OT_BATCH_SIZE]
        # Keyed by a hash of the batch's actual ENSG content, not its index --
        # the input gene list changes across pipeline reruns, so an
        # index-keyed cache would silently serve a stale, unrelated batch's
        # scores once the batch composition shifts (same bug as l4_uniprot).
        batch_hash = hashlib.md5("|".join(batch).encode()).hexdigest()[:16]
        cache_path = CACHE_DIR / f"ot_batch_{batch_hash}.json"

        print(f"  [L3] batch {i+1}/{n_batches} ({len(batch)} genes) ...",
              end=" ", flush=True)
        try:
            resp = post_json(
                OT_ENDPOINT,
                {"query": _build_query(batch)},
                cache_path=cache_path,
                timeout=OT_TIMEOUT,
                retries=OT_RETRIES,
                backoff=OT_BACKOFF,
            )
            parsed = _parse_batch(resp)
            all_scores.update(parsed)
            from_cache = cache_path.exists() and (
                # crude: if we just wrote it, size > 0 after the call
                True
            )
            print(f"ok ({len(parsed)} hits)", flush=True)
        except HTTPError as exc:
            print(f"FAILED ({exc}) — filling with 0.0", flush=True)

    # -- build per-hit columns -----------------------------------------------
    rows = []
    for _, hit in hits.iterrows():
        ensg = hit.get("ENSG_ID", "")
        sc   = all_scores.get(ensg, {})
        ad   = sc.get("ad", 0.0)
        pd_  = sc.get("pd", 0.0)
        ftd  = sc.get("ftd", 0.0)
        als  = sc.get("als", 0.0)
        dts  = sc.get("datatypes", {})
        mean_other = (pd_ + ftd + als) / 3.0
        rows.append({
            "ad_ot_score":   round(ad, 4),
            "ad_ot_label":   _ot_label(ad),
            "ad_ot_genetic": round(dts.get("genetic_association", 0.0), 4),
            "ad_ot_lit":     round(dts.get("literature", 0.0), 4),
            "ad_specificity": round(ad / (mean_other + 0.01), 3),
        })

    l3 = pd.DataFrame(rows, index=hits.index)
    result = pd.concat([hits.reset_index(drop=True), l3.reset_index(drop=True)], axis=1)

    n_supported = (result["ad_ot_label"] == "supported").sum()
    n_emerging  = (result["ad_ot_label"] == "emerging").sum()
    print(f"[L3] OT labels — supported: {n_supported}, emerging: {n_emerging}, "
          f"novel: {len(result) - n_supported - n_emerging}", flush=True)

    return result

"""FDA Drugs@FDA validation for target-linked drug records.

The NDC Directory is intentionally not used: FDA states that NDC assignment
does not establish approval.  This module queries the official openFDA
Drugs@FDA endpoint and stores results in a small, reusable local cache.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

import requests

from master_surveyor.config import CACHE_DIR

_API = "https://api.fda.gov/drug/drugsfda.json"
_CACHE = CACHE_DIR / "fda_drugsfda_name_cache.json"


def _key(name: str) -> str:
    return " ".join((name or "").lower().split())


def _load_cache() -> dict[str, bool]:
    try:
        return {str(k): bool(v) for k, v in json.loads(_CACHE.read_text()).items()}
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _lookup(name: str) -> tuple[str, bool]:
    """Verify a drug name against an official FDA-approved application."""
    key = _key(name)
    # Drugs@FDA records product names and active ingredients.  Searching both
    # catches common cases where a target database supplies a salt form.
    queries = [f'products.brand_name.exact:"{name}"',
               f'products.active_ingredients.name.exact:"{name}"']
    for query in queries:
        try:
            response = requests.get(f"{_API}?search={quote(query)}&limit=1", timeout=4)
            if response.status_code == 200 and response.json().get("results"):
                return key, True
        except requests.RequestException:
            continue
    return key, False


def approved_names(names: set[str]) -> set[str]:
    """Return source names independently verified in FDA Drugs@FDA.

    Results, including negative lookups, are cached by normalised name so plot
    regeneration does not repeatedly call the FDA service.
    """
    cache = _load_cache()
    todo = [name for name in names if _key(name) not in cache]
    if todo:
        with ThreadPoolExecutor(max_workers=12) as pool:
            for key, approved in pool.map(_lookup, todo):
                cache[key] = approved
        _CACHE.parent.mkdir(parents=True, exist_ok=True)
        _CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n")
    return {name for name in names if cache.get(_key(name), False)}

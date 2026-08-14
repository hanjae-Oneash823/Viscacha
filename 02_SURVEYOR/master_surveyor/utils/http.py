"""Lightweight HTTP client with JSON file caching and retry backoff.

Same pattern as assistant_surveyor/utils/http.py and
junior_surveyor/utils/http.py -- callers pass a cache_path; if the file
exists it is returned immediately without hitting the network.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import requests


class HTTPError(Exception):
    def __init__(self, msg: str, status: int | None = None):
        super().__init__(msg)
        self.status = status


def _log(msg: str) -> None:
    print(f"  {msg}", file=sys.stderr, flush=True)


def get_json(
    url: str,
    *,
    params: dict | None = None,
    cache_path: Path | None = None,
    timeout: float = 30,
    retries: int = 3,
    backoff: float = 3.0,
) -> Any:
    if cache_path and cache_path.exists():
        return json.loads(cache_path.read_text())
    result = _request("GET", url, params=params, timeout=timeout,
                       retries=retries, backoff=backoff)
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(result))
    return result


def get_text(
    url: str,
    *,
    params: dict | None = None,
    timeout: float = 30,
    retries: int = 3,
    backoff: float = 3.0,
) -> str:
    return _request("GET", url, params=params, as_text=True,
                    timeout=timeout, retries=retries, backoff=backoff)


def _request(
    method: str,
    url: str,
    *,
    params: dict | None = None,
    json_body: Any = None,
    as_text: bool = False,
    timeout: float = 30,
    retries: int = 3,
    backoff: float = 3.0,
) -> Any:
    session = requests.Session()
    session.headers["User-Agent"] = "MASTER_SURVEYOR/1.0 (VISCACHA pipeline)"
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            resp = session.request(
                method, url, params=params, json=json_body, timeout=timeout,
            )
            if resp.status_code == 200:
                return resp.text if as_text else resp.json()
            if resp.status_code == 404:
                raise HTTPError("HTTP 404", 404)
            if resp.status_code in (429, 500, 502, 503, 504):
                last_err = HTTPError(f"HTTP {resp.status_code}", resp.status_code)
                _log(f"transient {resp.status_code} on {url} (attempt {attempt}/{retries})")
                time.sleep(backoff * attempt)
                continue
            raise HTTPError(f"HTTP {resp.status_code} for {url}", resp.status_code)
        except (requests.RequestException, ValueError) as exc:
            last_err = exc
            _log(f"{type(exc).__name__} on {url} (attempt {attempt}/{retries})")
            time.sleep(backoff * attempt)
    raise HTTPError(f"exhausted retries for {url}: {last_err}")

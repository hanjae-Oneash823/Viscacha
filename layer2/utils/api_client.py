"""
Surveyor cached API client (plan 6, 7).

A thin wrapper over requests that adds:
  - on-disk response caching keyed by MD5(method + url + params + body)
  - per-database TTL (from config.cache.ttl_days)
  - retry with backoff on transient failures
  - audit logging of every call (cache hit vs. miss, status)

GET (JSON/text) and POST (JSON, for the OpenTargets GraphQL endpoint) are
supported. Returns parsed payloads; raises APIError on unrecoverable failure so
callers can apply the plan's blocking/non-blocking failure policy.
"""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from layer2.config import CONFIG
from layer2.utils.audit import AuditLog


class APIError(Exception):
    """Raised when a request cannot be satisfied from cache or network."""
    def __init__(self, message: str, status: int | str | None = None):
        super().__init__(message)
        self.status = status


def _now() -> datetime:
    return datetime.now(timezone.utc)


class APIClient:
    def __init__(self, audit: AuditLog, config=CONFIG):
        self.audit = audit
        self.config = config
        self.cache_dir = config.cache_dir
        self.ttl_days = config.cache["ttl_days"]
        self.force_refresh = config.cache.get("force_refresh", False)
        self.timeout = config.execution["request_timeout_seconds"]
        self.max_retries = config.execution["max_retries"]
        self.backoff = config.execution["retry_backoff_seconds"]
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": "Surveyor/2.1 (VISCACHA Layer 2)"})

    # -- cache plumbing ----------------------------------------------------
    def _cache_key(self, method: str, url: str, params: dict | None,
                   body: Any) -> str:
        blob = json.dumps(
            {"m": method, "u": url, "p": params or {}, "b": body},
            sort_keys=True, default=str,
        )
        return hashlib.md5(blob.encode()).hexdigest()

    def _cache_path(self, db: str, key: str) -> Path:
        return self.cache_dir / db / f"{key}.json"

    def _read_cache(self, db: str, path: Path) -> Any | None:
        if self.force_refresh or not path.exists():
            return None
        try:
            wrapper = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            return None
        fetched = datetime.fromisoformat(wrapper["_meta"]["fetched_at"])
        ttl = self.ttl_days.get(db, 30)
        if (_now() - fetched).days > ttl:
            return None  # expired
        return wrapper["_payload"]

    def _write_cache(self, db: str, path: Path, url: str, method: str,
                     status: int, payload: Any):
        path.parent.mkdir(parents=True, exist_ok=True)
        wrapper = {
            "_meta": {"fetched_at": _now().isoformat(), "url": url,
                      "method": method, "status": status},
            "_payload": payload,
        }
        path.write_text(json.dumps(wrapper))

    # -- request core ------------------------------------------------------
    def _request(self, db: str, method: str, url: str, *,
                 params: dict | None = None, json_body: Any = None,
                 headers: dict | None = None, as_text: bool = False,
                 timeout: float | None = None, max_retries: int | None = None) -> Any:
        key = self._cache_key(method, url, params, json_body)
        path = self._cache_path(db, key)

        cached = self._read_cache(db, path)
        if cached is not None:
            self.audit.api_call(db, url, "cache", cache_hit=True)
            return cached

        eff_timeout = timeout if timeout is not None else self.timeout
        eff_retries = max_retries if max_retries is not None else self.max_retries
        last_err: Exception | None = None
        for attempt in range(1, eff_retries + 1):
            try:
                resp = self._session.request(
                    method, url, params=params, json=json_body,
                    headers=headers, timeout=eff_timeout,
                )
                if resp.status_code == 200:
                    payload = resp.text if as_text else resp.json()
                    self._write_cache(db, path, url, method, 200, payload)
                    self.audit.api_call(db, url, 200, cache_hit=False)
                    return payload
                # Retry on transient server / rate-limit codes; fail fast otherwise
                if resp.status_code in (429, 500, 502, 503, 504):
                    last_err = APIError(f"transient HTTP {resp.status_code}",
                                        resp.status_code)
                    self.audit.warning(db, f"HTTP {resp.status_code} on {url} "
                                           f"(attempt {attempt}/{eff_retries})")
                    time.sleep(self.backoff * attempt)
                    continue
                self.audit.error(db, f"non-retryable response for {url}",
                                 resp.status_code)
                raise APIError(f"HTTP {resp.status_code} for {url}",
                               resp.status_code)
            except (requests.RequestException, ValueError) as e:
                last_err = e
                self.audit.warning(db, f"{type(e).__name__} on {url} "
                                       f"(attempt {attempt}/{eff_retries})")
                time.sleep(self.backoff * attempt)

        self.audit.error(db, f"exhausted retries for {url}: {last_err}")
        raise APIError(f"exhausted retries for {url}: {last_err}")

    # -- public ------------------------------------------------------------
    def get_json(self, db: str, url: str, params: dict | None = None,
                 headers: dict | None = None, timeout: float | None = None,
                 max_retries: int | None = None) -> Any:
        h = {"Content-Type": "application/json", **(headers or {})}
        return self._request(db, "GET", url, params=params, headers=h,
                             timeout=timeout, max_retries=max_retries)

    def get_text(self, db: str, url: str, params: dict | None = None,
                 headers: dict | None = None) -> str:
        return self._request(db, "GET", url, params=params, headers=headers,
                             as_text=True)

    def post_json(self, db: str, url: str, json_body: Any,
                  headers: dict | None = None) -> Any:
        h = {"Content-Type": "application/json", **(headers or {})}
        return self._request(db, "POST", url, json_body=json_body, headers=h)

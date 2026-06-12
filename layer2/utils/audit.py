"""
Surveyor audit log (plan 3.6, 8).

Timestamped, append-only record of every API call, cache hit/miss, warning,
runtime-computed value, and database version seen during a run. Every field in
a dossier should be traceable to an entry here.

The log is held in memory and flushed to disk with .write(). Entries are also
echoed to stderr at construction time when verbose=True so a live run is
followable.
"""

from __future__ import annotations

import sys
import threading
from datetime import datetime, timezone
from pathlib import Path


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class AuditLog:
    def __init__(self, path: Path, verbose: bool = True):
        self.path = path
        self.verbose = verbose
        self._entries: list[str] = []
        self._db_versions: dict[str, str] = {}
        self._lock = threading.Lock()  # M04-M07 run on a thread pool
        self._emit(f"=== Surveyor run started {_ts()} ===", tag="RUN")

    # -- core --------------------------------------------------------------
    def _emit(self, msg: str, tag: str = "INFO"):
        line = f"{_ts()} [{tag}] {msg}"
        with self._lock:
            self._entries.append(line)
        if self.verbose:
            print(line, file=sys.stderr)

    # -- typed helpers -----------------------------------------------------
    def api_call(self, db: str, url: str, status: int | str,
                 cache_hit: bool, note: str = ""):
        tag = "CACHE" if cache_hit else "API"
        extra = f" — {note}" if note else ""
        self._emit(f"{db}: HTTP {status} {url}{extra}", tag=tag)

    def db_version(self, db: str, version: str):
        with self._lock:
            self._db_versions[db] = version
        self._emit(f"{db} version: {version}", tag="VERSION")

    def warning(self, where: str, msg: str):
        self._emit(f"{where}: {msg}", tag="WARN")

    def error(self, where: str, msg: str, status: int | str | None = None):
        s = f" (HTTP {status})" if status is not None else ""
        self._emit(f"{where}: {msg}{s}", tag="ERROR")

    def value(self, name: str, value):
        self._emit(f"{name} = {value}", tag="VALUE")

    def info(self, msg: str):
        self._emit(msg, tag="INFO")

    # -- accessors ---------------------------------------------------------
    @property
    def database_versions(self) -> dict[str, str]:
        with self._lock:
            return dict(self._db_versions)

    def write(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._emit(f"=== Surveyor run finished {_ts()} ===", tag="RUN")
        with self._lock:
            self.path.write_text("\n".join(self._entries) + "\n")
        if self.verbose:
            print(f"[audit] written to {self.path}", file=sys.stderr)

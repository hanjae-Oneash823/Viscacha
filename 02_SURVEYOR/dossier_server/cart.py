"""Single global cart, JSON-file-backed (outputs/master_surveyor/cart.json).

No multi-user sessions -- this is a single-user internal tool, matching the
codebase's existing file-based cache/state convention. A threading.Lock
guards read-modify-write since Flask's dev server is multi-threaded by
default.
"""

from __future__ import annotations

import json
import threading

from dossier_server.config import CART_JSON

_lock = threading.Lock()


def _key(gene: str, cell_type: str, hit_enst: str) -> str:
    return f"{gene}|{cell_type}|{hit_enst}"


def _load() -> dict:
    if CART_JSON.exists():
        return json.loads(CART_JSON.read_text())
    return {"items": []}


def _save(cart: dict) -> None:
    CART_JSON.parent.mkdir(parents=True, exist_ok=True)
    CART_JSON.write_text(json.dumps(cart, indent=2))


def list_items() -> list[dict]:
    with _lock:
        return _load()["items"]


def add_item(gene: str, cell_type: str, hit_enst: str, selected_drugs: list[str]) -> list[dict]:
    """Upserts by (gene, cell_type, hit_enst) -- re-adding the same hit with
    a different drug selection replaces the previous selection rather than
    creating a duplicate cart row.
    """
    with _lock:
        cart = _load()
        key = _key(gene, cell_type, hit_enst)
        cart["items"] = [i for i in cart["items"] if _key(i["gene"], i["cell_type"], i["hit_enst"]) != key]
        cart["items"].append({
            "gene": gene, "cell_type": cell_type, "hit_enst": hit_enst,
            "selected_drugs": selected_drugs,
        })
        _save(cart)
        return cart["items"]


def remove_item(gene: str, cell_type: str, hit_enst: str) -> list[dict]:
    with _lock:
        cart = _load()
        key = _key(gene, cell_type, hit_enst)
        cart["items"] = [i for i in cart["items"] if _key(i["gene"], i["cell_type"], i["hit_enst"]) != key]
        _save(cart)
        return cart["items"]


def clear() -> None:
    with _lock:
        _save({"items": []})

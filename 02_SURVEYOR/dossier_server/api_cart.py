"""Cart CRUD -- add/remove/list the hits (and per-hit selected drugs) the
user has picked for export.
"""

from __future__ import annotations

from flask import Blueprint, jsonify, request

from dossier_server import cart

bp = Blueprint("api_cart", __name__)


@bp.route("/api/cart", methods=["GET"])
def get_cart():
    return jsonify({"items": cart.list_items()})


@bp.route("/api/cart/items", methods=["POST"])
def add_cart_item():
    body = request.get_json(force=True)
    for field in ("gene", "cell_type", "hit_enst"):
        if not body.get(field):
            return jsonify({"error": f"missing required field: {field}"}), 400
    items = cart.add_item(
        body["gene"], body["cell_type"], body["hit_enst"],
        body.get("selected_drugs", []),
    )
    return jsonify({"items": items})


@bp.route("/api/cart/items", methods=["DELETE"])
def remove_cart_item():
    body = request.get_json(force=True)
    for field in ("gene", "cell_type", "hit_enst"):
        if not body.get(field):
            return jsonify({"error": f"missing required field: {field}"}), 400
    items = cart.remove_item(body["gene"], body["cell_type"], body["hit_enst"])
    return jsonify({"items": items})


@bp.route("/api/cart", methods=["DELETE"])
def clear_cart():
    cart.clear()
    return jsonify({"items": []})

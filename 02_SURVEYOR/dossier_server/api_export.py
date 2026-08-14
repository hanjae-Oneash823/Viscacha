"""POST /api/export (snapshot the cart, submit a background job),
GET /api/jobs/<id> (poll status), GET /api/export/<id>/download/<idx> (zip
one exported hit folder for transfer off this box -- Discovery Studio runs
on the user's own desktop, not here).
"""

from __future__ import annotations

from flask import Blueprint, jsonify, request, send_file

from dossier_server import cart, jobs
from master_surveyor.m3_export import zip_folder

bp = Blueprint("api_export", __name__)


@bp.route("/api/export", methods=["POST"])
def start_export():
    items = cart.list_items()
    if not items:
        return jsonify({"error": "cart is empty"}), 400
    body = request.get_json(silent=True) or {}
    thresholds_used = body.get("thresholds_used", {})
    job_id = jobs.submit_export_job(items, thresholds_used)
    return jsonify({"job_id": job_id})


@bp.route("/api/jobs/<job_id>")
def job_status(job_id: str):
    status = jobs.get_job(job_id)
    if status is None:
        return jsonify({"error": "job not found"}), 404
    return jsonify(status)


@bp.route("/api/export/<job_id>/download/<int:item_idx>")
def download_item(job_id: str, item_idx: int):
    status = jobs.get_job(job_id)
    if status is None:
        return jsonify({"error": "job not found"}), 404
    try:
        item = status["items"][item_idx]
    except IndexError:
        return jsonify({"error": "item index out of range"}), 404
    result = item.get("result")
    if not result or item.get("item_status") != "done":
        return jsonify({"error": "this item has not completed export yet"}), 409

    from pathlib import Path
    folder = Path(result["folder"])
    zip_path = zip_folder(folder)
    return send_file(zip_path, as_attachment=True, download_name=zip_path.name)

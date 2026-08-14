"""DOSSIER_SERVER — Flask app factory + entry point.

Usage:
    python -m dossier_server.app
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from flask import Flask, request

from dossier_server import api_cart, api_export, api_hits, api_structure
from dossier_server.config import HOST, PORT
from dossier_server.detail_page import render_detail
from dossier_server.index_page import INDEX_HTML


def create_app() -> Flask:
    app = Flask(__name__)
    app.register_blueprint(api_hits.bp)
    app.register_blueprint(api_cart.bp)
    app.register_blueprint(api_export.bp)
    app.register_blueprint(api_structure.bp)

    @app.route("/")
    def index():
        return INDEX_HTML

    @app.route("/hit/<gene>/<cell_type>")
    def hit_detail_page(gene: str, cell_type: str):
        hit_enst = request.args.get("hit_enst")
        try:
            return render_detail(gene, cell_type, hit_enst)
        except ValueError as exc:
            return str(exc), 404

    return app


if __name__ == "__main__":
    create_app().run(host=HOST, port=PORT, threaded=True)

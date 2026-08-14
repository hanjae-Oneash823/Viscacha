"""GET /api/structures/<gene>/<cell_type>?hit_enst=... -- lists which
structures (canonical + every ranked alt-isoform, per fold method) are
already available directly from STRUCTURE_CACHE_DIR/<seq_hash>/. No export
required: the viewer can show anything ColabFold/ESMFold has produced the
moment it's cached, independent of the cart/export flow.

GET /api/structure_pdb?seq_hash=...&method=...&file=... -- fetches one
selected structure's raw PDB text.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from flask import Blueprint, jsonify, request

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier import data
from master_surveyor import m2_structures
from master_surveyor.config import STRUCTURE_CACHE_DIR

bp = Blueprint("api_structure", __name__)


def _colabfold_label(path: Path) -> str:
    tag = path.name.split("_unrelaxed_")[-1].replace(".pdb", "")
    return f"ColabFold ({tag})"


def _canonical_options(seq: str | float) -> list[dict]:
    if not isinstance(seq, str) or not seq:
        return []
    h = m2_structures.seq_hash(seq)
    out_dir = STRUCTURE_CACHE_DIR / h
    options = []
    if (out_dir / "canonical_afdb.pdb").exists():
        options.append({"seq_hash": h, "method": "afdb", "label": "AlphaFold DB"})
    cf_dir = out_dir / "canonical_colabfold"
    if cf_dir.exists():
        for p in sorted(cf_dir.glob("*_unrelaxed_rank_*.pdb")):
            options.append({"seq_hash": h, "method": "colabfold", "file": p.name, "label": _colabfold_label(p)})
    if (out_dir / "esmfold.pdb").exists():
        options.append({"seq_hash": h, "method": "esmfold", "label": "ESMFold"})
    return options


def _alt_options(seq: str | float) -> list[dict]:
    if not isinstance(seq, str) or not seq:
        return []
    h = m2_structures.seq_hash(seq)
    out_dir = STRUCTURE_CACHE_DIR / h
    options = []
    cf_dir = out_dir / "alt_colabfold"
    if cf_dir.exists():
        for p in sorted(cf_dir.glob("*_unrelaxed_rank_*.pdb")):
            options.append({"seq_hash": h, "method": "colabfold", "file": p.name, "label": _colabfold_label(p)})
    if (out_dir / "esmfold.pdb").exists():
        options.append({"seq_hash": h, "method": "esmfold", "label": "ESMFold"})
    return options


@bp.route("/api/structures/<gene>/<cell_type>")
def list_structures(gene: str, cell_type: str):
    hit_enst = request.args.get("hit_enst")
    try:
        rows = data.get_hit_rows(gene, cell_type, hit_enst)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404

    first = rows.iloc[0]
    canonical = {
        "uniprot_acc": first.get("uniprot_acc"),
        "options": _canonical_options(first.get("canonical_protein_seq")),
    }

    alt_isoforms = []
    for _, r in rows.iterrows():
        options = _alt_options(r.get("alt_protein_seq"))
        if not options:
            continue
        alt_isoforms.append({
            "alt_rank": int(r["alt_rank"]) if pd.notna(r.get("alt_rank")) else None,
            "alt_ENST_ID": r.get("alt_ENST_ID"),
            "transcript_name": r.get("alt_transcript_name"),
            "is_gate_driver": bool(r.get("alt_is_gate_driver", False)),
            "usage_delta": r.get("alt_usage_delta") if pd.notna(r.get("alt_usage_delta")) else None,
            "options": options,
        })

    return jsonify({"gene": gene, "cell_type": cell_type, "canonical": canonical, "alt_isoforms": alt_isoforms})


@bp.route("/api/structure_pdb")
def get_structure_pdb():
    seq_hash = request.args.get("seq_hash")
    method = request.args.get("method")
    file = request.args.get("file")
    if not seq_hash or not method:
        return jsonify({"error": "seq_hash and method are required"}), 400

    out_dir = STRUCTURE_CACHE_DIR / seq_hash
    if method == "afdb":
        path = out_dir / "canonical_afdb.pdb"
    elif method == "esmfold":
        path = out_dir / "esmfold.pdb"
    elif method == "colabfold" and file:
        # file's containing dir (canonical_colabfold vs alt_colabfold) isn't
        # known from seq_hash alone, so both are checked.
        path = next(
            (c for c in (out_dir / "canonical_colabfold" / file, out_dir / "alt_colabfold" / file) if c.exists()),
            None,
        )
    else:
        path = None

    if path is None or not path.exists():
        return jsonify({"error": "structure not found"}), 404
    return jsonify({"pdb": path.read_text()})

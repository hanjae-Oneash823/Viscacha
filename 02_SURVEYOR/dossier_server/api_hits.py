"""GET /api/hits (filtered/sorted list) and GET /api/hits/<gene>/<cell_type>
(single hit detail), backed by m0_select's live threshold filtering and
dossier/manifest.py's shared record formatting.
"""

from __future__ import annotations

import sys
from pathlib import Path

from flask import Blueprint, jsonify, request

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier import manifest as manifest_mod
from dossier.config import MASTER_SURVEYOR_GROUPS
from dossier.data import load_hits
from master_surveyor import m0_select, m2_structures
from master_surveyor.config import STRUCTURE_CACHE_DIR

bp = Blueprint("api_hits", __name__)


def _seq_folded(seq, *, kind: str) -> dict:
    """kind is 'canonical' or 'alt' -- selects which ColabFold subdir to check.
    ESMFold shares one cache slot per sequence hash regardless of kind, so no
    kind distinction is needed there (canonical and alt sequences hash
    differently, so there's no collision).
    """
    if not isinstance(seq, str) or not seq:
        return {"colabfold": False, "esmfold": False}
    out_dir = STRUCTURE_CACHE_DIR / m2_structures.seq_hash(seq)
    if kind == "canonical":
        colabfold = (out_dir / "canonical_afdb.pdb").exists() or (
            (out_dir / "canonical_colabfold").exists()
            and any((out_dir / "canonical_colabfold").glob("*_unrelaxed_rank_*.pdb"))
        )
    else:
        colabfold = (out_dir / "alt_colabfold").exists() and any((out_dir / "alt_colabfold").glob("*_unrelaxed_rank_*.pdb"))
    esmfold = (out_dir / "esmfold.pdb").exists()
    return {"colabfold": colabfold, "esmfold": esmfold}


def _alt_fold_counts_by_hit(scoped_df) -> dict[tuple, dict]:
    """(gene_name, cell_type, hit_ENST_ID) -> {total, colabfold, esmfold} --
    total/numerators counted over the same non-canonical ranked-alt rows
    manifest.n_alts_by_hit uses (excludes the is_canonical=True row that
    trial_failure's long format carries for the canonical transcript's own
    usage rank), so this denominator always matches the n_alts already shown.
    """
    out: dict[tuple, dict] = {}
    for (gene, ct, hit_enst), g in scoped_df.groupby(["gene_name", "cell_type", "hit_ENST_ID"]):
        is_tf = g["master_group"].iloc[0] == "trial_failure_candidate"
        rows = g[g["is_canonical"] != True] if is_tf else g  # noqa: E712
        cf = esm = 0
        for _, r in rows.iterrows():
            status = _seq_folded(r.get("alt_protein_seq"), kind="alt")
            cf += int(status["colabfold"])
            esm += int(status["esmfold"])
        out[(gene, ct, hit_enst)] = {"total": len(rows), "colabfold": cf, "esmfold": esm}
    return out


def _thresholds_from_request() -> dict:
    args = request.args
    thresholds = {}
    if "tf_min_abs_delta_usage" in args:
        thresholds["tf_min_abs_delta_usage"] = float(args["tf_min_abs_delta_usage"])
    if "tf_require_domain_overlap" in args:
        thresholds["tf_require_domain_overlap"] = args["tf_require_domain_overlap"].lower() == "true"
    if "dr_min_chembl_or_ot_phase" in args:
        thresholds["dr_min_chembl_or_ot_phase"] = int(args["dr_min_chembl_or_ot_phase"])
    if "dr_require_structural_change" in args:
        thresholds["dr_require_structural_change"] = args["dr_require_structural_change"].lower() == "true"
    return thresholds


@bp.route("/api/hits")
def list_hits():
    thresholds = _thresholds_from_request()
    raw = load_hits()
    scoped = raw[raw["master_group"].isin(MASTER_SURVEYOR_GROUPS)]
    n_alts_map = manifest_mod.n_alts_by_hit(scoped)
    alt_fold_map = _alt_fold_counts_by_hit(scoped)

    filtered = m0_select.run(raw, thresholds)   # scope filter + collapse + thresholds, one row per hit

    search = (request.args.get("q") or "").strip().lower()
    if search:
        filtered = filtered[filtered["gene_name"].str.lower().str.contains(search)]

    records = []
    for _, r in filtered.iterrows():
        key = (r["gene_name"], r["cell_type"], r["hit_ENST_ID"])
        record = manifest_mod.format_record(r, n_alts_map.get(key, 1), f"{r['gene_name']}_{r['cell_type']}.html")
        canon = _seq_folded(r.get("canonical_protein_seq"), kind="canonical")
        record["canonical_colabfold"] = canon["colabfold"]
        record["canonical_esmfold"] = canon["esmfold"]
        alt_counts = alt_fold_map.get(key, {"total": 1, "colabfold": 0, "esmfold": 0})
        record["alt_total"] = alt_counts["total"]
        record["alt_colabfold_folded"] = alt_counts["colabfold"]
        record["alt_esmfold_folded"] = alt_counts["esmfold"]
        records.append(record)

    sort_key = request.args.get("sort")
    if sort_key and records and sort_key in records[0]:
        reverse = request.args.get("dir", "desc") == "desc"
        records.sort(key=lambda r: (r.get(sort_key) is None, r.get(sort_key)), reverse=reverse)

    return jsonify({"hits": records, "total": len(records), "thresholds_used": {**m0_select.DEFAULT_THRESHOLDS, **thresholds}})


@bp.route("/api/hits/<gene>/<cell_type>")
def hit_detail(gene: str, cell_type: str):
    from dossier.data import get_hit_rows

    hit_enst = request.args.get("hit_enst")
    try:
        rows = get_hit_rows(gene, cell_type, hit_enst)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 404

    row = m0_select.representative_row(rows)
    record = manifest_mod.format_record(row, len(rows), f"{gene}_{cell_type}.html")
    record["changed_aa_start"] = int(row.get("changed_aa_start") or 0)
    record["changed_aa_end"] = int(row.get("changed_aa_end") or 0)
    record["affected_domain"] = row.get("affected_domain")
    record["uniprot_acc"] = row.get("uniprot_acc")
    return jsonify(record)

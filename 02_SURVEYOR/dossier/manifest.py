"""Shared "format one hit row into a display record" logic.

Extracted out of generate_index.py's build_manifest() so dossier_server's
live /api/hits endpoint can produce identically-shaped JSON records without
duplicating this formatting a second time (it previously only lived inline
in generate_index.py's loop body).
"""

from __future__ import annotations

import math

import pandas as pd


def has_text(val) -> bool:
    return pd.notna(val) and str(val).strip() != ""


def drug_evidence(row: pd.Series) -> dict[str, bool]:
    return {
        "ChEMBL": bool(has_text(row.get("drug_names")) or (row.get("chembl_max_phase") or 0) >= 1),
        "Open Targets": bool(has_text(row.get("ot_drug_names")) or (row.get("ot_max_phase") or 0) >= 1),
        "Pharos": bool((row.get("pharos_n_drugs") or 0) > 0),
        "DGIdb": bool((row.get("dgidb_interactions") or 0) > 0),
    }


def change_bucket(length_diff: float) -> str:
    if length_diff < 0:
        return "shorter"
    if length_diff > 0:
        return "longer"
    return "same"


def format_record(r: pd.Series, n_alts: int, filename: str) -> dict:
    """r is one representative row for a hit (m0_select.representative_row's
    output, or generate_index.py's equivalent) -- n_alts and filename are
    passed in separately since both depend on information (the full ranked-
    alt group, or the ambiguous-pair filename-suffix rule) that isn't on the
    representative row alone.
    """
    evidence = drug_evidence(r)
    padj = float(r.get("chi_padj")) if pd.notna(r.get("chi_padj")) else None
    return {
        "gene": r["gene_name"],
        "cell_type": r["cell_type"],
        "hit_enst": r["hit_ENST_ID"],
        "hit_transcript": r["hit_transcript_name"],
        "master_group": r["master_group"],
        "protein_change_type": r.get("protein_change_type") or "no_sequence",
        "change_bucket": change_bucket(r.get("protein_length_diff", 0) or 0),
        "padj": padj,
        "neg_log10_padj": round(min(-math.log10(padj), 60), 2) if padj and padj > 0 else 0,
        "control_pct": round(float(r.get("Control") or 0) * 100, 2),
        "ad_pct": round(float(r.get("AD") or 0) * 100, 2),
        "delta_pp": round((float(r.get("AD") or 0) - float(r.get("Control") or 0)) * 100, 2),
        "n_alts": n_alts,
        "domains_lost": r.get("domains_lost") if has_text(r.get("domains_lost")) else None,
        "domains_gained": r.get("domains_gained") if has_text(r.get("domains_gained")) else None,
        "evidence": evidence,
        "has_evidence": any(evidence.values()),
        # protein-level
        "protein_length_diff": int(r.get("protein_length_diff") or 0),
        "pct_identity": round(float(r.get("pct_identity")), 4) if pd.notna(r.get("pct_identity")) else None,
        "changed_aa_fraction": round(float(r.get("changed_aa_fraction")), 4) if pd.notna(r.get("changed_aa_fraction")) else None,
        "n_domains": int(r.get("n_domains") or 0),
        # drug database values -- 0 is a real, meaningful "no evidence" value here, not missing data
        "chembl_max_phase": int(r.get("chembl_max_phase") or 0),
        "n_drugs": int(r.get("n_drugs") or 0),
        "chembl_bioactive_compounds": int(r.get("chembl_bioactive_compounds") or 0),
        "chembl_best_pchembl": round(float(r.get("chembl_best_pchembl")), 2) if pd.notna(r.get("chembl_best_pchembl")) and r.get("chembl_best_pchembl") else None,
        "ot_max_phase": int(r.get("ot_max_phase") or 0),
        "ot_n_drugs": int(r.get("ot_n_drugs") or 0),
        "ot_trials_total": int(r.get("ot_trials_total") or 0),
        "ot_trials_terminated": int(r.get("ot_trials_terminated") or 0),
        "has_failed_trial": bool((r.get("ot_trials_terminated") or 0) > 0),
        "ot_trial_stop_reasons": r.get("ot_trial_stop_reasons") if has_text(r.get("ot_trial_stop_reasons")) else None,
        "ot_trial_stop_example": r.get("ot_trial_stop_example") if has_text(r.get("ot_trial_stop_example")) else None,
        "pharos_n_ligands": int(r.get("pharos_n_ligands") or 0),
        "pharos_n_drugs": int(r.get("pharos_n_drugs") or 0),
        "dgidb_interactions": int(r.get("dgidb_interactions") or 0),
        "pharos_tdl": r.get("pharos_tdl") if has_text(r.get("pharos_tdl")) else "unknown",
        # drug names, needed by the live API's drug picker (not previously
        # exposed by build_manifest, which only showed evidence booleans).
        # NaN is truthy in Python (`nan or ""` stays nan), so this must
        # check has_text explicitly rather than relying on `or ""`.
        "drug_names": [d for d in (r.get("drug_names") if has_text(r.get("drug_names")) else "").split("|") if d],
        "ot_drug_names": [d for d in (r.get("ot_drug_names") if has_text(r.get("ot_drug_names")) else "").split("|") if d],
        "file": filename,
    }


def n_alts_by_hit(scoped_df: pd.DataFrame) -> dict[tuple[str, str, str], int]:
    """(gene_name, cell_type, hit_ENST_ID) -> count of non-canonical ranked
    alts for trial_failure hits (1 for every other group) -- computed from
    the ungrouped/uncollapsed scoped df, since this count is lost once a
    hit is collapsed to its single representative row.
    """
    out: dict[tuple[str, str, str], int] = {}
    for (gene, ct, hit_enst), g in scoped_df.groupby(["gene_name", "cell_type", "hit_ENST_ID"]):
        is_tf = g["master_group"].iloc[0] == "trial_failure_candidate"
        out[(gene, ct, hit_enst)] = int((g["is_canonical"] != True).sum()) if is_tf else 1  # noqa: E712
    return out

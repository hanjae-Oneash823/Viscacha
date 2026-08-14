"""J3 — Drug/binding evidence enrichment via ChEMBL + DGIdb + Open Targets + Pharos.

For each unique gene:
  - ChEMBL REST (by UniProt accession): find target entry, then every drug
    molecule with a curated mechanism against it (no cap -- every molecule
    the API returns is kept), each carrying phase, molecule type, approval
    year, withdrawal/black-box flags, route of administration, and the
    mechanism of action + action type for that specific molecule-target
    pair. Also separately counts raw bioactivity (any assayed compound with
    a measured potency against the target, regardless of clinical status --
    a much broader recall net than "has a drug mechanism").
  - DGIdb GraphQL (by gene symbol): every drug-gene interaction, each
    carrying DGIdb's own `approved` flag and interaction type(s)
    (agonist/inhibitor/blocker/etc.).
  - Open Targets GraphQL (by ENSG_ID): drugAndClinicalCandidates -- an
    independently-curated drug/clinical-candidate list per target, not just
    a re-export of ChEMBL/DGIdb's own data. Each row now also pulls the
    drug's type, mechanism(s) of action, every indication it's been tried
    or approved for (disease name + phase reached), and drugWarnings
    (withdrawal / black-box-warning records with a reason, when present).
    Trial-stop evidence (terminated/withdrawn/suspended clinicalReports) is
    read per drug row, not pooled across every drug the gene has -- the
    previous version summed clinicalReports across all of a gene's drugs
    before checking status, so "this drug failed" could never be
    distinguished from "some other drug against this same target failed."
  - Pharos/TCRD GraphQL (by gene symbol): target development level (tdl:
    Tclin/Tchem/Tbio/Tdark) + ligandCounts, plus the actual named ligands
    flagged isdrug=true (previously only a count, no names at all).

Schema verified live (2026-07-29) via GraphQL introspection against
api.platform.opentargets.org, dgidb.org/api/graphql, and
pharos-api.ncats.io/graphql, and against the ChEMBL REST molecule/mechanism
endpoints -- field names below are not guessed.

ChEMBL bioactivity and Pharos tdl/ligand counts remain informational
evidence-gathering columns only -- they are NOT wired into J4's
drug_evidence gate. That gate's trial_failure branch requires "an actual
known drug to test against the AD-shifted alternate" (see j4_gate.py); a
bioactive probe compound or Tchem ligand is real binding evidence but not
necessarily a testable drug, so folding it into the same gate would change
what "drug_evidence" means for that branch. Kept separate so both signals
are visible without silently loosening the existing gate semantics.

New columns added (legacy aggregate columns, unchanged in shape, still feed
J4's gate and the static dossier's per-source cards):
    chembl_target_id, chembl_max_phase, n_drugs, drug_names, is_druggable,
    chembl_bioactive_compounds, chembl_best_pchembl,
    dgidb_interactions, dgidb_drug_names,
    ot_max_phase, ot_n_drugs, ot_drug_names, ot_trials_total,
    ot_trials_terminated, ot_trial_stop_reasons, ot_trial_stop_example,
    pharos_tdl, pharos_n_ligands, pharos_n_drugs

New this pass:
    drug_records_json      — the real payload: a JSON-encoded list of
                              per-drug dicts, deduped by normalized name
                              across all four sources, each carrying
                              {name, sources[], phase, drug_type,
                              mechanisms[], indications[], approved_indications[],
                              is_withdrawn, withdrawal_reasons[],
                              black_box_warning, trial_stopped,
                              trial_stop_reasons[], trial_stop_example,
                              chembl_id, route[], status}. `status` is one
                              of approved / withdrawn / failed_trial /
                              in_trials / investigational -- see
                              _derive_status() for the precedence rules.
                              No cap on the number of drugs kept per gene.
"""

from __future__ import annotations

import json
import re
import sys
import time
import math
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from junior_surveyor.config import (
    CACHE_DIR, CHEMBL_BASE, CHEMBL_BACKOFF, CHEMBL_RETRIES, CHEMBL_TIMEOUT,
    DGIDB_BATCH, DGIDB_BACKOFF, DGIDB_ENDPOINT, DGIDB_RETRIES, DGIDB_TIMEOUT,
    OT_BACKOFF, OT_DRUG_BATCH, OT_ENDPOINT, OT_RETRIES, OT_TIMEOUT,
    PHAROS_BACKOFF, PHAROS_BATCH, PHAROS_ENDPOINT, PHAROS_RETRIES, PHAROS_TIMEOUT,
)

_CHEMBL_WORKERS = 15   # concurrent ChEMBL request threads
_progress_lock  = threading.Lock()


def _log(msg: str) -> None:
    print(f"  [j3] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# ChEMBL
# ---------------------------------------------------------------------------

def _chembl_target_by_uniprot(uniprot_acc: str) -> dict:
    """Return target info dict for a UniProt accession."""
    import requests

    url = (f"{CHEMBL_BASE}/target"
           f"?target_components__accession={uniprot_acc}"
           f"&target_type=SINGLE+PROTEIN&format=json")
    try:
        resp = requests.get(url, timeout=CHEMBL_TIMEOUT)
        if resp.status_code != 200:
            return {}
        targets = resp.json().get("targets", [])
        if not targets:
            return {}
        return targets[0]   # take first SINGLE_PROTEIN hit
    except Exception:
        return {}


def _chembl_route(mol: dict) -> list[str]:
    routes = []
    if mol.get("oral"):
        routes.append("oral")
    if mol.get("parenteral"):
        routes.append("parenteral")
    if mol.get("topical"):
        routes.append("topical")
    return routes


def _chembl_drugs_for_target(target_id: str) -> tuple[int, int, list[dict]]:
    """Return (n_drugs, max_phase, drug_records) for a ChEMBL target.

    Fetches every mechanism row (no cap beyond the API's own limit=100 per
    request), keeping each molecule's mechanism_of_action + action_type
    (previously fetched then discarded, leaving only the bare molecule id),
    then retrieves full molecule detail for every distinct molecule in one
    batch GET using molecule_chembl_id__in= rather than one call per
    molecule. max_phase comes back as a float string ("4.0") so we cast via
    float->int.
    """
    import requests

    mech_url = (f"{CHEMBL_BASE}/mechanism"
                f"?target_chembl_id={target_id}"
                f"&format=json&limit=100")
    try:
        resp = requests.get(mech_url, timeout=CHEMBL_TIMEOUT)
        if resp.status_code != 200:
            return 0, 0, []
        mechs = resp.json().get("mechanisms", [])
        if not mechs:
            return 0, 0, []

        # First-seen mechanism_of_action/action_type per molecule -- a
        # molecule can have >1 mechanism row against the same target
        # (rare), we keep the first.
        mech_by_mol: dict[str, dict] = {}
        for m in mechs:
            mol_id = m.get("molecule_chembl_id", "")
            if mol_id and mol_id not in mech_by_mol:
                mech_by_mol[mol_id] = {
                    "mechanism_of_action": m.get("mechanism_of_action") or "",
                    "action_type": m.get("action_type") or "",
                }

        molecule_ids = list(mech_by_mol.keys())
        if not molecule_ids:
            return 0, 0, []

        ids_str  = ",".join(molecule_ids)
        mol_url  = (f"{CHEMBL_BASE}/molecule"
                    f"?molecule_chembl_id__in={ids_str}"
                    f"&format=json"
                    f"&fields=molecule_chembl_id,pref_name,max_phase,molecule_type,"
                    f"first_approval,withdrawn_flag,black_box_warning,"
                    f"oral,parenteral,topical,usan_stem_definition"
                    f"&limit={len(molecule_ids)}")
        mr = requests.get(mol_url, timeout=CHEMBL_TIMEOUT)
        if mr.status_code != 200:
            return len(molecule_ids), 0, []

        molecules = mr.json().get("molecules", [])
        drug_records: list[dict] = []
        for mol in molecules:
            mol_id = mol.get("molecule_chembl_id", "")
            phase  = int(float(mol.get("max_phase") or 0))
            name   = mol.get("pref_name") or mol_id
            mech   = mech_by_mol.get(mol_id, {})
            drug_records.append({
                "name":                  name,
                "chembl_id":             mol_id,
                "phase":                 phase,
                "molecule_type":         mol.get("molecule_type") or "",
                "first_approval":        mol.get("first_approval"),
                "withdrawn":             bool(mol.get("withdrawn_flag")),
                "black_box_warning":     bool(mol.get("black_box_warning")),
                "route":                 _chembl_route(mol),
                "usan_stem_definition":  mol.get("usan_stem_definition") or "",
                "mechanism_of_action":   mech.get("mechanism_of_action", ""),
                "action_type":           mech.get("action_type", ""),
            })

        if not drug_records:
            return len(molecule_ids), 0, []

        drug_records.sort(key=lambda r: -r["phase"])
        return len(drug_records), drug_records[0]["phase"], drug_records

    except Exception:
        return 0, 0, []


def _chembl_bioactivity_for_target(target_id: str) -> tuple[int, float]:
    """Return (n_bioactive_compounds, best_pchembl) for a ChEMBL target.

    Broader recall net than _chembl_drugs_for_target: any assayed compound
    with a measured potency (pchembl_value = -log10 of Ki/Kd/IC50/EC50 in
    molar), not just molecules with a curated drug mechanism.
    pchembl_value >= 5 means <=10 uM potency. Results are capped at the top
    200 most potent activities (order_by=-pchembl_value) -- a lower bound on
    the true distinct-compound count for very promiscuous targets, but
    sufficient for an informational evidence signal. This cap is about raw
    bioactivity screening hits, a different (much noisier) thing than named
    drugs, so it is untouched by the "no cap on drugs" request.
    """
    import requests

    url = (f"{CHEMBL_BASE}/activity"
           f"?target_chembl_id={target_id}"
           f"&pchembl_value__gte=5"
           f"&format=json&limit=200&order_by=-pchembl_value")
    try:
        resp = requests.get(url, timeout=CHEMBL_TIMEOUT)
        if resp.status_code != 200:
            return 0, 0.0
        activities = resp.json().get("activities", [])
        if not activities:
            return 0, 0.0
        molecule_ids = {a.get("molecule_chembl_id") for a in activities
                         if a.get("molecule_chembl_id")}
        best = float(activities[0]["pchembl_value"])
        return len(molecule_ids), round(best, 2)
    except Exception:
        return 0, 0.0


def _fetch_one_gene(row: dict) -> tuple[str, dict]:
    """Fetch ChEMBL target + drug + bioactivity info for a single
    {gene_name, uniprot_acc} dict.

    Returns (gene_name, result_dict).  Designed for concurrent use via
    ThreadPoolExecutor — each call is fully independent.
    """
    gene  = row["gene_name"]
    unipr = row.get("uniprot_acc", "")
    empty = {"chembl_target_id": "", "chembl_max_phase": 0,
             "n_drugs": 0, "drug_names": "", "is_druggable": False,
             "chembl_bioactive_compounds": 0, "chembl_best_pchembl": 0.0,
             "chembl_drug_records": []}

    if not unipr or unipr != unipr:
        return gene, empty

    target = _chembl_target_by_uniprot(unipr)
    if not target:
        return gene, empty

    target_id = target.get("target_chembl_id", "")
    n_drugs, max_phase, records = _chembl_drugs_for_target(target_id)
    n_bioactive, best_pchembl = _chembl_bioactivity_for_target(target_id)
    names = [r["name"] for r in records]
    return gene, {
        "chembl_target_id":           target_id,
        "chembl_max_phase":           max_phase,
        "n_drugs":                    n_drugs,
        "drug_names":                 "|".join(names),
        "is_druggable":               True,
        "chembl_bioactive_compounds": n_bioactive,
        "chembl_best_pchembl":        best_pchembl,
        "chembl_drug_records":        records,
    }


# ---------------------------------------------------------------------------
# DGIdb
# ---------------------------------------------------------------------------

def _dgidb_query(gene_names: list[str]) -> dict[str, dict]:
    """Return {gene_name: {"n": interaction_count, "records": [drug dicts]}}
    for a batch of gene names.

    Now also pulls DGIdb's own `approved` flag and interactionTypes
    (agonist/inhibitor/blocker/etc.) per drug -- previously only the name
    was kept. Records are deduped by case-insensitive name (merging
    interaction types across duplicate claims); "n" stays the full
    (uncapped) interaction-claim count, matching the old semantics.
    """
    import requests

    names_gql = json.dumps(gene_names)
    query = f"""
    {{
      genes(names: {names_gql}) {{
        nodes {{
          name
          interactions {{
            drug {{ name conceptId approved }}
            interactionTypes {{ type }}
          }}
        }}
      }}
    }}
    """
    for attempt in range(1, DGIDB_RETRIES + 1):
        try:
            resp = requests.post(
                DGIDB_ENDPOINT,
                json={"query": query},
                timeout=DGIDB_TIMEOUT,
            )
            if resp.status_code == 200:
                nodes = resp.json()["data"]["genes"]["nodes"]
                out = {}
                for n in nodes:
                    interactions = n["interactions"]
                    records: list[dict] = []
                    by_key: dict[str, dict] = {}
                    for i in interactions:
                        drug = i.get("drug") or {}
                        name = drug.get("name")
                        if not name:
                            continue
                        itypes = [t["type"] for t in (i.get("interactionTypes") or [])
                                  if t.get("type")]
                        key = name.strip().lower()
                        if key in by_key:
                            existing = by_key[key]
                            existing["interaction_types"] = list(dict.fromkeys(
                                existing["interaction_types"] + itypes
                            ))
                            continue
                        rec = {
                            "name":              name,
                            "concept_id":        drug.get("conceptId") or "",
                            "approved":          bool(drug.get("approved")),
                            "interaction_types": itypes,
                        }
                        by_key[key] = rec
                        records.append(rec)
                    out[n["name"]] = {"n": len(interactions), "records": records}
                return out
            time.sleep(DGIDB_BACKOFF * attempt)
        except Exception as exc:
            _log(f"DGIdb error (attempt {attempt}): {exc}")
            time.sleep(DGIDB_BACKOFF * attempt)
    return {}


# ---------------------------------------------------------------------------
# Open Targets
# ---------------------------------------------------------------------------

def _parse_ot_phase(stage: str | None) -> int:
    """Map maxClinicalStage (e.g. "PHASE_2_3", "APPROVAL") to a 0-4 scale.

    Parses digits out rather than hardcoding the full enum, since Open
    Targets' ClinicalStage values include combined stages like PHASE_1_2.
    """
    if not stage:
        return 0
    s = stage.upper()
    if "APPROV" in s:
        return 4
    nums = [int(n) for n in re.findall(r"(\d+)", s)]
    return max(nums) if nums else 0


_STOPPED_STATUSES = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}


def _ot_query(ensg_ids: list[str]) -> dict[str, dict]:
    """Return per-gene drug + clinical-trial-report summary for a batch.

    One aliased target(ensemblId: ...) block per gene in a single request,
    mirroring assistant_surveyor/l3_opentargets.py's disease-side batching.
    Each drug row now also pulls drugType, mechanismsOfAction, indications
    (every disease it's been tried/approved for, with the phase reached),
    and drugWarnings (withdrawal / black-box records with a reason).

    clinicalReports carries individual trial records *per drug row* -- read
    here per-row so "this drug's trials stopped" can be distinguished from
    "some other drug against this same target stopped" (the previous
    version pooled clinicalReports across every drug row for the gene
    before checking status, which meant per-drug trial-failure attribution
    was impossible; that pooled view is kept below as
    ot_trial_stop_reasons/ot_trial_stop_example for backward compatibility
    with the static dossier's existing gene-level trial-stop note, but the
    per-drug detail now lives on each record in ot_drug_records).
    """
    import requests

    parts = []
    for i, ensg in enumerate(ensg_ids):
        parts.append(
            f'g{i}: target(ensemblId: "{ensg}") {{ '
            f'drugAndClinicalCandidates {{ rows {{ maxClinicalStage '
            f'drug {{ name drugType description '
            f'mechanismsOfAction {{ rows {{ mechanismOfAction actionType }} }} '
            f'indications {{ rows {{ maxClinicalStage disease {{ name }} }} }} '
            f'drugWarnings {{ warningType description toxicityClass efoTerm year country }} '
            f'}} '
            f'clinicalReports {{ trialOverallStatus trialWhyStopped trialStopReasonCategories }} '
            f'}} }} }}'
        )
    query = "{ " + " ".join(parts) + " }"

    for attempt in range(1, OT_RETRIES + 1):
        try:
            resp = requests.post(OT_ENDPOINT, json={"query": query}, timeout=OT_TIMEOUT)
            if resp.status_code != 200:
                time.sleep(OT_BACKOFF * attempt)
                continue

            data = resp.json().get("data") or {}
            out: dict[str, dict] = {}
            for i, ensg in enumerate(ensg_ids):
                block = data.get(f"g{i}") or {}
                rows = (block.get("drugAndClinicalCandidates") or {}).get("rows") or []

                drug_records: list[dict] = []
                for r in rows:
                    drug = r.get("drug") or {}
                    name = drug.get("name") or ""
                    if not name:
                        continue
                    phase = _parse_ot_phase(r.get("maxClinicalStage"))

                    reports = r.get("clinicalReports") or []
                    stopped = [c for c in reports
                               if (c.get("trialOverallStatus") or "").upper() in _STOPPED_STATUSES]
                    stop_categories = list(dict.fromkeys(
                        cat for c in stopped for cat in (c.get("trialStopReasonCategories") or [])
                    ))
                    stop_example = next(
                        (c["trialWhyStopped"] for c in stopped if c.get("trialWhyStopped")), ""
                    )

                    mechanisms = [
                        {"mechanism_of_action": m.get("mechanismOfAction") or "",
                         "action_type": m.get("actionType") or ""}
                        for m in ((drug.get("mechanismsOfAction") or {}).get("rows") or [])
                    ]

                    ind_by_disease: dict[str, int] = {}
                    for ind in (drug.get("indications") or {}).get("rows") or []:
                        dname = (ind.get("disease") or {}).get("name") or ""
                        if not dname:
                            continue
                        iphase = _parse_ot_phase(ind.get("maxClinicalStage"))
                        ind_by_disease[dname] = max(ind_by_disease.get(dname, 0), iphase)
                    indications = sorted(
                        ({"disease": d, "phase": p} for d, p in ind_by_disease.items()),
                        key=lambda x: -x["phase"],
                    )

                    warnings = [
                        {
                            "warning_type":  w.get("warningType") or "",
                            "description":   w.get("description") or "",
                            "toxicity_class": w.get("toxicityClass") or "",
                            "efo_term":      w.get("efoTerm") or "",
                            "year":          w.get("year"),
                            "country":       w.get("country") or "",
                        }
                        for w in (drug.get("drugWarnings") or [])
                    ]
                    is_withdrawn = any(w["warning_type"] == "Withdrawn" for w in warnings)
                    black_box    = any(w["warning_type"] == "Black Box Warning" for w in warnings)

                    drug_records.append({
                        "name":               name,
                        "phase":              phase,
                        "drug_type":          drug.get("drugType") or "",
                        "description":        drug.get("description") or "",
                        "mechanisms":         mechanisms,
                        "indications":        indications,
                        "warnings":           warnings,
                        "is_withdrawn":       is_withdrawn,
                        "black_box_warning":  black_box,
                        "trial_stopped":      bool(stopped),
                        "trial_stop_reasons": stop_categories,
                        "trial_stop_example": stop_example,
                        "n_trials":           len(reports),
                        "n_trials_stopped":   len(stopped),
                    })

                drug_records.sort(key=lambda d: -d["phase"])
                names = list(dict.fromkeys(d["name"] for d in drug_records))

                # Pooled gene-level view, kept for backward compat with the
                # static dossier's existing _trial_stop_note().
                reports_all = [c for r in rows for c in (r.get("clinicalReports") or [])]
                stopped_all = [c for c in reports_all
                               if (c.get("trialOverallStatus") or "").upper() in _STOPPED_STATUSES]
                stop_categories_all = list(dict.fromkeys(
                    cat for c in stopped_all for cat in (c.get("trialStopReasonCategories") or [])
                ))
                stop_example_all = next(
                    (c["trialWhyStopped"] for c in stopped_all if c.get("trialWhyStopped")), ""
                )

                out[ensg] = {
                    "ot_max_phase":         drug_records[0]["phase"] if drug_records else 0,
                    "ot_n_drugs":           len(drug_records),
                    "ot_drug_names":        "|".join(names),
                    "ot_trials_total":      len(reports_all),
                    "ot_trials_terminated": len(stopped_all),
                    "ot_trial_stop_reasons": "|".join(stop_categories_all),
                    "ot_trial_stop_example": stop_example_all,
                    "ot_drug_records":      drug_records,
                }
            return out
        except Exception as exc:
            _log(f"Open Targets error (attempt {attempt}): {exc}")
            time.sleep(OT_BACKOFF * attempt)
    return {}


# ---------------------------------------------------------------------------
# Pharos / TCRD
# ---------------------------------------------------------------------------

def _pharos_query(gene_names: list[str]) -> dict[str, dict]:
    """Return {gene_name: {pharos_tdl, pharos_n_ligands, pharos_n_drugs,
    pharos_drug_records}}.

    One aliased target(q: {sym: ...}) block per gene in a single request,
    same aliased-batch pattern as _ot_query. Now also requests
    ligands(isdrug: true) -- the actual named drug-flagged ligands, which
    were previously invisible (only ligandCounts, a bare number, was kept).
    Unknown symbols come back as a null block (handled below), not an
    error.
    """
    import requests

    parts = []
    for i, gene in enumerate(gene_names):
        sym = gene.replace('"', '').replace("\\", "")
        parts.append(
            f'g{i}: target(q: {{sym: "{sym}"}}) {{ tdl ligandCounts {{ name value }} '
            f'ligands(isdrug: true) {{ name description activities {{ type value }} }} }}'
        )
    query = "{ " + " ".join(parts) + " }"

    for attempt in range(1, PHAROS_RETRIES + 1):
        try:
            resp = requests.post(PHAROS_ENDPOINT, json={"query": query}, timeout=PHAROS_TIMEOUT)
            if resp.status_code != 200:
                time.sleep(PHAROS_BACKOFF * attempt)
                continue

            data = resp.json().get("data") or {}
            out: dict[str, dict] = {}
            for i, gene in enumerate(gene_names):
                block = data.get(f"g{i}")
                if not block:
                    out[gene] = {"pharos_tdl": "", "pharos_n_ligands": 0, "pharos_n_drugs": 0,
                                 "pharos_drug_records": []}
                    continue
                counts = {c["name"]: c["value"] for c in (block.get("ligandCounts") or [])}
                records = []
                for lig in (block.get("ligands") or []):
                    name = lig.get("name")
                    if not name:
                        continue
                    activity_values = [a.get("value") for a in (lig.get("activities") or [])
                                        if isinstance(a.get("value"), (int, float))]
                    records.append({
                        "name":          name,
                        "description":   lig.get("description") or "",
                        "best_activity": max(activity_values) if activity_values else None,
                    })
                out[gene] = {
                    "pharos_tdl":          block.get("tdl") or "",
                    "pharos_n_ligands":    counts.get("ligand", 0),
                    "pharos_n_drugs":      counts.get("drug", 0),
                    "pharos_drug_records": records,
                }
            return out
        except Exception as exc:
            _log(f"Pharos error (attempt {attempt}): {exc}")
            time.sleep(PHAROS_BACKOFF * attempt)
    return {}


# ---------------------------------------------------------------------------
# Cross-source merge + status derivation
# ---------------------------------------------------------------------------

STATUS_APPROVED       = "approved"
STATUS_WITHDRAWN      = "withdrawn"
STATUS_FAILED_TRIAL   = "failed_trial"
STATUS_IN_TRIALS      = "in_trials"
STATUS_INVESTIGATIONAL = "investigational"


def _norm_name(name: str) -> str:
    """Exact-match dedup key: lowercase + collapsed whitespace. Deliberately
    does NOT strip salts ("donepezil hydrochloride" vs "donepezil") -- doing
    that safely needs a real chemistry-aware normalizer, and a wrong merge
    (different salts, different profiles) is worse than a duplicate card."""
    return re.sub(r"\s+", " ", (name or "").strip().lower())


def _new_merge_record() -> dict:
    return {
        "name": None, "sources": [], "phase": 0, "drug_type": "",
        "chembl_id": "", "route": [], "mechanisms": [], "indications": [],
        "is_withdrawn": False, "withdrawal_reasons": [], "black_box_warning": False,
        "trial_stopped": False, "trial_stop_reasons": [], "trial_stop_example": "",
        "_approved_flag": False,
    }


def _merge_drug_records(
    chembl_records: list[dict],
    ot_records: list[dict],
    dgidb_records: list[dict],
    pharos_records: list[dict],
) -> list[dict]:
    """Dedup by normalized name across all four sources into one structured,
    uncapped list, each carrying a derived `status`. See module docstring
    for the exact shape."""
    merged: dict[str, dict] = {}

    def _get(name: str) -> dict:
        key = _norm_name(name)
        return merged.setdefault(key, _new_merge_record())

    for r in chembl_records:
        rec = _get(r["name"])
        rec["name"] = rec["name"] or r["name"]
        if "chembl" not in rec["sources"]:
            rec["sources"].append("chembl")
        rec["phase"] = max(rec["phase"], r["phase"])
        rec["chembl_id"] = rec["chembl_id"] or r.get("chembl_id", "")
        rec["drug_type"] = rec["drug_type"] or r.get("molecule_type", "")
        rec["route"] = rec["route"] or r.get("route", [])
        if r.get("mechanism_of_action") or r.get("action_type"):
            rec["mechanisms"].append({
                "mechanism_of_action": r.get("mechanism_of_action", ""),
                "action_type": r.get("action_type", ""),
            })
        if r.get("withdrawn"):
            rec["is_withdrawn"] = True
        if r.get("black_box_warning"):
            rec["black_box_warning"] = True

    for r in ot_records:
        rec = _get(r["name"])
        rec["name"] = rec["name"] or r["name"]
        if "ot" not in rec["sources"]:
            rec["sources"].append("ot")
        rec["phase"] = max(rec["phase"], r["phase"])
        rec["drug_type"] = rec["drug_type"] or r.get("drug_type", "")
        rec["mechanisms"].extend(r.get("mechanisms", []))
        rec["indications"].extend(r.get("indications", []))
        if r.get("is_withdrawn"):
            rec["is_withdrawn"] = True
            for w in r.get("warnings", []):
                if w.get("warning_type") == "Withdrawn" and w.get("description"):
                    label = w.get("efo_term") or w.get("toxicity_class") or ""
                    reason = f"{label}: {w['description']}" if label else w["description"]
                    rec["withdrawal_reasons"].append(reason)
        if r.get("black_box_warning"):
            rec["black_box_warning"] = True
        if r.get("trial_stopped"):
            rec["trial_stopped"] = True
            rec["trial_stop_reasons"].extend(r.get("trial_stop_reasons", []))
            rec["trial_stop_example"] = rec["trial_stop_example"] or r.get("trial_stop_example", "")

    for r in dgidb_records:
        rec = _get(r["name"])
        rec["name"] = rec["name"] or r["name"]
        if "dgidb" not in rec["sources"]:
            rec["sources"].append("dgidb")
        if r.get("approved"):
            rec["_approved_flag"] = True
        existing_types = {m.get("action_type", "").lower() for m in rec["mechanisms"]}
        for t in r.get("interaction_types", []):
            if t and t.lower() not in existing_types:
                rec["mechanisms"].append({"mechanism_of_action": "", "action_type": t})
                existing_types.add(t.lower())

    for r in pharos_records:
        rec = _get(r["name"])
        rec["name"] = rec["name"] or r["name"]
        if "pharos" not in rec["sources"]:
            rec["sources"].append("pharos")

    out: list[dict] = []
    for rec in merged.values():
        by_disease: dict[str, int] = {}
        for ind in rec["indications"]:
            d = ind["disease"]
            by_disease[d] = max(by_disease.get(d, 0), ind["phase"])
        indications_sorted = sorted(by_disease.items(), key=lambda kv: -kv[1])
        rec["indications"] = [{"disease": d, "phase": p} for d, p in indications_sorted]
        rec["approved_indications"] = [d for d, p in indications_sorted if p >= 4]

        seen_mech = set()
        mechs = []
        for m in rec["mechanisms"]:
            mk = (m.get("mechanism_of_action", ""), m.get("action_type", ""))
            if mk == ("", "") or mk in seen_mech:
                continue
            seen_mech.add(mk)
            mechs.append(m)
        rec["mechanisms"] = mechs
        rec["withdrawal_reasons"] = list(dict.fromkeys(rec["withdrawal_reasons"]))
        rec["trial_stop_reasons"] = list(dict.fromkeys(rec["trial_stop_reasons"]))

        is_approved = rec["phase"] >= 4 or rec["_approved_flag"]
        if rec["is_withdrawn"]:
            status = STATUS_WITHDRAWN
        elif is_approved:
            status = STATUS_APPROVED
        elif rec["trial_stopped"]:
            status = STATUS_FAILED_TRIAL
        elif rec["phase"] >= 1:
            status = STATUS_IN_TRIALS
        else:
            status = STATUS_INVESTIGATIONAL
        rec["status"] = status
        del rec["_approved_flag"]
        out.append(rec)

    out.sort(key=lambda r: (-r["phase"], (r["name"] or "").lower()))
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(df: pd.DataFrame) -> pd.DataFrame:
    """Enrich df with ChEMBL + DGIdb + Open Targets + Pharos drug-target
    columns, plus the unified drug_records_json column."""

    # Deduplicate by gene for API efficiency
    gene_df = (df[["gene_name", "uniprot_acc"]]
               .drop_duplicates("gene_name")
               .to_dict("records"))
    unique_genes = [r["gene_name"] for r in gene_df]
    _log(f"{len(unique_genes):,} unique genes to query")

    # ── ChEMBL ──────────────────────────────────────────────────────────────
    # _v2 cache: schema expanded (full molecule detail + mechanism info per
    # drug, not just name+phase) -- old j3_chembl.json cache lacks these
    # fields, so it must not be reused.
    chembl_cache = CACHE_DIR / "j3_chembl_v2.json"
    if chembl_cache.exists():
        _log("ChEMBL: loading from cache")
        chembl_map: dict[str, dict] = json.loads(chembl_cache.read_text())
    else:
        _log(f"querying ChEMBL ({_CHEMBL_WORKERS} workers) …")
        chembl_map = {}
        done = 0
        total = len(gene_df)

        with ThreadPoolExecutor(max_workers=_CHEMBL_WORKERS) as pool:
            futures = {pool.submit(_fetch_one_gene, row): row["gene_name"]
                       for row in gene_df}
            for fut in as_completed(futures):
                gene, result = fut.result()
                chembl_map[gene] = result
                done += 1
                if done % 50 == 0 or done == total:
                    _log(f"  ChEMBL {done:,}/{total:,} genes")

        chembl_cache.parent.mkdir(parents=True, exist_ok=True)
        chembl_cache.write_text(json.dumps(chembl_map))

    n_druggable = sum(1 for v in chembl_map.values() if v.get("is_druggable"))
    _log(f"ChEMBL: {n_druggable:,}/{len(chembl_map):,} genes have a target entry")

    # ── DGIdb ───────────────────────────────────────────────────────────────
    # _v3 cache: schema expanded (approved flag + interactionTypes per drug).
    dgidb_cache = CACHE_DIR / "j3_dgidb_v3.json"
    if dgidb_cache.exists():
        _log("DGIdb: loading from cache")
        dgidb_map: dict[str, dict] = json.loads(dgidb_cache.read_text())
    else:
        _log("querying DGIdb …")
        dgidb_map = {}
        n_batches = math.ceil(len(unique_genes) / DGIDB_BATCH)
        for i in range(n_batches):
            batch = unique_genes[i * DGIDB_BATCH:(i + 1) * DGIDB_BATCH]
            dgidb_map.update(_dgidb_query(batch))
            _log(f"  DGIdb batch {i+1}/{n_batches}")
            time.sleep(0.5)
        dgidb_cache.write_text(json.dumps(dgidb_map))

    n_dgidb = sum(1 for v in dgidb_map.values() if v.get("n", 0) > 0)
    _log(f"DGIdb: {n_dgidb:,} genes have interaction records")

    # ── Open Targets ────────────────────────────────────────────────────────
    gene_to_ensg = (
        df[["gene_name", "ENSG_ID"]]
        .dropna(subset=["ENSG_ID"])
        .drop_duplicates("gene_name")
        .set_index("gene_name")["ENSG_ID"]
        .to_dict()
    )
    unique_ensg = sorted(set(gene_to_ensg.values()))

    # _v3 cache: schema expanded (drugType, mechanisms, indications,
    # drugWarnings per drug row; clinicalReports now read per-row).
    ot_cache = CACHE_DIR / "j3_opentargets_v3.json"
    if ot_cache.exists():
        _log("Open Targets: loading from cache")
        ot_map: dict[str, dict] = json.loads(ot_cache.read_text())
    else:
        _log(f"querying Open Targets (full drug detail + per-drug clinicalReports), "
             f"{len(unique_ensg):,} genes …")
        ot_map = {}
        n_batches = math.ceil(len(unique_ensg) / OT_DRUG_BATCH)
        for i in range(n_batches):
            batch = unique_ensg[i * OT_DRUG_BATCH:(i + 1) * OT_DRUG_BATCH]
            ot_map.update(_ot_query(batch))
            _log(f"  Open Targets batch {i+1}/{n_batches}")
            time.sleep(0.3)
        ot_cache.parent.mkdir(parents=True, exist_ok=True)
        ot_cache.write_text(json.dumps(ot_map))

    n_ot = sum(1 for v in ot_map.values() if v.get("ot_max_phase", 0) > 0)
    n_stopped = sum(1 for v in ot_map.values() if v.get("ot_trials_terminated", 0) > 0)
    _log(f"Open Targets: {n_ot:,}/{len(ot_map):,} genes have a "
         f"drug/clinical candidate, {n_stopped:,} have a terminated/"
         f"withdrawn/suspended trial")

    # ── Pharos ──────────────────────────────────────────────────────────────
    # _v2 cache: schema expanded (named ligands(isdrug:true), not just counts).
    pharos_cache = CACHE_DIR / "j3_pharos_v2.json"
    if pharos_cache.exists():
        _log("Pharos: loading from cache")
        pharos_map: dict[str, dict] = json.loads(pharos_cache.read_text())
    else:
        _log(f"querying Pharos (tdl + ligandCounts + named drug ligands), "
             f"{len(unique_genes):,} genes …")
        pharos_map = {}
        n_batches = math.ceil(len(unique_genes) / PHAROS_BATCH)
        for i in range(n_batches):
            batch = unique_genes[i * PHAROS_BATCH:(i + 1) * PHAROS_BATCH]
            pharos_map.update(_pharos_query(batch))
            _log(f"  Pharos batch {i+1}/{n_batches}")
            time.sleep(0.3)
        pharos_cache.parent.mkdir(parents=True, exist_ok=True)
        pharos_cache.write_text(json.dumps(pharos_map))

    n_tchem_or_better = sum(1 for v in pharos_map.values()
                             if v.get("pharos_tdl") in {"Tclin", "Tchem"})
    _log(f"Pharos: {n_tchem_or_better:,}/{len(pharos_map):,} genes are "
         f"Tclin/Tchem (known chemical ligand)")

    # ── Attach legacy aggregate columns to dataframe ───────────────────────
    for col in ["chembl_target_id", "chembl_max_phase", "n_drugs", "drug_names",
                "is_druggable", "chembl_bioactive_compounds", "chembl_best_pchembl"]:
        _default = {"is_druggable": False, "chembl_max_phase": 0, "n_drugs": 0,
                    "chembl_bioactive_compounds": 0, "chembl_best_pchembl": 0.0}.get(col, "")
        df[col] = df["gene_name"].map(
            lambda g, _col=col, _def=_default: chembl_map.get(g, {}).get(_col, _def)
        )

    df["dgidb_interactions"] = df["gene_name"].map(
        lambda g: dgidb_map.get(g, {}).get("n", 0)
    )
    df["dgidb_drug_names"] = df["gene_name"].map(
        lambda g: "|".join(r["name"] for r in dgidb_map.get(g, {}).get("records", []))
    )

    _ot_numeric = {"ot_max_phase", "ot_n_drugs", "ot_trials_total", "ot_trials_terminated"}
    for col in ["ot_max_phase", "ot_n_drugs", "ot_drug_names",
                "ot_trials_total", "ot_trials_terminated",
                "ot_trial_stop_reasons", "ot_trial_stop_example"]:
        df[col] = df["gene_name"].map(
            lambda g, _col=col: ot_map.get(gene_to_ensg.get(g, ""), {}).get(
                _col, 0 if _col in _ot_numeric else "")
        )

    for col in ["pharos_tdl", "pharos_n_ligands", "pharos_n_drugs"]:
        df[col] = df["gene_name"].map(
            lambda g, _col=col: pharos_map.get(g, {}).get(
                _col, 0 if _col in {"pharos_n_ligands", "pharos_n_drugs"} else "")
        )

    # ── Unified, structured, uncapped drug list ────────────────────────────
    gene_drug_records = {
        gene: _merge_drug_records(
            chembl_map.get(gene, {}).get("chembl_drug_records", []),
            ot_map.get(gene_to_ensg.get(gene, ""), {}).get("ot_drug_records", []),
            dgidb_map.get(gene, {}).get("records", []),
            pharos_map.get(gene, {}).get("pharos_drug_records", []),
        )
        for gene in unique_genes
    }
    df["drug_records_json"] = df["gene_name"].map(
        lambda g: json.dumps(gene_drug_records.get(g, []), separators=(",", ":"))
    )

    n_with_drugs = sum(1 for v in gene_drug_records.values() if v)
    n_total_drugs = sum(len(v) for v in gene_drug_records.values())
    status_counts: dict[str, int] = {}
    for recs in gene_drug_records.values():
        for r in recs:
            status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
    _log(f"drug_records_json: {n_with_drugs:,}/{len(unique_genes):,} genes have "
         f">=1 drug, {n_total_drugs:,} drug records total -- {status_counts}")

    return df

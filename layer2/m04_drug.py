"""
M04 — Drug Target and Known Compounds (plan 5, M04; ChEMBL-first per build decision).

Database: ChEMBL (open, no key). DrugBank is an optional enrichment, gated on
DRUGBANK_API_KEY — not implemented here; ChEMBL covers drug-target status,
AD-direct relevance (via indications), bioactivities, and docking candidates.

Per gene:
  1. Resolve the gene symbol to a human SINGLE PROTEIN ChEMBL target.
  2. Mechanisms -> drugs targeting the gene (action type, mechanism of action).
  3. Drug indications -> AD-direct relevance (Alzheimer / dementia).
  4. Bioactivities (Kd/Ki) -> sub-micromolar docking candidates with SMILES.
  5. RDKit ETKDG -> a 3D conformer (lowest-energy) SDF per docking candidate.

All ChEMBL access is non-blocking (plan 7): on failure the drug section is
marked unavailable and the run continues. drug_target_status is 'known' if any
mechanism/drug is found, else 'novel'.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from layer2.config import CONFIG
from layer2.inputs import GeneCandidate
from layer2.utils.api_client import APIClient, APIError
from layer2.utils.audit import AuditLog

_AD_TERMS = ("alzheimer", "dementia")


@dataclass
class Drug:
    molecule_chembl_id: str
    name: str
    max_phase: float | None
    action_type: str
    mechanism: str
    ad_direct: bool


@dataclass
class DockingCandidate:
    molecule_chembl_id: str
    smiles: str
    assay_type: str          # Kd | Ki
    value_nM: float
    sdf_path: str | None
    conformer_ok: bool


@dataclass
class M04Result:
    gene_id: str
    target_chembl_id: str | None
    drug_target_status: str           # known | novel
    ad_relevance: str                 # direct | indirect | none
    drugs: list[Drug] = field(default_factory=list)
    docking_candidates: list[DockingCandidate] = field(default_factory=list)
    n_direct: int = 0
    n_indirect: int = 0
    available: bool = True
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# ChEMBL helpers
# ---------------------------------------------------------------------------

# Tight budget: ChEMBL is sometimes slow/down — fail fast to 'unavailable'
# rather than hang the whole run. Cached responses are unaffected.
_CHEMBL_TIMEOUT = 15
_CHEMBL_RETRIES = 1


def _chembl(client: APIClient, path: str, params: dict) -> dict:
    url = f"{CONFIG.apis['chembl'].rstrip('/')}/{path}"
    return client.get_json("chembl", url, params={"format": "json", **params},
                           timeout=_CHEMBL_TIMEOUT, max_retries=_CHEMBL_RETRIES)


def _resolve_target(client: APIClient, gene: str) -> str | None:
    data = _chembl(client, "target/search", {"q": gene, "limit": 10})
    for t in data.get("targets", []):
        if (t.get("target_type") == "SINGLE PROTEIN"
                and t.get("organism") == "Homo sapiens"):
            return t.get("target_chembl_id")
    return None


def _to_nM(value, units) -> float | None:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    u = (units or "").lower()
    if u == "nm":
        return v
    if u == "um" or u == "µm":
        return v * 1000.0
    if u == "pm":
        return v / 1000.0
    if u == "m":
        return v * 1e9
    return None


# ---------------------------------------------------------------------------
# RDKit conformer generation
# ---------------------------------------------------------------------------

def _generate_conformer_sdf(smiles: str, out_path: Path, config) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = config.rdkit["random_seed"]
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=config.rdkit["n_conformers"],
                                          params=params)
        if not cids:
            return False
        energies = AllChem.MMFFOptimizeMoleculeConfs(mol)
        best = min(range(len(energies)), key=lambda i: energies[i][1])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(out_path))
        writer.write(mol, confId=cids[best])
        writer.close()
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(gene: GeneCandidate, client: APIClient, audit: AuditLog,
        config=CONFIG) -> M04Result:
    potency_nM = config.thresholds["chembl_potency_uM"] * 1000.0

    try:
        target_id = _resolve_target(client, gene.gene_id)
    except APIError as e:
        audit.warning("chembl", f"{gene.gene_id}: target resolution failed: {e}")
        return M04Result(gene.gene_id, None, "novel", "none", available=False,
                         warnings=[f"ChEMBL unavailable: {e}"])

    if not target_id:
        audit.info(f"M04 {gene.gene_id}: no human ChEMBL target — novel")
        return M04Result(gene.gene_id, None, "novel", "none")

    res = M04Result(gene.gene_id, target_id, "novel", "none")

    # Mechanisms -> drugs
    try:
        mechs = _chembl(client, "mechanism",
                        {"target_chembl_id": target_id, "limit": 100}).get("mechanisms", [])
    except APIError as e:
        audit.warning("chembl", f"{gene.gene_id}: mechanism query failed: {e}")
        res.available = False
        res.warnings.append(f"ChEMBL mechanism unavailable: {e}")
        mechs = []

    for m in mechs:
        mol_id = m.get("molecule_chembl_id")
        if not mol_id:
            continue
        name, max_phase, ad_direct = mol_id, None, False
        try:
            mol = _chembl(client, f"molecule/{mol_id}", {})
            name = mol.get("pref_name") or mol_id
            max_phase = mol.get("max_phase")
            inds = _chembl(client, "drug_indication",
                           {"molecule_chembl_id": mol_id, "limit": 50}
                           ).get("drug_indications", [])
            ad_direct = any(
                any(term in (str(i.get(k) or "").lower())
                    for k in ("mesh_heading", "efo_term"))
                for i in inds for term in _AD_TERMS)
        except APIError:
            pass
        res.drugs.append(Drug(mol_id, name, max_phase, m.get("action_type", ""),
                              m.get("mechanism_of_action", "") or "", ad_direct))

    # Bioactivities -> sub-uM docking candidates
    try:
        acts = _chembl(client, "activity",
                       {"target_chembl_id": target_id,
                        "standard_type__in": "Kd,Ki", "limit": 200}
                       ).get("activities", [])
    except APIError as e:
        audit.warning("chembl", f"{gene.gene_id}: activity query failed: {e}")
        acts = []

    seen_mol: set[str] = set()
    for a in acts:
        nM = _to_nM(a.get("standard_value"), a.get("standard_units"))
        mol_id = a.get("molecule_chembl_id")
        if nM is None or nM >= potency_nM or not mol_id or mol_id in seen_mol:
            continue
        seen_mol.add(mol_id)
        smiles = ""
        try:
            mol = _chembl(client, f"molecule/{mol_id}", {})
            smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles", "") or ""
        except APIError:
            pass
        sdf_path = conf_ok = None
        if smiles:
            out = config.alphafold_dir / "_compounds" / gene.gene_id / f"{mol_id}.sdf"
            conf_ok = _generate_conformer_sdf(smiles, out, config)
            sdf_path = str(out) if conf_ok else None
        res.docking_candidates.append(DockingCandidate(
            mol_id, smiles, a.get("standard_type", ""), round(nM, 2),
            sdf_path, bool(conf_ok)))

    # Summaries
    res.n_direct = sum(1 for d in res.drugs if d.ad_direct)
    res.drug_target_status = "known" if (res.drugs or res.docking_candidates) else "novel"
    if res.n_direct:
        res.ad_relevance = "direct"
    elif res.drugs:
        res.ad_relevance = "indirect"   # has drugs; indirect AD link resolved with M06 in M08
    else:
        res.ad_relevance = "none"

    audit.info(f"M04 {gene.gene_id} ({target_id}): {len(res.drugs)} drugs "
               f"({res.n_direct} AD-direct), {len(res.docking_candidates)} docking candidates, "
               f"status={res.drug_target_status}")
    return res


if __name__ == "__main__":
    from layer2.inputs import load_inputs
    audit = AuditLog(CONFIG.audit_log_path, verbose=False)
    client = APIClient(audit)
    inp = load_inputs()
    print(f"{'gene':9} {'target':16} {'status':6} {'ad_rel':8} {'drugs':>5} "
          f"{'direct':>6} {'dock':>4}")
    print("-" * 70)
    for gene in inp.gene_work_list:
        r = run(gene, client, audit)
        print(f"{r.gene_id:9} {str(r.target_chembl_id):16} {r.drug_target_status:6} "
              f"{r.ad_relevance:8} {len(r.drugs):5} {r.n_direct:6} {len(r.docking_candidates):4}")

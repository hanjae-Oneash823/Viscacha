"""M1 — SMILES + 3D conformers for user-selected drug candidates.

Only fetches/embeds the drugs the user actually picked for a given hit (via
dossier_server's cart), not every candidate name in hits_deep.csv -- keeps
this fast enough to run synchronously inside an export request instead of
needing the background job queue that m2_structures.py needs.

Extends junior_surveyor/j3_drug_targets.py's ChEMBL molecule-lookup pattern
(same REST base, same batched-fetch idea) with canonical_smiles, and follows
its whole-cache-file convention (outputs/junior_surveyor/cache/j3_*.json)
under master_surveyor's own cache dir instead.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem

from master_surveyor.config import CACHE_DIR, CHEMBL_BASE, CHEMBL_BACKOFF, CHEMBL_RETRIES, CHEMBL_TIMEOUT
from master_surveyor.utils.http import HTTPError, get_json

SMILES_CACHE_PATH = CACHE_DIR / "m1_chembl_smiles.json"


def _log(msg: str) -> None:
    print(f"  [m1] {msg}", file=sys.stderr, flush=True)


def _load_cache() -> dict:
    if SMILES_CACHE_PATH.exists():
        return json.loads(SMILES_CACHE_PATH.read_text())
    return {}


def _save_cache(cache: dict) -> None:
    SMILES_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    SMILES_CACHE_PATH.write_text(json.dumps(cache, indent=None))


def _chembl_molecule_by_name(drug_name: str) -> dict | None:
    """Exact (case-insensitive) pref_name lookup -- the tier-2 fallback the
    plan doc calls for when a drug name resolved by J3's mechanism-based
    lookup doesn't carry a SMILES (e.g. Pharos/DGIdb-only evidence).
    """
    url = f"{CHEMBL_BASE}/molecule"
    params = {
        "pref_name__iexact": drug_name,
        "format": "json",
        "fields": "molecule_chembl_id,pref_name,max_phase,molecule_structures",
        "limit": 1,
    }
    # Network/transient failures propagate (HTTPError) rather than returning
    # None here -- the caller must not cache a transient failure the same
    # way as a confirmed "not in ChEMBL", or a later retry would be
    # permanently blocked by its own cache.
    result = get_json(url, params=params, timeout=CHEMBL_TIMEOUT,
                       retries=CHEMBL_RETRIES, backoff=CHEMBL_BACKOFF)

    molecules = result.get("molecules", [])
    if not molecules:
        return None
    mol = molecules[0]
    structures = mol.get("molecule_structures") or {}
    smiles = structures.get("canonical_smiles")
    if not smiles:
        return None
    return {
        "chembl_id": mol.get("molecule_chembl_id", ""),
        "pref_name": mol.get("pref_name") or drug_name,
        "max_phase": int(float(mol.get("max_phase") or 0)),
        "smiles": smiles,
    }


def fetch_smiles(drug_names: list[str]) -> dict[str, dict | None]:
    """Return {drug_name: {chembl_id, pref_name, max_phase, smiles} | None}.

    None means no SMILES could be resolved for that name (not found in
    ChEMBL by exact pref_name, or found with no structure on record) --
    callers should surface this as "ligand unavailable", not fail the whole
    export. A transient network failure also comes back as None for this
    call, but is deliberately NOT written to the cache (only a confirmed
    result -- found-with-SMILES or confirmed-not-found -- is), so the next
    run retries instead of being stuck with a permanent false negative.
    """
    cache = _load_cache()
    out: dict[str, dict | None] = {}
    dirty = False
    for name in drug_names:
        key = name.strip()
        if not key:
            continue
        if key in cache:
            out[key] = cache[key]
            continue
        _log(f"fetching SMILES for {key!r} …")
        try:
            record = _chembl_molecule_by_name(key)
        except HTTPError as e:
            _log(f"{key!r}: ChEMBL lookup failed, not caching ({e})")
            out[key] = None
            continue
        cache[key] = record
        out[key] = record
        dirty = True
    if dirty:
        _save_cache(cache)
    return out


def embed_conformer(smiles: str, seed: int = 0xC0FFEE) -> Chem.Mol | None:
    """RDKit ETKDG embed + MMFF optimize. Returns None (not an exception) on
    any embedding/optimization failure -- callers treat a missing conformer
    the same way as a missing SMILES: skip that one ligand, keep going.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass   # keep the unoptimized (but successfully embedded) conformer
    return mol


def build_ligands_sdf(drug_names: list[str], out_path: Path) -> list[dict]:
    """Write one multi-molecule SDF (one entry per resolvable drug) to
    out_path. Returns a per-drug status list (name, chembl_id, max_phase,
    smiles, ok, reason) for the export manifest / job status -- callers
    should surface failures (no SMILES, embedding failed) to the user
    rather than silently dropping ligands from the count.
    """
    smiles_map = fetch_smiles(drug_names)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    statuses: list[dict] = []

    for name in drug_names:
        record = smiles_map.get(name.strip())
        if record is None:
            statuses.append({"name": name, "ok": False, "reason": "no SMILES resolved via ChEMBL"})
            continue

        mol = embed_conformer(record["smiles"])
        if mol is None:
            statuses.append({
                "name": name, "ok": False, "reason": "RDKit conformer embedding failed",
                "chembl_id": record["chembl_id"], "smiles": record["smiles"],
            })
            continue

        mol.SetProp("_Name", record["pref_name"])
        mol.SetProp("chembl_id", record["chembl_id"])
        mol.SetProp("max_phase", str(record["max_phase"]))
        mol.SetProp("smiles", record["smiles"])
        writer.write(mol)
        statuses.append({
            "name": name, "ok": True, "chembl_id": record["chembl_id"],
            "max_phase": record["max_phase"], "smiles": record["smiles"],
        })

    writer.close()
    n_ok = sum(1 for s in statuses if s["ok"])
    _log(f"wrote {n_ok}/{len(drug_names)} ligands -> {out_path}")
    return statuses


if __name__ == "__main__":
    test_drugs = ["Donepezil", "Memantine", "Galantamine"]
    result = fetch_smiles(test_drugs)
    for name, rec in result.items():
        print(name, "->", rec)

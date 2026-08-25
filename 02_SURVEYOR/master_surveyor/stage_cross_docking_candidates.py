#!/usr/bin/env python3
"""Stage ligand conformers and canonical receptors for cross-docking rows."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"
SYSTEMS = CAMPAIGN / "systems"

DRUGS = {
    "GABRA2_AZD7325": (
        "AZD7325",
        "CCCNC(=O)C1=NN=C2C(=C1N)C=CC=C2C3=C(C=CC=C3F)OC",
    ),
    "CACNA1D_isradipine": (
        "isradipine",
        "CC1=C(C(C(=C(N1)C)C(=O)OC(C)C)C2=CC=CC3=NON=C32)C(=O)OC",
    ),
    "PDE9A_BI409306": (
        "BI409306",
        "C1COCCC1N2C3=C(C=N2)C(=O)NC(=N3)CC4=CC=CC=N4",
    ),
}


class ProteinChains(Select):
    def __init__(self, chains: set[str]) -> None:
        self.chains = chains

    def accept_chain(self, chain) -> bool:
        return chain.id in self.chains

    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=False))


def generate_ligand(candidate: str, name: str, smiles: str) -> dict[str, object]:
    prepared = SYSTEMS / candidate / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260825
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"RDKit embedding failed for {name}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        force_field = "MMFF94"
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        force_field = "UFF"
    mol.SetProp("_Name", name)
    mol.SetProp("SMILES", smiles)
    sdf = prepared / f"{name}.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    return {
        "candidate": candidate,
        "drug": name,
        "smiles": smiles,
        "conformer_seed": 20260825,
        "optimization": force_field,
        "sdf": str(sdf.relative_to(ROOT)),
    }


def extract_protein(structure, chains: set[str], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(destination), ProteinChains(chains))


def main() -> None:
    records = [generate_ligand(candidate, *drug) for candidate, drug in DRUGS.items()]

    gabra_base = SYSTEMS / "GABRA2_AZD7325"
    gabra = MMCIFParser(QUIET=True).get_structure("9CSB", gabra_base / "inputs" / "9CSB.cif")
    gabra_receptor = gabra_base / "prepared" / "gabra_9CSB_native_pentamer.pdb"
    extract_protein(gabra, {"A", "B", "C", "D", "E"}, gabra_receptor)
    records.append({
        "candidate": "GABRA2_AZD7325",
        "canonical_template": "9CSB",
        "assembly": "beta3-alpha1-beta2-alpha2-gamma2",
        "alpha2_chain": "D",
        "gamma2_chain": "E",
        "receptor_pdb": str(gabra_receptor.relative_to(ROOT)),
    })

    pde_base = SYSTEMS / "PDE9A_BI409306"
    pde = PDBParser(QUIET=True).get_structure("4GH6", pde_base / "inputs" / "4GH6.pdb")
    pde_receptor = pde_base / "prepared" / "pde9a_4GH6_chain_A_protein.pdb"
    extract_protein(pde, {"A"}, pde_receptor)
    luo = pde[0]["A"][("H_LUO", 601, " ")]
    coords = np.asarray([a.coord for a in luo if (a.element or "").upper() != "H"])
    records.append({
        "candidate": "PDE9A_BI409306",
        "canonical_template": "4GH6",
        "protein_chain": "A",
        "template_ligand": "LUO A601",
        "center_A": [round(float(x), 3) for x in coords.mean(axis=0)],
        "box_size_A": [round(float(x + 12.0), 3) for x in np.ptp(coords, axis=0)],
        "receptor_pdb": str(pde_receptor.relative_to(ROOT)),
    })

    output = CAMPAIGN / "analysis" / "metadata" / "cross_docking_stage_metadata.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()

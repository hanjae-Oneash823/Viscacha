#!/usr/bin/env python3
"""Stage matched canonical and KIT-223 kinase-pocket receptors for docking.

Both sequence-derived models are structurally aligned to the experimental
c-KIT--imatinib complex (1T46).  The imatinib coordinates define the common
ATP-site grid; masitinib is then prepared independently for docking.  This is
an exploratory comparative workflow, not a masitinib co-crystal validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select, Superimposer
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
from rdkit import Chem
from rdkit.Chem import AllChem


MASITINIB_SMILES = "Cc1ccc(cc1Nc2nc(cs2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"


def residues(chain):
    return [r for r in chain if is_aa(r, standard=False) and "CA" in r]


class NearbyProtein(Select):
    def __init__(self, center: np.ndarray, radius: float) -> None:
        self.center, self.radius = center, radius

    def accept_residue(self, residue) -> bool:
        if not is_aa(residue, standard=False):
            return False
        return any(np.linalg.norm(atom.coord - self.center) <= self.radius for atom in residue.get_atoms())


def ligand_coords(chain, resname: str, resseq: int) -> np.ndarray:
    residue = chain[("H_" + resname, resseq, " ")]
    return np.asarray([a.coord for a in residue.get_atoms() if a.element != "H"], dtype=float)


def align_and_write(template_chain, model_path: Path, center: np.ndarray, out: Path, deletion_at: int | None = None) -> dict:
    structure = PDBParser(QUIET=True).get_structure(model_path.stem, model_path)
    model_chain = next(structure.get_chains())
    template_res, model_res = residues(template_chain), residues(model_chain)
    model_by_resnum = {r.id[1]: r for r in model_res}
    pairs = []
    for residue in template_res:
        canonical_resnum = residue.id[1]
        model_resnum = canonical_resnum if deletion_at is None or canonical_resnum < deletion_at else canonical_resnum - 1
        if model_resnum in model_by_resnum:
            pairs.append((residue, model_by_resnum[model_resnum]))
    fixed = [template_residue["CA"] for template_residue, _ in pairs]
    moving = [model_residue["CA"] for _, model_residue in pairs]
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    superimposer.apply(list(structure.get_atoms()))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out), NearbyProtein(center, radius=20.0))
    plddt = np.asarray([a.bfactor for a in structure.get_atoms()])
    return {"matched_ca": len(pairs), "alignment_rmsd_A": round(float(superimposer.rms), 3), "mean_model_plddt": round(float(plddt.mean()), 2)}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--template", type=Path, required=True)
    p.add_argument("--canonical-model", type=Path, required=True)
    p.add_argument("--alternate-model", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    template = PDBParser(QUIET=True).get_structure("1T46", args.template)
    chain = template[0]["A"]
    center = ligand_coords(chain, "STI", 3).mean(axis=0)
    extent = np.ptp(ligand_coords(chain, "STI", 3), axis=0)

    canonical = align_and_write(chain, args.canonical_model, center, args.outdir / "kit_canonical_kinase_aligned.pdb")
    alternate = align_and_write(chain, args.alternate_model, center, args.outdir / "kit_223_kinase_aligned.pdb", deletion_at=715)

    mol = Chem.AddHs(Chem.MolFromSmiles(MASITINIB_SMILES))
    params = AllChem.ETKDGv3(); params.randomSeed = 20260824
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("RDKit could not generate a masitinib conformer")
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("_Name", "masitinib")
    writer = Chem.SDWriter(str(args.outdir / "masitinib.sdf")); writer.write(mol); writer.close()

    metadata = {
        "template_pdb": "1T46", "template_ligand": "STI (imatinib)",
        "masitinib_smiles": MASITINIB_SMILES,
        "common_grid_center_A": [round(float(x), 3) for x in center],
        "common_grid_size_A": [round(float(x + 10), 3) for x in extent],
        "canonical": canonical, "kit_223": alternate,
        "interpretation": "KIT-223 differs by one residue at canonical position 715. Comparative docking is exploratory and requires model-quality/pose inspection; it does not alone demonstrate an affinity difference.",
    }
    (args.outdir / "stage_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

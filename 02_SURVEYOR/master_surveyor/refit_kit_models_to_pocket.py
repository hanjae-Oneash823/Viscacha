#!/usr/bin/env python3
"""Locally refit existing KIT isoform models to the experimental 1T46 pocket.

The original whole-domain superposition left the predicted ligand pockets
outside the experimental grid.  Local superposition is the appropriate
coordinate transfer for matched cross-docking and is quality-checked here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.PDB import PDBIO, PDBParser, Select, Superimposer
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "outputs" / "docking_campaign" / "KIT_masitinib"


def residues_seq(chain):
    residues = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    seq = "".join(protein_letters_3to1.get(r.resname, "X") for r in residues)
    return residues, seq


def pairs(fixed_chain, moving_chain):
    fixed, a = residues_seq(fixed_chain)
    moving, b = residues_seq(moving_chain)
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"; aligner.match_score = 2; aligner.mismatch_score = -1
    aligner.open_gap_score = -8; aligner.extend_gap_score = -0.5
    alignment = aligner.align(a, b)[0]
    output = []
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        for offset in range(min(a1 - a0, b1 - b0)):
            output.append((fixed[a0 + offset], moving[b0 + offset]))
    return output


class NearPocket(Select):
    def __init__(self, center: np.ndarray, radius: float = 22.0) -> None:
        self.center, self.radius = center, radius

    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=False)) and any(
            np.linalg.norm(atom.coord - self.center) <= self.radius for atom in residue
        )


def refit(label: str, input_name: str, template_chain, center: np.ndarray) -> dict[str, object]:
    model = PDBParser(QUIET=True).get_structure(label, BASE / "prepared" / input_name)
    chain = next(model.get_chains())
    all_pairs = pairs(template_chain, chain)
    local_pairs = [
        (a, b) for a, b in all_pairs
        if any(np.linalg.norm(atom.coord - center) <= 12.0 for atom in a)
    ]
    sup = Superimposer()
    sup.set_atoms([a["CA"] for a, _ in local_pairs], [b["CA"] for _, b in local_pairs])
    sup.apply(list(model.get_atoms()))
    out = BASE / "prepared" / f"kit_{label}_pocket_refit.pdb"
    io = PDBIO(); io.set_structure(model); io.save(str(out), NearPocket(center))
    plddt = [float(b["CA"].bfactor) for _, b in local_pairs]
    return {
        "label": label, "input": input_name, "output": str(out.relative_to(ROOT)),
        "matched_local_CA": len(local_pairs), "local_CA_RMSD_A": round(float(sup.rms), 3),
        "mean_local_plddt": round(float(np.mean(plddt)), 2), "minimum_local_plddt": round(float(np.min(plddt)), 2),
    }


def main() -> None:
    template = PDBParser(QUIET=True).get_structure("1T46", BASE / "inputs" / "1T46.pdb")
    ligand = template[0]["A"][("H_STI", 3, " ")]
    center = np.mean([a.coord for a in ligand if (a.element or "").upper() != "H"], axis=0)
    records = [
        refit("canonical", "kit_canonical_kinase_aligned.pdb", template[0]["A"], center),
        refit("223", "kit_223_kinase_aligned.pdb", template[0]["A"], center),
    ]
    (BASE / "prepared" / "pocket_refit_quality.json").write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()

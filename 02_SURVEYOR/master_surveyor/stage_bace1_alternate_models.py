#!/usr/bin/env python3
"""Align AlphaFold BACE1 deletion-isoform models to the 5HU1 pocket frame."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "outputs" / "docking_campaign" / "BACE1_verubecestat"


def residues_and_sequence(chain):
    residues = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    sequence = "".join(protein_letters_3to1.get(r.resname, "X") for r in residues)
    return residues, sequence


def pair_residues(fixed_chain, moving_chain):
    fixed, fixed_seq = residues_and_sequence(fixed_chain)
    moving, moving_seq = residues_and_sequence(moving_chain)
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(fixed_seq, moving_seq)[0]
    pairs = []
    for (f0, f1), (m0, m1) in zip(*alignment.aligned):
        for offset in range(min(f1 - f0, m1 - m0)):
            pairs.append((fixed[f0 + offset], moving[m0 + offset]))
    return pairs


def stage(isoform: str) -> dict[str, object]:
    experimental = PDBParser(QUIET=True).get_structure("5HU1", BASE / "inputs" / "5HU1.pdb")
    crystal_chain = experimental[0]["A"]
    model_path = next((BASE / "models" / f"BACE1_{isoform}").glob("*_unrelaxed_rank_001_*.pdb"))
    model = PDBParser(QUIET=True).get_structure(f"BACE1_{isoform}", model_path)
    model_chain = next(model.get_chains())

    pairs = pair_residues(crystal_chain, model_chain)
    superimposer = Superimposer()
    superimposer.set_atoms([a["CA"] for a, _ in pairs], [b["CA"] for _, b in pairs])
    superimposer.apply(list(model.get_atoms()))

    ligand = crystal_chain[("H_66F", 501, " ")]
    ligand_coords = np.asarray([a.coord for a in ligand if (a.element or "").upper() != "H"])
    contact_numbers = set()
    for residue in crystal_chain:
        if not is_aa(residue, standard=False):
            continue
        if any(np.min(np.linalg.norm(ligand_coords - atom.coord, axis=1)) <= 6.0 for atom in residue):
            contact_numbers.add(residue.id[1])

    pocket_pairs = [(a, b) for a, b in pairs if a.id[1] in contact_numbers]
    pocket_rmsd = float(np.sqrt(np.mean([
        np.sum((a["CA"].coord - b["CA"].coord) ** 2) for a, b in pocket_pairs
    ])))
    pocket_plddt = [float(b["CA"].bfactor) for _, b in pocket_pairs]
    all_plddt = [float(r["CA"].bfactor) for r in model_chain if "CA" in r]

    destination = BASE / "prepared" / f"bace1_{isoform}_alphafold_aligned_to_5HU1.pdb"
    io = PDBIO()
    io.set_structure(model)
    io.save(str(destination))
    return {
        "isoform": f"BACE1-{isoform}",
        "source_model": str(model_path.relative_to(ROOT)),
        "aligned_model": str(destination.relative_to(ROOT)),
        "matched_CA": len(pairs),
        "global_fit_CA_RMSD_A": round(float(superimposer.rms), 3),
        "pocket_contact_residues_in_5HU1": sorted(contact_numbers),
        "mapped_pocket_CA_count": len(pocket_pairs),
        "pocket_CA_RMSD_after_global_fit_A": round(pocket_rmsd, 3),
        "mean_model_plddt": round(float(np.mean(all_plddt)), 2),
        "mean_mapped_pocket_plddt": round(float(np.mean(pocket_plddt)), 2),
        "minimum_mapped_pocket_plddt": round(float(np.min(pocket_plddt)), 2),
    }


def main() -> None:
    records = [stage("476"), stage("457")]
    output = BASE / "prepared" / "alternate_model_quality.json"
    output.write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()

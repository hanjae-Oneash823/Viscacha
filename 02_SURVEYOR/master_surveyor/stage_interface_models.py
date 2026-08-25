#!/usr/bin/env python3
"""Transfer interface pockets and build explicit mixed-receptor hypotheses.

GABRA2: transfer the diazepam-defined alpha/gamma benzodiazepine site from
6X3X onto the native alpha2/gamma2-containing human pentamer 9CSB.

CHRNA7: replace either face of the 7EKP A/B encenicline site with one
experimentally observed CHRFAM7A-fusion extracellular domain from 9QTO.  The
two hybrids are topology hypotheses, not claimed biological stoichiometries.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, Superimposer
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"


def residues_and_sequence(chain):
    residues = [r for r in chain if is_aa(r, standard=False) and "CA" in r]
    sequence = "".join(protein_letters_3to1.get(r.resname, "X") for r in residues)
    return residues, sequence


def matched_residues(fixed_chain, moving_chain):
    fixed_res, fixed_seq = residues_and_sequence(fixed_chain)
    moving_res, moving_seq = residues_and_sequence(moving_chain)
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(fixed_seq, moving_seq)[0]
    fixed_pairs, moving_pairs = [], []
    for (f0, f1), (m0, m1) in zip(*alignment.aligned):
        length = min(f1 - f0, m1 - m0)
        for offset in range(length):
            fixed_pairs.append(fixed_res[f0 + offset])
            moving_pairs.append(moving_res[m0 + offset])
    return fixed_pairs, moving_pairs


def align_chain(fixed_chain, moving_chain, atoms_to_transform) -> dict[str, object]:
    fixed_res, moving_res = matched_residues(fixed_chain, moving_chain)
    superimposer = Superimposer()
    superimposer.set_atoms([r["CA"] for r in fixed_res], [r["CA"] for r in moving_res])
    superimposer.apply(list(atoms_to_transform))
    identities = sum(a.resname == b.resname for a, b in zip(fixed_res, moving_res))
    return {
        "matched_CA": len(fixed_res),
        "identical_residues": identities,
        "sequence_identity_over_match": round(identities / len(fixed_res), 4),
        "fit_CA_RMSD_A": round(float(superimposer.rms), 3),
    }


def write_residue(residue, path: Path, chain_id: str = "T") -> np.ndarray:
    lines, coords = [], []
    for serial, atom in enumerate(residue.get_atoms(), 1):
        element = (atom.element or atom.name[0]).strip().upper()
        x, y, z = atom.coord
        lines.append(
            f"HETATM{serial:5d} {atom.name:>4s} {residue.resname:>3s} {chain_id}{residue.id[1]:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{atom.bfactor:6.2f}          {element:>2s}\n"
        )
        if element != "H":
            coords.append([float(x), float(y), float(z)])
    lines.append("END\n")
    path.write_text("".join(lines))
    return np.asarray(coords)


class ProteinOnly(Select):
    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=False))


def build_hybrid(target, fusion, replace_chain: str, destination: Path) -> dict[str, object]:
    moving = copy.deepcopy(fusion)
    fusion_chain = moving[0]["A"]
    metrics = align_chain(target[0][replace_chain], fusion_chain, moving.get_atoms())
    fusion_copy = copy.deepcopy(fusion_chain)
    fusion_copy.id = replace_chain

    hybrid = Structure(destination.stem)
    model = Model(0)
    hybrid.add(model)
    for chain in target[0]:
        if chain.id in {"A", "B", "C", "D", "E"} and chain.id != replace_chain:
            model.add(copy.deepcopy(chain))
    model.add(fusion_copy)

    io = PDBIO()
    io.set_structure(hybrid)
    io.save(str(destination), ProteinOnly())
    metrics.update({
        "replaced_7EKP_chain": replace_chain,
        "fusion_source": "9QTO chain A extracellular-domain construct",
        "output": str(destination.relative_to(ROOT)),
    })
    return metrics


def main() -> None:
    records: dict[str, object] = {}

    # GABRA2 alpha2/gamma2 interface grid transfer.
    gabra_base = CAMPAIGN / "GABRA2_AZD7325"
    target = MMCIFParser(QUIET=True).get_structure("9CSB", gabra_base / "inputs" / "9CSB.cif")
    mobile = PDBParser(QUIET=True).get_structure("6X3X", gabra_base / "inputs" / "6X3X.pdb")
    alpha_metrics = align_chain(target[0]["D"], mobile[0]["D"], mobile.get_atoms())
    dzp = mobile[0]["D"][("H_DZP", 404, " ")]
    dzp_path = gabra_base / "prepared" / "diazepam_D404_transferred_to_9CSB.pdb"
    coords = write_residue(dzp, dzp_path)

    # Quantify whether the neighboring gamma2 subunit also overlays after the alpha fit.
    fixed_gamma, moving_gamma = matched_residues(target[0]["E"], mobile[0]["E"])
    gamma_rmsd = float(np.sqrt(np.mean([
        np.sum((a["CA"].coord - b["CA"].coord) ** 2)
        for a, b in zip(fixed_gamma, moving_gamma)
    ])))
    records["GABRA2_AZD7325"] = {
        "template": "6X3X diazepam D404 at alpha1(D)/gamma2(E) interface",
        "target": "9CSB alpha2(D)/gamma2(E) interface",
        "alpha_alignment": alpha_metrics,
        "neighbor_gamma2_CA_RMSD_A_after_alpha_fit": round(gamma_rmsd, 3),
        "transferred_ligand": str(dzp_path.relative_to(ROOT)),
        "center_A": [round(float(x), 3) for x in coords.mean(axis=0)],
        "box_size_A": [round(float(x + 12.0), 3) for x in np.ptp(coords, axis=0)],
    }

    # CHRNA7/CHRFAM7A one-fusion-subunit hypotheses at each face of site A/B.
    chr_base = CAMPAIGN / "CHRNA7_encenicline"
    chr_target = PDBParser(QUIET=True).get_structure("7EKP", chr_base / "inputs" / "7EKP.pdb")
    fusion = MMCIFParser(QUIET=True).get_structure("9QTO", chr_base / "inputs" / "9QTO.cif")
    hybrids = []
    for chain_id, label in [("A", "fusion_at_A_face"), ("B", "fusion_at_B_face")]:
        destination = chr_base / "prepared" / f"chrna7_chrFam7a_{label}.pdb"
        hybrids.append(build_hybrid(chr_target, fusion, chain_id, destination))
    records["CHRNA7_encenicline"] = {
        "site": "7EKP encenicline I33 A601, interface of chains A and B",
        "hybrids": hybrids,
        "interpretation": "Two one-fusion-subunit topology hypotheses bracket which face of the A/B site is altered. They do not establish in-sample stoichiometry.",
    }

    output = CAMPAIGN / "interface_model_metadata.json"
    output.write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()

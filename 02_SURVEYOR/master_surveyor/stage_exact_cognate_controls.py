#!/usr/bin/env python3
"""Extract receptors and cognate ligands for exact redocking controls.

Controls:
  * human BACE1--verubecestat, PDB 5HU1, ligand 66F
  * human alpha-7 nAChR--encenicline, PDB 7EKP, ligand I33

The ligand coordinates remain in the experimental receptor frame so the
resulting docked poses can be evaluated by fixed-frame heavy-atom RMSD.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "outputs" / "docking_campaign"


class ProteinChains(Select):
    def __init__(self, chain_ids: set[str]) -> None:
        self.chain_ids = chain_ids

    def accept_chain(self, chain) -> bool:
        return chain.id in self.chain_ids

    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=False))


def write_ligand(residue, path: Path) -> np.ndarray:
    lines: list[str] = []
    coords: list[list[float]] = []
    atom_serials: list[int] = []
    for serial, atom in enumerate(residue.get_atoms(), 1):
        element = (atom.element or atom.name[0]).strip().upper()
        x, y, z = atom.coord
        lines.append(
            f"HETATM{serial:5d} {atom.name:>4s} {residue.resname:>3s} A{residue.id[1]:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{atom.bfactor:6.2f}          {element:>2s}\n"
        )
        atom_serials.append(serial)
        if element != "H":
            coords.append([float(x), float(y), float(z)])
    lines.append("END\n")
    path.write_text("".join(lines))
    return np.asarray(coords, dtype=float)


def stage(
    candidate: str,
    pdb_id: str,
    protein_chains: set[str],
    ligand_chain: str,
    ligand_resname: str,
    ligand_resseq: int,
    receptor_name: str,
    ligand_name: str,
) -> dict[str, object]:
    base = CAMPAIGN / candidate
    prepared = base / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    structure = PDBParser(QUIET=True).get_structure(pdb_id, base / "inputs" / f"{pdb_id}.pdb")

    receptor_path = prepared / f"{receptor_name}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_path), ProteinChains(protein_chains))

    residue = structure[0][ligand_chain][("H_" + ligand_resname, ligand_resseq, " ")]
    ligand_path = prepared / f"{ligand_name}_crystal.pdb"
    coords = write_ligand(residue, ligand_path)
    extent = np.ptp(coords, axis=0)
    return {
        "candidate": candidate,
        "pdb_id": pdb_id,
        "protein_chains": sorted(protein_chains),
        "experimental_ligand": {
            "chain": ligand_chain,
            "resname": ligand_resname,
            "resseq": ligand_resseq,
            "heavy_atoms": len(coords),
        },
        "center_A": [round(float(x), 3) for x in coords.mean(axis=0)],
        "box_size_A": [round(float(x + 12.0), 3) for x in extent],
        "receptor_pdb": str(receptor_path.relative_to(ROOT)),
        "ligand_crystal_pdb": str(ligand_path.relative_to(ROOT)),
    }


def main() -> None:
    records = [
        stage(
            "BACE1_verubecestat", "5HU1", {"A"}, "A", "66F", 501,
            "bace1_5HU1_chain_A_protein", "verubecestat_66F_A501",
        ),
        stage(
            "CHRNA7_encenicline", "7EKP", {"A", "B", "C", "D", "E"}, "A", "I33", 601,
            "chrna7_7EKP_pentamer_protein", "encenicline_I33_A601",
        ),
    ]
    for record in records:
        base = CAMPAIGN / str(record["candidate"]) / "prepared"
        (base / "stage_metadata.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Stage a transparent FYN--saracatinib redocking control from PDB 10DJ.

The script intentionally extracts only experimental coordinates.  It does not
model the truncated alternate FYN ORF: that ORF ends before the kinase domain,
so an ATP-pocket docking calculation would not be interpretable.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa


class ProteinChainSelect(Select):
    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=False))


def write_ligand_pdb(structure, chain_id: str, resname: str, resseq: int, out: Path) -> np.ndarray:
    residue = structure[0][chain_id][("H_" + resname, resseq, " ")]
    lines: list[str] = []
    heavy_coords: list[list[float]] = []
    for serial, atom in enumerate(residue.get_atoms(), start=1):
        element = (atom.element or atom.get_name()[0]).strip().upper()
        x, y, z = atom.coord
        lines.append(
            f"HETATM{serial:5d} {atom.get_name():>4s} {resname:>3s} {chain_id:1s}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{atom.bfactor:6.2f}          {element:>2s}\n"
        )
        if element != "H":
            heavy_coords.append([float(x), float(y), float(z)])
    lines.append("END\n")
    out.write_text("".join(lines))
    return np.asarray(heavy_coords)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    structure = PDBParser(QUIET=True).get_structure("10DJ", args.pdb)

    receptor_path = args.outdir / "fyn_chain_A_protein.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(receptor_path), ProteinChainSelect("A"))

    ligand_path = args.outdir / "saracatinib_H8H_A601_crystal.pdb"
    coords = write_ligand_pdb(structure, "A", "H8H", 601, ligand_path)
    minimum, maximum = coords.min(axis=0), coords.max(axis=0)
    center = coords.mean(axis=0)
    extent = maximum - minimum
    metadata = {
        "pdb_id": "10DJ",
        "protein_chain": "A",
        "experimental_ligand": {"resname": "H8H", "chain": "A", "resseq": 601},
        "heavy_atom_count": int(len(coords)),
        "crystal_ligand_center_angstrom": [round(float(x), 3) for x in center],
        "crystal_ligand_extent_angstrom": [round(float(x), 3) for x in extent],
        "recommended_vina_box_size_angstrom": [round(float(x + 10.0), 3) for x in extent],
        "note": "Box is ligand extent plus 10 A total padding; validate and enlarge only if poses are clipped.",
    }
    (args.outdir / "stage_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

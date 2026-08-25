#!/usr/bin/env python3
"""Extract the experimental c-KIT--imatinib control from PDB 1T46."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa


class ProteinA(Select):
    def accept_chain(self, chain) -> bool:
        return chain.id == "A"

    def accept_residue(self, residue) -> bool:
        return bool(is_aa(residue, standard=False))


def main() -> None:
    root = Path("outputs/docking_campaign/KIT_masitinib")
    prepared = root / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    structure = PDBParser(QUIET=True).get_structure("1T46", root / "inputs" / "1T46.pdb")
    io = PDBIO(); io.set_structure(structure)
    io.save(str(prepared / "kit_1T46_chain_A_protein.pdb"), ProteinA())

    ligand = structure[0]["A"][("H_STI", 3, " ")]
    heavy: list[list[float]] = []
    lines: list[str] = []
    serial_map: dict[int, int] = {}
    for serial, atom in enumerate(ligand.get_atoms(), 1):
        serial_map[atom.serial_number] = serial
        element = (atom.element or atom.get_name()[0]).strip().upper()
        x, y, z = atom.coord
        lines.append(f"HETATM{serial:5d} {atom.get_name():>4s} STI A   3    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{atom.bfactor:6.2f}          {element:>2s}\n")
        if element != "H":
            heavy.append([float(x), float(y), float(z)])
    # Preserve crystallographic ligand bonds; otherwise PDB-to-SDF conversion
    # may infer incorrect aromatic/valence states for imatinib.
    for line in (root / "inputs" / "1T46.pdb").read_text().splitlines():
        if not line.startswith("CONECT"):
            continue
        fields = [int(x) for x in line[6:].split()]
        if fields and all(atom_serial in serial_map for atom_serial in fields):
            lines.append("CONECT" + "".join(f"{serial_map[atom_serial]:5d}" for atom_serial in fields) + "\n")
    lines.append("END\n")
    (prepared / "imatinib_STI_A3_crystal.pdb").write_text("".join(lines))
    coords = np.asarray(heavy)
    metadata = {
        "pdb_id": "1T46", "experimental_ligand": "STI (imatinib)",
        "center_A": [round(float(x), 3) for x in coords.mean(axis=0)],
        "box_size_A": [round(float(x + 10), 3) for x in np.ptp(coords, axis=0)],
        "heavy_atom_count": int(len(coords)),
    }
    (prepared / "imatinib_redock_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

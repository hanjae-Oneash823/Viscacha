#!/usr/bin/env python3
"""Stage the experimental CaV1.3--amiodarone control (PDB 8E59).

CACNA1D-214 truncates after residue 1625.  The deposited CaV1.3 chain in
8E59 is modelled only through residue 1589, so its ligand pocket is shared by
the canonical and alternate sequences.  This is deliberately a pocket-
preservation control, not an isoform-affinity comparison.
"""

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
    root = Path(__file__).resolve().parents[2]
    campaign = root / "outputs" / "docking_campaign" / "systems" / "CACNA1D_isradipine"
    pdb = campaign / "inputs" / "8E59.pdb"
    prepared = campaign / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    structure = PDBParser(QUIET=True).get_structure("8E59", pdb)
    io = PDBIO(); io.set_structure(structure)
    io.save(str(prepared / "cacna1d_8E59_chain_A_protein.pdb"), ProteinA())

    # In 8E59, BBI is amiodarone. 3PE is a phosphatidylethanolamine lipid and
    # must not be used as the drug reference.
    ligand = structure[0]["A"][("H_BBI", 2201, " ")]
    lines, coords = [], []
    for serial, atom in enumerate(ligand.get_atoms(), 1):
        element = (atom.element or atom.get_name()[0]).strip().upper()
        x, y, z = atom.coord
        lines.append(f"HETATM{serial:5d} {atom.get_name():>4s} BBI A2201    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{atom.bfactor:6.2f}          {element:>2s}\\n")
        if element != "H": coords.append([float(x), float(y), float(z)])
    (prepared / "amiodarone_BBI_A2201_crystal.pdb").write_text("".join(lines) + "END\\n")
    c = np.asarray(coords)
    metadata = {
        "pdb_id": "8E59", "protein_chain": "A", "experimental_ligand": "amiodarone (BBI)",
        "ligand_residue": 2201, "modelled_cacna1d_residue_range": [121, 1589],
        "CACNA1D_214_length": 1625,
        "interpretation": "The entire experimental ligand template lies before the CACNA1D-214 truncation; therefore the template supports a shared-pocket control, not a differential docking claim.",
        "center": [round(float(v), 3) for v in c.mean(axis=0)],
        "box_size": [round(float(v + 12), 3) for v in (c.max(axis=0) - c.min(axis=0))],
    }
    (prepared / "stage_metadata.json").write_text(json.dumps(metadata, indent=2) + "\\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

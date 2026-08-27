#!/usr/bin/env python3
"""Align and quality-check the observed BACE1-202/P56817-5 ensemble.

The observed transcript encodes the 401-aa UniProt isoform 5.  It replaces
canonical residues 1--20 and deletes residues 21--120, so this script measures
how many experimentally observed verubecestat contacts remain before deciding
whether comparative docking is meaningful.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from Bio import Align
from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1


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


def contact_numbers(crystal_chain, cutoff=6.0):
    ligand = crystal_chain[("H_66F", 501, " ")]
    ligand_xyz = np.asarray([a.coord for a in ligand if (a.element or "").upper() != "H"])
    contacts = set()
    for residue in crystal_chain:
        if not is_aa(residue, standard=False):
            continue
        protein_xyz = np.asarray([a.coord for a in residue if (a.element or "").upper() != "H"])
        if len(protein_xyz) and np.linalg.norm(protein_xyz[:, None] - ligand_xyz, axis=2).min() <= cutoff:
            contacts.add(residue.id[1])
    return sorted(contacts)


def score_json_for(model_path):
    return model_path.with_name(model_path.name.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crystal", type=Path, required=True)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    crystal = PDBParser(QUIET=True).get_structure("5HU1", args.crystal)
    crystal_chain = crystal[0]["A"]
    contacts = contact_numbers(crystal_chain)
    rows = []
    model_paths = sorted(args.models.glob("*_unrelaxed_rank_*.pdb"))
    if not model_paths:
        raise FileNotFoundError(f"No ColabFold models found under {args.models}")

    for model_path in model_paths:
        model = PDBParser(QUIET=True).get_structure(model_path.stem, model_path)
        moving_chain = next(model.get_chains())
        pairs = pair_residues(crystal_chain, moving_chain)
        superimposer = Superimposer()
        superimposer.set_atoms([a["CA"] for a, _ in pairs], [b["CA"] for _, b in pairs])
        superimposer.apply(list(model.get_atoms()))
        pair_by_number = {a.id[1]: b for a, b in pairs}
        mapped_contacts = [p for p in contacts if p in pair_by_number]
        missing_contacts = [p for p in contacts if p not in pair_by_number]
        deviations = [
            float(np.linalg.norm(crystal_chain[(" ", p, " ")]["CA"].coord - pair_by_number[p]["CA"].coord))
            for p in mapped_contacts
        ]
        pocket_plddt = [float(pair_by_number[p]["CA"].bfactor) for p in mapped_contacts]
        all_plddt = [float(r["CA"].bfactor) for r in moving_chain if "CA" in r]
        scores_path = score_json_for(model_path)
        scores = json.loads(scores_path.read_text()) if scores_path.exists() else {}
        aligned_path = args.outdir / model_path.name.replace("unrelaxed", "aligned")
        io = PDBIO()
        io.set_structure(model)
        io.save(str(aligned_path))
        rows.append({
            "model": model_path.name,
            "aligned_model": str(aligned_path),
            "matched_ca": len(pairs),
            "global_fit_ca_rmsd_A": float(superimposer.rms),
            "experimental_contacts_total": len(contacts),
            "experimental_contacts_retained": len(mapped_contacts),
            "experimental_contacts_missing": len(missing_contacts),
            "missing_contact_positions": ";".join(map(str, missing_contacts)),
            "retained_contact_ca_displacement_mean_A": float(np.mean(deviations)),
            "retained_contact_ca_displacement_max_A": float(np.max(deviations)),
            "mean_model_plddt": float(np.mean(all_plddt)),
            "mean_retained_contact_plddt": float(np.mean(pocket_plddt)),
            "minimum_retained_contact_plddt": float(np.min(pocket_plddt)),
            "ptm": scores.get("ptm"),
            "catalytic_Asp93_retained": 93 in mapped_contacts,
            "catalytic_Asp289_retained": 289 in mapped_contacts,
        })

    metrics_path = args.outdir / "BACE1_202_401_ensemble_metrics.csv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    top = rows[0]
    summary = {
        "candidate": "BACE1-202 / ENST00000392937 / UniProt P56817-5 (401 aa)",
        "sequence_change": "canonical residues 1-20 replaced; canonical residues 21-120 deleted",
        "dtu": {
            "cell_type": "Oligodendrocyte",
            "psi_AD": 0.185023286028755,
            "psi_control": 0.0853706235285183,
            "delta_psi": 0.0996526625002365,
            "gene_adjusted_p": 0.655428864064341,
            "empirical_FDR": 0.977051941395206,
            "passes_significance_gate": False,
        },
        "experimental_contact_cutoff_A": 6.0,
        "experimental_contact_positions": contacts,
        "ensemble": rows,
        "top_aligned_model": top["aligned_model"],
        "classification": "B structurally; transcript-evidence-limited",
    }
    (args.outdir / "BACE1_202_401_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

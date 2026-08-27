#!/usr/bin/env python3
"""Validate B/C candidate geometry using sequence-aware structural comparisons.

This script performs no docking. It measures whether predicted isoform changes alter
the experimentally defined ligand-binding region, and reports local confidence and
ensemble consistency so model-derived claims can be qualified appropriately.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1


GFRA1_CANONICAL = Path("outputs/master_surveyor/cache/structures/ea3c29477a12212b/canonical_afdb.pdb")
GFRA1_ALT_GLOB = "outputs/master_surveyor/cache/structures/3a0926d056318a13/alt_colabfold/*unrelaxed_rank_*.pdb"
CACNA1D_CANONICAL = Path("outputs/master_surveyor/cache/structures/d96a034227d14cdd/canonical_afdb.pdb")
CACNA1D_ALT_GLOB = "outputs/master_surveyor/cache/structures/a90d14f0a0086f3f/alt_colabfold/*unrelaxed_rank_*.pdb"
SORT1_CANONICAL = Path("outputs/master_surveyor/cache/structures/71c75d29add39385/canonical_afdb.pdb")
SORT1_ALT_GLOB = "outputs/master_surveyor/cache/structures/5c0b1b2b37b5e647/alt_colabfold/*unrelaxed_rank_*.pdb"

# Rat CaV1.3 residues reported around the dihydropyridine pocket. Human mapping is
# checked by sequence alignment; these values are used only as a local-region seed.
CACNA1D_POCKET_CANONICAL = [1078, 1081, 1082, 1085, 1154, 1156, 1194, 1198, 1205, 1209, 1212, 1489, 1492, 1493]

# The latozinemab patent places the discontinuous/conformational epitope within
# human SORT1 207--231 and calls out T218, Y222, S223, and S227 specifically.
# These are patent-reported epitope positions, not contacts inferred by docking.
SORT1_EPITOPE_CANONICAL = list(range(207, 232))
SORT1_EPITOPE_CALLOUTS = [218, 222, 223, 227]


def residues(path: Path, chain_id: str | None = None):
    structure = PDBParser(QUIET=True).get_structure(path.stem, path)
    chains = list(structure[0])
    chain = next((c for c in chains if c.id == chain_id), chains[0]) if chain_id else chains[0]
    return [r for r in chain if r.id[0] == " " and r.resname in protein_letters_3to1]


def sequence(rs):
    return "".join(protein_letters_3to1[r.resname] for r in rs)


def seq_map(ref_rs, query_rs):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5
    aln = aligner.align(sequence(ref_rs), sequence(query_rs))[0]
    mapping = {}
    for (r0, r1), (q0, q1) in zip(aln.aligned[0], aln.aligned[1]):
        for offset in range(min(r1 - r0, q1 - q0)):
            mapping[r0 + offset + 1] = q0 + offset + 1
    return mapping, float(aln.score)


def ca_array(rs, positions):
    return np.asarray([rs[p - 1]["CA"].coord for p in positions], dtype=float)


def kabsch_transform(mobile, target):
    mob_centroid = mobile.mean(axis=0)
    tgt_centroid = target.mean(axis=0)
    h = (mobile - mob_centroid).T @ (target - tgt_centroid)
    u, _, vt = np.linalg.svd(h)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    tran = tgt_centroid - mob_centroid @ rot.T
    moved = mobile @ rot.T + tran
    return rot, tran, float(np.sqrt(np.mean(np.sum((moved - target) ** 2, axis=1))))


def local_metrics(canonical_path, alt_paths, positions, flank=3):
    canon = residues(canonical_path)
    expanded = sorted({p + d for p in positions for d in range(-flank, flank + 1) if 1 <= p + d <= len(canon)})
    rows = []
    for alt_path in alt_paths:
        alt = residues(alt_path)
        mapping, score = seq_map(canon, alt)
        usable = [p for p in expanded if p in mapping and "CA" in canon[p - 1] and "CA" in alt[mapping[p] - 1]]
        c_xyz = ca_array(canon, usable)
        a_xyz = np.asarray([alt[mapping[p] - 1]["CA"].coord for p in usable], dtype=float)
        rot, tran, ca_rmsd = kabsch_transform(a_xyz, c_xyz)
        atom_sq = []
        for p in usable:
            cr = canon[p - 1]
            ar = alt[mapping[p] - 1]
            if cr.resname != ar.resname:
                continue
            shared = sorted(set(a.name for a in cr) & set(a.name for a in ar) - {"H"})
            for name in shared:
                moved = ar[name].coord @ rot.T + tran
                atom_sq.append(float(np.sum((moved - cr[name].coord) ** 2)))
        pocket_alt_positions = [mapping[p] for p in positions if p in mapping]
        plddt = [float(alt[p - 1]["CA"].bfactor) for p in pocket_alt_positions]
        rows.append({
            "model": alt_path.name,
            "canonical_length": len(canon),
            "alternate_length": len(alt),
            "alignment_score": score,
            "mapped_seed_residues": len(pocket_alt_positions),
            "local_ca_count": len(usable),
            "local_ca_rmsd_A": ca_rmsd,
            "local_shared_heavy_atom_rmsd_A": float(np.sqrt(np.mean(atom_sq))) if atom_sq else None,
            "pocket_mean_plddt": float(np.mean(plddt)) if plddt else None,
            "pocket_min_plddt": float(np.min(plddt)) if plddt else None,
            "position_map": ";".join(f"{p}:{mapping[p]}" for p in positions if p in mapping),
        })
    return rows


def gfra_interface_metrics(canonical_path, alt_paths, complex_path):
    canonical = residues(canonical_path)
    complex_gfra = residues(complex_path, "A")
    complex_gdnf = residues(complex_path, "B")
    can_to_complex, _ = seq_map(canonical, complex_gfra)

    interface_complex_positions = []
    min_dist_by_complex = {}
    gdnf_atoms = np.asarray([a.coord for r in complex_gdnf for a in r if not a.element.startswith("H")])
    for i, r in enumerate(complex_gfra, start=1):
        xyz = np.asarray([a.coord for a in r if not a.element.startswith("H")])
        dmin = float(np.sqrt(((xyz[:, None, :] - gdnf_atoms[None, :, :]) ** 2).sum(axis=2)).min())
        min_dist_by_complex[i] = dmin
        if dmin <= 5.0:
            interface_complex_positions.append(i)
    complex_to_can = {v: k for k, v in can_to_complex.items()}
    interface_can = sorted(complex_to_can[p] for p in interface_complex_positions if p in complex_to_can)

    # Anchor the domain with resolved positions surrounding the interface, then measure
    # interface-residue deviations. Positions 140-144 are the human deletion itself.
    anchor_can = sorted(p for p in can_to_complex if 145 <= p <= 350 and "CA" in canonical[p - 1])
    rows = []
    for alt_path in alt_paths:
        alt = residues(alt_path)
        can_to_alt, score = seq_map(canonical, alt)
        anchors = [p for p in anchor_can if p in can_to_alt and "CA" in alt[can_to_alt[p] - 1]]
        c_xyz = ca_array(canonical, anchors)
        a_xyz = np.asarray([alt[can_to_alt[p] - 1]["CA"].coord for p in anchors])
        rot, tran, domain_rmsd = kabsch_transform(a_xyz, c_xyz)
        interface_devs = []
        interface_plddt = []
        for p in interface_can:
            if p not in can_to_alt:
                continue
            moved = alt[can_to_alt[p] - 1]["CA"].coord @ rot.T + tran
            interface_devs.append(float(np.linalg.norm(moved - canonical[p - 1]["CA"].coord)))
            interface_plddt.append(float(alt[can_to_alt[p] - 1]["CA"].bfactor))
        deletion_neighbor = [p for p in range(135, 156) if p in can_to_alt]
        neighbor_plddt = [float(alt[can_to_alt[p] - 1]["CA"].bfactor) for p in deletion_neighbor]
        rows.append({
            "model": alt_path.name,
            "canonical_length": len(canonical),
            "alternate_length": len(alt),
            "alignment_score": score,
            "domain_anchor_ca_rmsd_A": domain_rmsd,
            "n_interface_residues": len(interface_devs),
            "interface_ca_displacement_mean_A": float(np.mean(interface_devs)),
            "interface_ca_displacement_max_A": float(np.max(interface_devs)),
            "interface_mean_plddt": float(np.mean(interface_plddt)),
            "deletion_neighbor_mean_plddt": float(np.mean(neighbor_plddt)),
        })
    meta = {
        "complex": str(complex_path),
        "complex_gfra_chain": "A",
        "complex_gdnf_chain": "B",
        "contact_cutoff_A": 5.0,
        "experimental_interface_human_positions": [int(p) for p in interface_can],
        "experimental_interface_position_count": len(interface_can),
        "deleted_human_positions": [140, 141, 142, 143, 144],
        "deleted_residues_resolved_in_complex": [int(p) for p in range(140, 145) if p in can_to_complex],
        "first_human_position_mapped_to_complex": int(min(can_to_complex)) if can_to_complex else None,
    }
    return rows, meta


def sort1_structure_meta(canonical_path, alt_paths, experimental_path):
    canonical = residues(canonical_path)
    experimental = residues(experimental_path, "A")
    can_to_exp, _ = seq_map(canonical, experimental)
    deleted = list(range(1, 138))
    mature_deleted = list(range(78, 138))
    epitope_sequences = {}
    for alt_path in alt_paths:
        alt = residues(alt_path)
        can_to_alt, _ = seq_map(canonical, alt)
        epitope_sequences[alt_path.name] = {
            "mapped_epitope_positions": sum(p in can_to_alt for p in SORT1_EPITOPE_CANONICAL),
            "mapped_callout_positions": sum(p in can_to_alt for p in SORT1_EPITOPE_CALLOUTS),
            "epitope_sequence_identical": all(
                p in can_to_alt and canonical[p - 1].resname == alt[can_to_alt[p] - 1].resname
                for p in SORT1_EPITOPE_CANONICAL
            ),
        }
    return {
        "experimental_structure": str(experimental_path),
        "experimental_chain": "A",
        "canonical_deleted_positions": deleted,
        "mature_beta_propeller_deleted_positions": mature_deleted,
        "deleted_positions_resolved_in_experiment": [p for p in mature_deleted if p in can_to_exp],
        "patent_epitope_interval": [207, 231],
        "patent_epitope_callouts": SORT1_EPITOPE_CALLOUTS,
        "alternate_epitope_mapping": epitope_sequences,
        "interpretation": (
            "The direct latozinemab epitope is retained in sequence. The 1-137 deletion is "
            "outside that epitope but removes the signal peptide and the first 60 residues of "
            "the mature beta-propeller; SORT1-209 is therefore a C/trafficking hypothesis, not B."
        ),
    }


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows):
    numeric = [k for k, v in rows[0].items() if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return {k: {"mean": float(np.mean([r[k] for r in rows if r[k] is not None])),
                "min": float(np.min([r[k] for r in rows if r[k] is not None])),
                "max": float(np.max([r[k] for r in rows if r[k] is not None]))}
            for k in numeric if any(r[k] is not None for r in rows)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gfra-complex", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("outputs/docking_campaign/analysis/bc_candidates"))
    args = parser.parse_args()
    root = Path.cwd()
    gfra_alts = sorted(root.glob(GFRA1_ALT_GLOB))
    cac_alts = sorted(root.glob(CACNA1D_ALT_GLOB))
    sort1_alts = sorted(root.glob(SORT1_ALT_GLOB))
    gfra_rows, gfra_meta = gfra_interface_metrics(GFRA1_CANONICAL, gfra_alts, args.gfra_complex)
    cac_rows = local_metrics(CACNA1D_CANONICAL, cac_alts, CACNA1D_POCKET_CANONICAL)
    sort1_rows = local_metrics(SORT1_CANONICAL, sort1_alts, SORT1_EPITOPE_CANONICAL, flank=2)
    sort1_meta = sort1_structure_meta(
        SORT1_CANONICAL,
        sort1_alts,
        Path("outputs/docking_campaign/systems/SORT1_latozinemab/inputs/4PO7.pdb"),
    )
    write_csv(args.outdir / "GFRA1_209_interface_model_metrics.csv", gfra_rows)
    write_csv(args.outdir / "CACNA1D_214_pocket_model_metrics.csv", cac_rows)
    write_csv(args.outdir / "SORT1_209_epitope_model_metrics.csv", sort1_rows)
    summary = {
        "GFRA1_209": {"metrics": summarize(gfra_rows), "experimental_mapping": gfra_meta},
        "CACNA1D_214": {"metrics": summarize(cac_rows), "pocket_seed_positions": CACNA1D_POCKET_CANONICAL},
        "SORT1_209": {"metrics": summarize(sort1_rows), "experimental_mapping": sort1_meta},
    }
    (args.outdir / "bc_structural_validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

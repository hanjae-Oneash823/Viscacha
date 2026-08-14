"""M2b — turn M2's raw structures into a confidence read for the
isoform-altered region specifically, not a whole-protein average.

Implements the isoform-specific QC protocol (see docs/MASTER_SURVEYOR_plan.md
and this session's design discussion):
  5. per-span pLDDT + PAE (region_confidence) -- a high WHOLE-PROTEIN pLDDT
     can hide a genuinely unreliable altered region.
  6. ensemble spread across the multi-seed ColabFold models -- the altered
     region is exactly where you'd expect model-to-model variability if the
     fold there is genuinely ambiguous; a tight ensemble is more trustworthy
     than one high-scoring model.
  4. ColabFold-vs-ESMFold local agreement over that same span -- two
     independent methods (MSA-based vs MSA-free) disagreeing is itself a
     confidence signal, not something to average away.
  7. alt-vs-canonical superposition + local RMSD across the sequence -- the
     actual "where does the isoform structurally diverge" signal.

Scope note on domain overlap: hits_deep.csv's `affected_domain` column
already encodes whether changed_aa_start:changed_aa_end overlaps a known
Pfam domain (computed in junior_surveyor/j2_protein_diff.py); this module
does not re-derive domain *coordinates* (that would mean re-running the
Pfam/HMMER scan), so `assess()` just carries that existing flag alongside
the new RMSD numbers rather than independently recomputing a geometric
overlap.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from statistics import mean

import numpy as np
from Bio.PDB import PDBParser, Superimposer

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from master_surveyor.align_utils import canonical_span_to_alt, match_anchors, residue_correspondence

_PARSER = PDBParser(QUIET=True)

# Starting heuristic thresholds for the low/medium/high verdict -- tune once
# real per-hit numbers exist (Phase 4 verification); intentionally kept as
# named constants here, not scattered magic numbers, so that tuning is a
# one-line change.
_HIGH_MAX_ENSEMBLE_RMSD   = 1.5
_HIGH_MAX_CF_ESMFOLD_RMSD = 3.0
_HIGH_MIN_REGION_PLDDT    = 70.0
_MED_MAX_ENSEMBLE_RMSD    = 3.0
_MED_MAX_CF_ESMFOLD_RMSD  = 5.0
_MED_MIN_REGION_PLDDT     = 50.0


def _log(msg: str) -> None:
    print(f"  [m2b] {msg}", file=sys.stderr, flush=True)


def _ca_atoms_by_residue(pdb_path: Path) -> dict[int, "Bio.PDB.Atom.Atom"]:
    """1-indexed sequential residue position (file order, matching the plain
    1..len(sequence) numbering ColabFold/ESMFold both write) -> CA atom.
    """
    structure = _PARSER.get_structure(pdb_path.stem, str(pdb_path))
    chain = next(iter(structure[0]))
    return {i: res["CA"] for i, res in enumerate(chain, start=1) if "CA" in res}


def _plddt_by_residue(pdb_path: Path) -> dict[int, float]:
    """pLDDT is stored as the CA atom's B-factor in both ColabFold's and
    ESMFold's output PDBs, but NOT on the same scale: ColabFold follows the
    standard AlphaFold convention (0-100), while transformers'
    EsmForProteinFolding.infer_pdb() writes it as a 0-1 fraction (verified
    empirically -- a test prediction came back with B-factors in [0.32,
    0.94], never above 1.0). Both are rescaled to 0-100 here so every
    downstream number (region_confidence, the low/medium/high verdict
    thresholds) means the same thing regardless of which engine produced
    the structure. Heuristic: real pLDDT is essentially never <=1.0 for
    every atom in a structure, so max(bfactor) <= 1.0 reliably identifies
    the fractional-scale case.
    """
    raw = {i: atom.get_bfactor() for i, atom in _ca_atoms_by_residue(pdb_path).items()}
    if raw and max(raw.values()) <= 1.0:
        return {i: v * 100.0 for i, v in raw.items()}
    return raw


def region_confidence(pdb_path: Path, start: int, end: int, pae_json_path: Path | None = None) -> dict:
    """start/end are 1-indexed, inclusive, in THIS structure's own residue
    numbering (i.e. already alt-local if pdb_path is an alt structure --
    see align_utils.canonical_span_to_alt for mapping from hits_deep.csv's
    canonical-coordinate changed_aa_start/changed_aa_end).
    """
    plddt = _plddt_by_residue(pdb_path)
    span_vals = [v for i, v in plddt.items() if start <= i <= end]
    whole_vals = list(plddt.values())

    result = {
        "region_mean_plddt": round(mean(span_vals), 2) if span_vals else None,
        "whole_protein_mean_plddt": round(mean(whole_vals), 2) if whole_vals else None,
        "region_length": len(span_vals),
    }

    if pae_json_path and pae_json_path.exists():
        payload = json.loads(pae_json_path.read_text())
        pae = payload.get("pae")
        if pae:
            pae_arr = np.array(pae, dtype=float)
            n = pae_arr.shape[0]
            region_idx = [i - 1 for i in range(start, end + 1) if 0 <= i - 1 < n]
            rest_idx = [i for i in range(n) if i not in set(region_idx)]
            if region_idx and rest_idx:
                block = pae_arr[np.ix_(region_idx, rest_idx)]
                result["region_vs_rest_mean_pae"] = round(float(block.mean()), 2)

    return result


def local_rmsd_across_models(model_paths: list[Path], start: int, end: int) -> float | None:
    """Mean pairwise CA RMSD (vs the first path) restricted to [start, end],
    for models that all describe the SAME sequence with the SAME residue
    numbering -- used both for the ColabFold seed-ensemble spread and for
    the ColabFold-vs-ESMFold agreement check (same computation, same
    same-sequence assumption, just a different pair of inputs).
    """
    if len(model_paths) < 2:
        return None

    span_atom_sets = []
    for path in model_paths:
        atoms = _ca_atoms_by_residue(path)
        span = [atoms[i] for i in range(start, end + 1) if i in atoms]
        span_atom_sets.append(span)

    lengths = {len(s) for s in span_atom_sets}
    if len(lengths) != 1 or 0 in lengths:
        _log(f"local_rmsd_across_models: span length mismatch across models {lengths} -- skipping")
        return None

    ref = span_atom_sets[0]
    rmsds = []
    for other in span_atom_sets[1:]:
        sup = Superimposer()
        sup.set_atoms(ref, other)
        rmsds.append(sup.rms)
    return round(float(np.mean(rmsds)), 3) if rmsds else None


def superpose_alt_vs_canonical(
    alt_pdb: Path, canonical_pdb: Path,
    canonical_seq: str, alt_seq: str,
    changed_aa_start: int, changed_aa_end: int,
) -> dict:
    """Fit alt onto canonical using ONLY exact-sequence-match ("=") anchor
    residues (canonical and alt diverge in length around indels, so a raw
    1:1 residue-index superposition breaks the moment the two sequences go
    out of register) -- then, under that same rigid-body transform, measure
    how far every aligned position's CA lands from its canonical
    counterpart. This is the actual "where does the isoform structurally
    diverge" signal (protocol point 7), not just a global RMSD number.
    """
    columns = residue_correspondence(canonical_seq, alt_seq)
    anchors = match_anchors(columns)   # [(canon_pos, alt_pos), ...] for "=" columns only
    if len(anchors) < 3:
        return {"error": "fewer than 3 exact-match anchor residues -- superposition not meaningful"}

    canon_atoms_all = _ca_atoms_by_residue(canonical_pdb)
    alt_atoms_all = _ca_atoms_by_residue(alt_pdb)

    anchor_canon_atoms = [canon_atoms_all[c] for c, a in anchors if c in canon_atoms_all and a in alt_atoms_all]
    anchor_alt_atoms = [alt_atoms_all[a] for c, a in anchors if c in canon_atoms_all and a in alt_atoms_all]
    if len(anchor_canon_atoms) < 3:
        return {"error": "fewer than 3 anchor residues resolved in both PDB structures"}

    sup = Superimposer()
    sup.set_atoms(anchor_canon_atoms, anchor_alt_atoms)
    anchor_fit_rmsd = sup.rms

    all_alt_atoms = list(alt_atoms_all.values())
    sup.apply(all_alt_atoms)   # rotate+translate every alt CA into the canonical frame, in place

    in_span, outside_span = [], []
    for col in columns:
        if col.canon_pos is None or col.alt_pos is None:
            continue   # pure indel position -- no canonical counterpart to measure against
        if col.canon_pos not in canon_atoms_all or col.alt_pos not in alt_atoms_all:
            continue
        dist = float(np.linalg.norm(
            canon_atoms_all[col.canon_pos].get_coord() - alt_atoms_all[col.alt_pos].get_coord()
        ))
        (in_span if changed_aa_start <= col.canon_pos <= changed_aa_end else outside_span).append(dist)

    return {
        "n_anchor_residues": len(anchor_canon_atoms),
        "anchor_fit_rmsd": round(anchor_fit_rmsd, 3),
        "changed_region_mean_local_dist": round(mean(in_span), 2) if in_span else None,
        "changed_region_max_local_dist": round(max(in_span), 2) if in_span else None,
        # Propagation signal: do SEQUENCE-IDENTICAL positions outside the
        # edited span still land far from canonical after the fit? That
        # would mean the change perturbs the fold beyond the primary edit
        # (e.g. an allosteric shift), which a naive "only look at the
        # changed positions" check would miss entirely.
        "outside_span_mean_local_dist": round(mean(outside_span), 2) if outside_span else None,
        "outside_span_max_local_dist": round(max(outside_span), 2) if outside_span else None,
    }


def _verdict(ensemble_rmsd: float | None, cf_esmfold_rmsd: float | None, region_plddt: float | None) -> str:
    if region_plddt is None:
        return "unknown"
    if (ensemble_rmsd is not None and ensemble_rmsd <= _HIGH_MAX_ENSEMBLE_RMSD
            and (cf_esmfold_rmsd is None or cf_esmfold_rmsd <= _HIGH_MAX_CF_ESMFOLD_RMSD)
            and region_plddt >= _HIGH_MIN_REGION_PLDDT):
        return "high"
    if (ensemble_rmsd is not None and ensemble_rmsd <= _MED_MAX_ENSEMBLE_RMSD
            and (cf_esmfold_rmsd is None or cf_esmfold_rmsd <= _MED_MAX_CF_ESMFOLD_RMSD)
            and region_plddt >= _MED_MIN_REGION_PLDDT):
        return "medium"
    return "low"


def assess(
    alt_result, canonical_path: Path,
    canonical_seq: str, alt_seq: str,
    changed_aa_start: int, changed_aa_end: int,
) -> dict:
    """Top-level orchestrator: takes an m2_structures.AltStructureResult plus
    the canonical structure path/sequences, returns one structure_confidence
    dict ready for m3_export.py's manifest.json.
    """
    columns = residue_correspondence(canonical_seq, alt_seq)
    alt_span = canonical_span_to_alt(columns, changed_aa_start, changed_aa_end)
    if alt_span is None:
        _log("no alt-local span could be mapped from changed_aa_start/end -- confidence assessment skipped")
        return {"error": "no alt-local span mapped from canonical coordinates"}
    alt_start, alt_end = alt_span

    top_model = alt_result.seed_models[0]
    top_scores = alt_result.scores_json[0] if alt_result.scores_json else None
    region = region_confidence(top_model, alt_start, alt_end, top_scores)

    ensemble_rmsd = local_rmsd_across_models(alt_result.seed_models, alt_start, alt_end)

    cf_esmfold_rmsd = None
    if alt_result.esmfold_model is not None:
        cf_esmfold_rmsd = local_rmsd_across_models([top_model, alt_result.esmfold_model], alt_start, alt_end)

    superposition = superpose_alt_vs_canonical(
        top_model, canonical_path, canonical_seq, alt_seq, changed_aa_start, changed_aa_end,
    )

    verdict = _verdict(ensemble_rmsd, cf_esmfold_rmsd, region.get("region_mean_plddt"))

    return {
        "verdict": verdict,
        "region_confidence": region,
        "ensemble_spread_rmsd": ensemble_rmsd,
        "n_ensemble_models": len(alt_result.seed_models),
        "colabfold_vs_esmfold_rmsd": cf_esmfold_rmsd,
        "esmfold_error": alt_result.esmfold_error,
        "superposition": superposition,
        "alt_local_span": [alt_start, alt_end],
    }

#!/usr/bin/env python3
"""Generate a concise, evidence-bounded preliminary docking handoff."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "docking_campaign"
SYSTEMS = OUT / "systems"


def main() -> None:
    fyn = json.loads((SYSTEMS / "FYN_saracatinib" / "analysis" / "summary.json").read_text())
    kit = json.loads((SYSTEMS / "KIT_masitinib" / "analysis" / "masitinib_1T46_summary.json").read_text())
    cac_result = json.loads((SYSTEMS / "CACNA1D_isradipine" / "runs" / "amiodarone_corrected_seed20260825_ex32" / "result.json").read_text())
    cac_top = cac_result["results"][0]
    best = fyn["best_top_pose"]

    report = f"""# Wednesday preliminary docking handoff

## Take-home result

The AutoDock Vina workflow has one successful, reproducible crystal-ligand re-docking example: **saracatinib in the canonical FYN kinase domain**. Across {fyn['completed_runs']} independently seeded runs, the top pose was within 2.0 Å of the crystal pose in {fyn['top_pose_recovered_under_2A']}/{fyn['completed_runs']} runs (median fixed-frame heavy-atom RMSD {fyn['top_pose_rmsd_median_A']:.2f} Å; best {best['rmsd_to_crystal_heavy_atom_uncorrected_angstrom']:.2f} Å; best Vina score {best['vina_affinity_kcal_mol']:.3f} kcal/mol). This validates the *canonical FYN* receptor/ligand/grid protocol, not cellular activity or an isoform effect.

![FYN preliminary result](../../systems/FYN_saracatinib/figures/fyn_saracatinib_preliminary.png)

## Canonical c-KIT–masitinib baseline

Using the experimental imatinib-bound c-KIT kinase-domain pocket (PDB 1T46), masitinib gave a stable canonical baseline across {kit['n_independent_vina_runs']} Vina seeds: mean {kit['score_mean_kcal_mol']:.3f} kcal/mol, median {kit['score_median_kcal_mol']:.3f} kcal/mol, SD {kit['score_sd_kcal_mol']:.3f} kcal/mol, range {kit['score_range_kcal_mol'][0]:.3f} to {kit['score_range_kcal_mol'][1]:.3f} kcal/mol.

No KIT-223 comparison score is reported. The KIT-223 single-residue deletion is at canonical residue 715, while the relevant kinase-insert segment (residues 690–761) is unresolved in 1T46. A meaningful comparison requires a rebuilt, locally validated loop model; the initial predicted-model pocket runs were therefore excluded rather than over-interpreted.

![KIT preliminary result](../../systems/KIT_masitinib/figures/kit_masitinib_preliminary.png)

## Full candidate triage and CACNA1D control

The current AD-enriched coding-isoform data support three candidates from the original list: FYN, KIT, and CACNA1D-214. The remaining proposed comparisons are excluded from quantitative docking in this handoff:

- **GABRA2-206:** not an AD-enriched coding isoform in the available hit table (the matching gene row is GABRA2-202 and is control-enriched).
- **BACE1-476 / BACE1-457:** neither alternate is present in the available hit table, so the expression prerequisite for the proposed comparison is absent.
- **CHRNA7/CHRFAM7A:** no CHRFAM7A transcript/genotype evidence is available here; a mixed-pentamer calculation would be unanchored.
- **PDE9A:** the candidate list does not specify an AD-enriched coding-altered transcript; no such alternate is available to model.

**CACNA1D-214** is AD-enriched and encodes a C-terminally truncated protein (1,625 aa versus canonical 2,161 aa). The experimental human CaV1.3 template 8E59 contains the amiodarone-bound channel through residue 1,589, entirely before the truncation. This supports a shared-pocket/pocket-retention interpretation, rather than an isoform-specific affinity claim. A single 16-core amiodarone re-docking attempt did not validate the template with Vina (top score {cac_top['vina_affinity_kcal_mol']:.3f} kcal/mol; top fixed-frame RMSD {cac_top['rmsd_to_crystal_heavy_atom_uncorrected_angstrom']:.2f} Å), so no CACNA1D drug score is carried forward.

## Methods to state on Wednesday

- AutoDock Vina scoring; seeded independent runs; 16 CPU cores per run and no remaining concurrent docking jobs.
- FYN validation used human FYN–saracatinib crystal structure PDB 10DJ, chain A; grid center (-11.255, 14.853, -9.445) Å and box 20 × 20 × 26 Å; exhaustiveness 32.
- FYN RMSD is fixed-frame heavy-atom, based on preserved PDBQT atom order and not symmetry-corrected.
- c-KIT results are a canonical-pocket feasibility/baseline result, not a validated protein-isoform affinity comparison.
- The CaV1.3 8E59/amiodarone re-docking control failed Vina pose recovery and is explicitly excluded from quantitative interpretation.
- GNINA was staged but cannot currently run because its CUDA binary requires `libcudnn.so.9`; Vina is the reported engine for this preliminary handoff.

## Files

- FYN machine-readable summary: `../../systems/FYN_saracatinib/analysis/summary.json`
- FYN figure (PNG/PDF): `../../systems/FYN_saracatinib/figures/fyn_saracatinib_preliminary.*`
- c-KIT machine-readable summary: `../../systems/KIT_masitinib/analysis/masitinib_1T46_summary.json`
- c-KIT figure (PNG/PDF): `../../systems/KIT_masitinib/figures/kit_masitinib_preliminary.*`
- CaV1.3 corrected staging and control: `../../systems/CACNA1D_isradipine/prepared/stage_metadata.json` and `../../systems/CACNA1D_isradipine/runs/amiodarone_corrected_seed20260825_ex32/result.json`

## Next defensible analysis

1. Build and validate a KIT kinase-insert loop model that explicitly includes residue 715 before attempting KIT-223 docking.
2. Prepare an explicit canonical/alternate FYN structural comparison only after confirming the alternate transcript’s translated product and biological relevance.
3. Add BACE1, CHRFAM7A, GABRA2-206, or PDE9A only after their required alternate transcript/genotype evidence is supplied.
4. If the CUDA/cuDNN environment is repaired, rescore retained Vina poses with GNINA CNNscore/CNNaffinity as an orthogonal ranking, while retaining Vina redocking as the pose-validation benchmark.
"""
    destination = OUT / "docs" / "archive" / "WEDNESDAY_PRELIMINARY_DOCKING_REPORT.md"
    destination.write_text(report)
    print(destination)


if __name__ == "__main__":
    main()

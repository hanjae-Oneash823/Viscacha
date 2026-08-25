# Expanded docking campaign: key findings

## One-sentence result

Every candidate row was attempted, but only BACE1 produced a new, quality-controlled canonical-versus-splice-isoform comparison; several other rows are useful structural exclusions or exploratory controls rather than direct affinity comparisons.

## Strongest new result: BACE1–verubecestat

The canonical control worked extremely well. Across six seeds, Vina placed verubecestat back into its experimental pose with a mean RMSD of **0.955 Å**, and all six top poses passed the 2 Å criterion.

| Protein | Mean top score ± SD | Mean top-pose RMSD | Change from canonical |
|---|---:|---:|---:|
| BACE1 canonical | −9.510 ± 0.021 kcal/mol | 0.955 Å | reference |
| BACE1-476 | −9.171 ± 0.013 kcal/mol | 7.122 Å | +0.339 kcal/mol, less favorable |
| BACE1-457 | −8.285 ± 0.016 kcal/mol | 7.013 Å | +1.225 kcal/mol, less favorable |

In easy language: both deletion isoforms pushed verubecestat into a very different pose. The 457-aa isoform also produced a clearly weaker Vina score than canonical BACE1. This is the strongest new structural observation, especially because BACE1-457 deletes five residues that contact verubecestat in the experimental structure.

This result is still model-based. It supports a follow-up hypothesis; it is not proof that either isoform changes drug response in patients.

## CHRNA7–encenicline: a strong computational effect with a major biological caveat

Canonical α7 redocking was excellent: six of six top poses recovered the experimental ligand orientation, with mean RMSD **0.297 Å** and mean score **−9.864 kcal/mol**.

Two explicit one-fusion-subunit hypotheses gave:

| Receptor hypothesis | Mean top score ± SD | Mean RMSD to canonical 7EKP pose | Change from canonical |
|---|---:|---:|---:|
| Canonical α7 | −9.864 ± 0.009 | 0.297 Å | reference |
| Fusion at A face | −7.537 ± 0.013 | 2.430 Å | +2.327 kcal/mol |
| Fusion at B face | −6.599 ± 0.348 | 8.450 Å | +3.264 kcal/mol |

In easy language: the result changes a lot depending on which side of the binding site contains the fusion subunit. That makes the geometry interesting, but it also proves that stoichiometry and orientation must be known before making a biological claim. The current samples do not yet provide that support.

## Other rows

### FYN–saracatinib

This remains the cleanest preliminary presentation result. Canonical redocking was validated, while the 115-aa AD-associated alternate lacks the entire kinase pocket. The correct interpretation is **loss of the drug target domain**, not a weaker docking score.

### KIT–masitinib

The experimental canonical baseline is stable at **−12.817 ± 0.032 kcal/mol**. A valid KIT-223 comparison was not obtained. Ser715 lies in a segment missing from the experimental 1T46 structure, and predicted-model attempts produced severe ligand clashes and positive Vina scores. Those numbers were rejected.

The defensible conclusion is that KIT-223 was attempted but remains unresolved—not that it binds masitinib better or worse.

### GABRA2–AZD7325

The accepted native α2/γ2 interface passed two structural-transfer checks and gave a reproducible exploratory canonical score of **−7.375 ± 0.009 kcal/mol**. GABRA2-206 is only 73 aa and cannot form the extracellular/transmembrane pentamer, so an alternate docking score would be meaningless.

There is also a data-label warning: the local gate table does not reproduce the “AD-only 9.3% versus 0%” statement in the candidate file. The transcript direction must be rechecked before presenting this as an AD-enriched alternate.

### CACNA1D–isradipine

The human experimental structure ends before the alternate's C-terminal truncation, so canonical CaV1.3 and CACNA1D-214 have identical modeled pocket coordinates. The paired pilot runs were exactly identical, as expected. Six exploratory seeds gave **−6.547 ± 0.185 kcal/mol** for the shared pocket.

The corrected amiodarone control found a near-native pose at rank 7 rather than rank 1, so this protocol is not fully validated for CaV1.3. Present this as a pocket-preservation negative control, not an affinity prediction.

### PDE9A–BI 409306

Canonical cross-docking was stable at **−8.082 ± 0.033 kcal/mol**. No alternate comparison was run because the candidate file does not name the coding-altered PDE9A transcript. A specific sequence is required before this row can be finished.

## What is ready to present

1. **FYN:** validated canonical redocking plus complete kinase-pocket loss in the alternate.
2. **BACE1:** validated canonical redocking plus reproducible pose shifts in both deletion isoforms, with the largest score penalty in BACE1-457.
3. **Campaign-status figure:** demonstrates that all eight rows were attempted and explains why some do not have alternate scores.
4. **CHRNA7 hypothesis figure:** useful as a clearly labeled exploratory result showing topology dependence.
5. **GABRA2 and CACNA1D:** useful structural controls, with the stated caveats.

## What not to claim

- Do not compare absolute Vina scores between different proteins or drugs.
- Do not call the BACE1 or CHRFAM7A results experimental binding affinities.
- Do not claim a KIT-223 affinity difference from the rejected predicted-model runs.
- Do not claim GABRA2-206 is AD-enriched until the usage-direction discrepancy is resolved.
- Do not choose a PDE9A alternate without an explicit transcript/protein sequence.
- Do not interpret a pocket-absent alternate as a low-affinity numerical docking result.

## Best new figures

- `standalone_figures/expanded_campaign/matched_comparative_docking.png`
- `standalone_figures/expanded_campaign/all_candidate_status.png`
- `standalone_figures/expanded_campaign/canonical_crossdock_stability.png`
- `standalone_figures/expanded_campaign/3D_BACE1_variant_pocket_overlay.png`
- `standalone_figures/expanded_campaign/3D_CHRNA7_topology_site_overlay.png`

# Protein–drug docking campaign: detailed AutoDock Vina and GNINA findings

> **Biological qualification added 2026-08-26:** the BACE1-476 and BACE1-457 calculations below are structural method-development comparisons; those exact isoforms were not confirmed as the AD-increased transcript in the available DTU table. For presentation candidate selection, use [BC_CANDIDATE_ASSESSMENT.md](BC_CANDIDATE_ASSESSMENT.md), which analyzes the observed 401-aa BACE1-202 isoform and the statistically supported CACNA1D-214 candidate.

**Scope:** all eight candidate rows
**Accepted expanded replicate records:** 54
**GNINA rescored run ensembles / poses:** 45 / 738
**Maximum CPU allocation:** 16 cores
**Document status:** current authoritative findings report

## 1. Executive conclusion

All eight proposed rows were attempted, but the scientifically appropriate endpoint differed by candidate. The campaign produced three strong canonical pose-recovery controls (FYN–saracatinib, BACE1–verubecestat, and CHRNA7–encenicline), one quality-controlled modeled splice-isoform comparison (BACE1), one topology-dependent receptor hypothesis (CHRNA7/CHRFAM7A), three exploratory canonical/shared-pocket controls (GABRA2, CACNA1D, and PDE9A), and three important structural or definition-based exclusions (FYN alternate, GABRA2-206, and the unspecified PDE9A alternate). KIT-223 was explicitly modeled and docked but rejected because the altered local structure is unresolved experimentally and the predicted pockets generated severe steric incompatibility.

The strongest new comparative finding is that both modeled BACE1 deletion isoforms reproducibly displaced verubecestat from the canonical crystallographic binding mode. BACE1-457 also showed a paired Vina score penalty of +1.225 ± 0.023 kcal/mol relative to canonical BACE1. This supports a structural hypothesis that the BACE1-457 deletion perturbs verubecestat recognition, but it remains model-derived and is not an experimental affinity measurement.

The strongest canonical validation remains FYN–saracatinib: all nine independent top-ranked poses were below 2.0 Å RMSD. The most precise new control is CHRNA7–encenicline, with a mean top-pose RMSD of 0.297 Å across six runs.

GNINA 1.3.3 rescoring added an orthogonal scoring view over 738 retained Vina poses from 45 run ensembles. It strongly reinforced both BACE1 deletion hypotheses and the severe CHRFAM7A B-face effect. It preserved acceptable canonical pose recovery for BACE1 (6/6) and CHRNA7 (6/6), but not for FYN (0/9 after CNN reranking). The FYN disagreement is reported explicitly and prevents describing the three canonical controls as uniformly supported by both scoring methods.

## 2. Confidence framework

The findings are separated by what the underlying calculation can support.

| Confidence category | Required evidence | Findings in this campaign |
|---|---|---|
| Validated canonical protocol | Exact cognate redocking with reproducible top-pose recovery below 2 Å | FYN, BACE1, CHRNA7 |
| Moderate modeled comparison | Validated canonical protocol plus alternate model with acceptable local pocket confidence and a fully matched docking design | BACE1-476 and BACE1-457 |
| Hypothesis only | Canonical validation exists, but alternate assembly/topology is assumed rather than measured | CHRNA7/CHRFAM7A hybrids |
| Exploratory | Cross-docking or a control without successful top-ranked cognate pose recovery | GABRA2, CACNA1D, PDE9A, canonical KIT baseline |
| Structural exclusion | The alternate cannot form the relevant pocket/receptor | FYN 115-aa product, GABRA2-206 |
| Unresolved | The alternate sequence or experimental local structure is inadequate | KIT-223 and PDE9A alternate |

## 3. Complete candidate disposition

| Row | System | Work completed | Retained conclusion | Confidence |
|---:|---|---|---|---|
| 1 | FYN–saracatinib | Nine-run exact canonical redocking; alternate domain analysis | Canonical protocol validated; alternate kinase pocket absent | High for canonical pose recovery and pocket-loss classification |
| 2 | KIT–masitinib | Nine-run experimental canonical baseline; whole-domain and local-refit alternate attempts | Canonical score reproducible; no valid KIT-223 comparison | Limited/exploratory |
| 3 | GABRA2–AZD7325 | Native pentamer selection, two-subunit interface-transfer validation, six canonical runs | Canonical cross-docking stable; 73-aa alternate cannot form receptor | Exploratory plus structural exclusion |
| 4 | CACNA1D–isradipine | Corrected amiodarone control; identical canonical/alternate pilot; six shared-pocket runs | Resolved pocket retained; no isoform-specific difference can be modeled from 8E59 | Exploratory negative control |
| 5 | BACE1-476–verubecestat | ColabFold model, local pocket QC, six paired canonical/alternate runs | Reproducible pose displacement; small score penalty | Moderate, model-dependent |
| 6 | BACE1-457–verubecestat | ColabFold model, local pocket QC, six paired canonical/alternate runs | Reproducible pose displacement and larger score penalty | Moderate, model-dependent |
| 7 | CHRNA7/CHRFAM7A–encenicline | Exact canonical control; two experimental-domain hybrid hypotheses; six paired runs per topology | Strong topology dependence, but no biological stoichiometry conclusion | Hypothesis only |
| 8 | PDE9A–BI 409306 | Six canonical cross-docking runs | Stable canonical exploratory baseline; alternate not specified | Incomplete comparison |

## 4. Canonical validation results

### 4.1 FYN–saracatinib

FYN was validated by redocking crystallographic saracatinib from PDB 10DJ into the chain-A kinase pocket. Nine independently seeded top poses all met the predefined 2.0 Å threshold.

| Metric | Result |
|---|---:|
| Independent runs | 9 |
| Top poses below 2.0 Å | 9/9 (100%) |
| Mean top-pose RMSD | 1.877 Å |
| Median top-pose RMSD | 1.896 Å |
| Best top-pose RMSD | 1.754 Å |
| Mean top score | −9.175 kcal/mol |
| Top-score SD | 0.025 kcal/mol |
| Best-geometry pose score | −9.182 kcal/mol |

This validates the specific canonical receptor preparation, saracatinib preparation, ligand-centered grid, and Vina search protocol for pose recovery. It does not validate a thermodynamic affinity, cellular response, or disease effect.

GNINA did not reproduce the FYN pose ranking. The Vina rank-1 poses had mean CNNscore 0.599 ± 0.027 and mean CNNaffinity 7.570 ± 0.027, but GNINA selected poses with mean original Vina rank 8.11 and mean RMSD 2.884 ± 0.248 Å. None of the nine CNN-selected poses remained below 2 Å. The selected poses are still in the same general pocket, but the ranking fails the predefined cognate threshold. This is a genuine method disagreement: FYN remains a successful Vina redocking example, not a Vina/GNINA consensus example.

The proposed alternate ends at residue 115. It lacks the entire SH2 and kinase domains and does not contain a complete SH3 domain. Because saracatinib binds the kinase ATP site, the alternate does not possess the molecular object being docked against. The correct result is **kinase pocket absent**, not a weak Vina score. This conclusion is structurally strong but biologically conditional on the novel transcript producing a stable protein.

### 4.2 BACE1–verubecestat

Canonical BACE1 exact redocking used PDB 5HU1 and verubecestat `66F A501`. All six accepted Open Babel-prepared canonical runs recovered the crystallographic pose.

| Metric | Result |
|---|---:|
| Independent runs | 6 |
| Top poses below 2.0 Å | 6/6 (100%) |
| Mean top-pose RMSD | 0.955 ± 0.006 Å |
| Mean top score | −9.510 ± 0.021 kcal/mol |
| Score range | −9.536 to −9.488 kcal/mol |

A separate Meeko-prepared pilot also passed (0.958 Å top-pose RMSD). The Open Babel series was retained because the same preparation method could be applied consistently to both alternate models. Thus, the matched comparison was anchored to a canonical receptor that passed independently under the chosen preparation route.

### 4.3 CHRNA7–encenicline

Canonical α7 exact redocking used the 7EKP homopentamer and encenicline `I33 A601` at the A/B intersubunit site.

| Metric | Result |
|---|---:|
| Independent runs | 6 |
| Top poses below 2.0 Å | 6/6 (100%) |
| Mean top-pose RMSD | 0.297 ± 0.012 Å |
| Mean top score | −9.864 ± 0.009 kcal/mol |
| Score range | −9.876 to −9.852 kcal/mol |

The Meeko-prepared pilot also passed at 0.310 Å. The Open Babel-prepared receptor series was selected for the matched hybrid comparison. The extremely small canonical RMSD indicates that the selected grid and preparation preserve the experimentally observed encenicline orientation under the rigid-receptor Vina protocol.

## 5. BACE1 splice-isoform comparison

### 5.1 Model evidence

BACE1-476 and BACE1-457 were predicted independently from reviewed UniProt sequences using one AlphaFold2-ptm model and one seed. Whole-model confidence was moderate, but the binding-site residues had substantially higher confidence.

| Alternate | Mean model pLDDT | pTM | Mapped 5HU1 contacts | Pocket Cα RMSD | Mean/min pocket pLDDT |
|---|---:|---:|---:|---:|---:|
| BACE1-476 | 77.75 | 0.722 | 28/28 | 1.484 Å | 86.87 / 73.19 |
| BACE1-457 | 78.70 | 0.750 | 23/28 | 1.937 Å | 85.27 / 69.81 |

The BACE1-457 deletion removes five of the 28 residues in the 6 Å experimental verubecestat contact shell. That is a direct sequence-to-pocket link and is the principal structural reason the 457-aa result deserves follow-up.

### 5.2 Matched docking results

All three receptors used the same verubecestat PDBQT, 5HU1 coordinate frame, box, exhaustiveness, pose count, and six seeds.

| Receptor | Mean top score ± SD | Score range | Mean top-pose RMSD ± SD | Paired score change from canonical |
|---|---:|---:|---:|---:|
| Canonical BACE1 | −9.510 ± 0.021 | −9.536 to −9.488 | 0.955 ± 0.006 Å | Reference |
| BACE1-476 | −9.171 ± 0.013 | −9.184 to −9.149 | 7.122 ± 0.010 Å | +0.339 ± 0.030 kcal/mol |
| BACE1-457 | −8.285 ± 0.016 | −8.305 to −8.265 | 7.013 ± 0.016 Å | +1.225 ± 0.023 kcal/mol |

Positive alternate-minus-canonical deltas indicate less favorable Vina scores under the matched protocol. Both alternate models consistently selected top poses approximately 7 Å from the canonical crystal orientation. Because the modeled receptors were aligned into the 5HU1 frame, this RMSD quantifies ligand displacement relative to the canonical pose; it is not a validation RMSD for an alternate experimental structure.

The BACE1-476 score change is small relative to the intrinsic limitations of docking scores, even though its pose displacement is reproducible. BACE1-457 shows both a reproducible pose shift and a larger score penalty. The appropriate conclusion is that BACE1-457 is the stronger mechanistic hypothesis, not that its experimental affinity is exactly 1.225 kcal/mol weaker.

GNINA independently reinforced the direction and ordering of the BACE1 result when the same Vina rank-1 pose was rescored for each paired seed. Relative to canonical BACE1, the mean alternate-minus-canonical CNNaffinity change was −1.674 ± 0.012 for BACE1-476 and −3.290 ± 0.023 for BACE1-457. The corresponding CNNscore changes were −0.804 ± 0.007 and −0.879 ± 0.007. GNINA selected the canonical near-crystal pose in 6/6 runs (0.955 ± 0.006 Å), while neither alternate produced a CNN-selected pose below 2 Å. These model outputs are not kcal/mol and are not experimental pKd values; their value is the concordant within-system ordering **canonical > BACE1-476 > BACE1-457**.

## 6. CHRNA7/CHRFAM7A topology hypotheses

### 6.1 Structural construction

The 9QTO fusion extracellular domain was aligned to each side of the 7EKP A/B binding site. Each alignment used 206 identical Cα positions and had 1.147 Å RMSD. No interfacial atom pair was closer than 1.5 Å. These checks made the models geometrically dockable but did not establish biological assembly.

### 6.2 Matched docking results

| Receptor hypothesis | Mean top score ± SD | Score range | Mean pose displacement ± SD | Paired score change |
|---|---:|---:|---:|---:|
| Canonical α7 | −9.864 ± 0.009 | −9.876 to −9.852 | 0.297 ± 0.012 Å | Reference |
| Fusion at A face | −7.537 ± 0.013 | −7.552 to −7.514 | 2.430 ± 1.003 Å | +2.327 ± 0.018 kcal/mol |
| Fusion at B face | −6.599 ± 0.348 | −6.927 to −6.217 | 8.450 ± 0.152 Å | +3.264 ± 0.344 kcal/mol |

The A-face hypothesis generally retained a near-canonical pose, although one seed generated a larger displacement and increased the RMSD SD. The B-face hypothesis consistently displaced encenicline by approximately 8.45 Å and showed greater score variability. Therefore, the modeled effect depends strongly on which face of the binding site contains the fusion domain.

This topology dependence is itself the key finding: an unspecified mixed-pentamer model cannot support a unique affinity prediction. Sample-specific CHRFAM7A genotype, transcript/protein expression, subunit stoichiometry, orientation, assembly, and surface localization are required before either numerical hypothesis can be assigned biological meaning.

GNINA also preserved the topology ordering. Rescoring the matched Vina rank-1 pose gave alternate-minus-canonical CNNaffinity changes of −0.433 ± 0.266 for the A-face model and −2.686 ± 0.295 for the B-face model. CNNscore changed by −0.195 ± 0.109 and −0.859 ± 0.043, respectively. GNINA selected a sub-2 Å pose for canonical CHRNA7 in 6/6 runs (mean 1.847 Å), but in 0/6 runs for either hybrid. The B-face perturbation is therefore supported by both Vina and GNINA; the A-face effect remains smaller and more variable.

## 7. KIT–masitinib: reproducible canonical baseline, rejected alternate

Masitinib cross-docking into the experimental imatinib-defined 1T46 pocket was numerically stable across nine canonical runs.

| Metric | Result |
|---|---:|
| Mean top score | −12.817 kcal/mol |
| Median top score | −12.814 kcal/mol |
| Sample SD | 0.032 kcal/mol |
| Range | −12.872 to −12.759 kcal/mol |

This is not cognate redocking because 1T46 contains imatinib, not masitinib. Consequently, the low score variance demonstrates repeatability only; no experimental masitinib pose was available for validation.

KIT-223 deletes Ser715 within a 1T46 segment unresolved from residues 690–761. Whole-domain predicted models and locally refitted pockets were both attempted. The local refits looked superficially good by Cα RMSD (0.692 Å canonical; 0.690 Å KIT-223) and mean local pLDDT (86.31 and 85.34), but the minimum local pLDDT values fell near 42. More decisively, refit docking gave positive top scores of +0.183 kcal/mol for the predicted canonical receptor and +16.529 kcal/mol for KIT-223, consistent with severe clashes or invalid pocket packing.

Those positive values are rejection diagnostics, not comparative affinity estimates. Reporting their difference would convert a model-preparation failure into a false biological result. The retained conclusion is: **canonical experimental baseline reproducible; KIT-223 affinity unresolved**.

## 8. GABRA2–AZD7325

The accepted 9CSB α2/γ2 site passed both sides of the interface-transfer check: 0.985 Å Cα RMSD for the α-chain fit and 0.740 Å for the neighboring γ2 chain after the same transformation. Six AZD7325 cross-docking runs gave:

| Metric | Result |
|---|---:|
| Mean top score | −7.375 ± 0.009 kcal/mol |
| Range | −7.389 to −7.363 kcal/mol |

The score is exploratory because AZD7325 was generated from SMILES and the site was transferred from a diazepam-bound homologous interface. It must not be described as validated pose recovery.

The proposed GABRA2-206 product is only 73 aa and cannot form the extracellular binding domain, transmembrane helices, or pentameric receptor architecture. Its endpoint is **receptor/pocket absent**, not a docking score.

There is also a transcript-label warning: the local gate table reviewed during triage did not reproduce the candidate file's stated AD-only direction. The disease-enrichment label should be reconciled against the source expression table before this row is presented as AD-specific.

## 9. CACNA1D–isradipine shared-pocket control

The original CaV1.3 control contained a ligand-code error: `3PE` is a phosphatidylethanolamine lipid, whereas amiodarone in 8E59 is `BBI A2201`. The control was restaged and rerun with `BBI`.

The corrected top-ranked amiodarone pose scored −1.041 kcal/mol but had 14.949 Å RMSD. A near-native pose was present at rank 7 with 1.560 Å RMSD. Thus, Vina sampled the experimental orientation but did not rank it correctly. The protocol failed the top-ranked pose-recovery criterion and was not considered fully validated.

The resolved 8E59 coordinates span CACNA1D residues 121–1589, while CACNA1D-214 ends at residue 1625. The deposited pocket therefore lies entirely before the truncation. Canonical and alternate pilot docking used exactly the same receptor coordinates and returned the same isradipine result (−6.707 kcal/mol for seed 20260825), as expected.

Across six shared-pocket seeds:

| Metric | Result |
|---|---:|
| Mean top score | −6.547 ± 0.185 kcal/mol |
| Range | −6.727 to −6.365 kcal/mol |

This is a pocket-preservation negative control. It indicates that the chosen experimental template cannot encode an isoform-specific pocket difference; it does not prove identical channel pharmacology, because the removed C terminus could affect gating, trafficking, regulation, or allosteric state outside a rigid local docking model.

## 10. PDE9A–BI 409306

BI 409306 was cross-docked into the 4GH6 catalytic pocket defined by ligand `LUO A601`. The protein-only receptor preparation excluded the crystallographic Zn²⁺ and Mg²⁺ ions, so the result does not represent a fully reconstructed PDE catalytic site. Six canonical runs gave:

| Metric | Result |
|---|---:|
| Mean top score | −8.082 ± 0.033 kcal/mol |
| Range | −8.137 to −8.045 kcal/mol |

The calculation is a stable canonical exploratory baseline, not cognate pose validation. Its missing catalytic metals further reduce mechanistic interpretability. The candidate row does not name the AD-enriched coding-altered PDE9A transcript or provide its protein sequence. Because PDE9A N-terminal variants may leave the catalytic domain unchanged, choosing an arbitrary isoform could manufacture a comparison unrelated to the proposed mechanism. No alternate score should be calculated until the exact sequence and catalytic-domain alteration are specified.

## 11. Cross-system interpretation

Raw Vina scores cannot be ranked across these proteins. Each system differs in receptor composition, pocket size and polarity, ligand size/flexibility, preparation, and search-space volume. For example, the canonical KIT score (−12.817 kcal/mol) must not be interpreted as stronger biological binding than BACE1 (−9.510 kcal/mol) or CHRNA7 (−9.864 kcal/mol). The defensible comparisons are within-system, preparation-matched differences and cognate pose-recovery metrics.

Likewise, score SD and RMSD answer different questions:

- low score SD indicates stochastic convergence under the specified search;
- low cognate RMSD indicates recovery of an experimental pose;
- an alternate-to-canonical RMSD indicates pose displacement in a shared frame, not validation of the alternate;
- neither metric is an experimental dissociation constant or free energy.

GNINA adds two more non-equivalent quantities. `CNNscore` ranks pose plausibility, whereas `CNNaffinity` is a machine-learning affinity prediction. Higher values are more favorable within the GNINA model, but neither is on the Vina kcal/mol scale. The primary consensus comparison rescored the same Vina rank-1 pose for each seed; the secondary analysis allowed GNINA to rerank the retained Vina ensemble. Because GNINA did not generate an independent pose ensemble, agreement constitutes scoring consensus rather than independent conformational validation.

## 12. Supported and unsupported claims

### 12.1 Supported

1. The canonical FYN, BACE1, and CHRNA7 protocols recover their deposited ligand poses reproducibly.
2. The proposed 115-aa FYN product lacks the saracatinib-binding kinase domain.
3. Both modeled BACE1 deletion isoforms reproducibly change the selected verubecestat pose; BACE1-457 shows the larger paired score penalty.
4. The modeled CHRFAM7A effect is strongly dependent on which face of the encenicline site contains the fusion domain.
5. Canonical KIT–masitinib cross-docking is numerically reproducible, but the KIT-223 comparison failed structural validity checks.
6. The proposed 73-aa GABRA2 alternate cannot form the AZD7325 receptor site.
7. The resolved CaV1.3/isradipine pocket is identical for canonical and CACNA1D-214 in the selected template.
8. PDE9A alternate docking remains undefined without a specific coding-altered sequence.
9. GNINA supports the direction and severity ordering of the BACE1 deletion effects and the CHRFAM7A topology dependence.
10. GNINA and Vina disagree on the preferred FYN pose; only the Vina protocol passes the FYN top-pose RMSD criterion.

### 12.2 Unsupported

- A Vina score is an experimental binding affinity or ΔG.
- BACE1-457 has exactly 1.225 kcal/mol weaker experimental binding.
- CHRFAM7A necessarily weakens encenicline binding in vivo.
- KIT-223 binds masitinib more weakly based on the rejected positive-score models.
- FYN or GABRA2 pocket absence proves that the predicted transcript is translated or stable.
- Identical CACNA1D pocket coordinates imply identical electrophysiological drug response.
- Any docking result demonstrates Alzheimer’s disease efficacy, brain exposure, clinical benefit, or causal mechanism.

## 13. Recommended follow-up experiments

### Highest priority: BACE1

1. Confirm BACE1-476 and BACE1-457 transcript enrichment in the relevant samples.
2. Verify protein expression, maturation, localization, and catalytic activity.
3. Generate structural ensembles or perform explicit-solvent molecular dynamics around the deletion-affected pocket.
4. Apply symmetry-aware RMSD, protein–ligand interaction fingerprints, and an orthogonal scoring method.
5. Measure verubecestat inhibition or binding experimentally for purified/expressed isoforms.

### CHRNA7/CHRFAM7A

1. Establish CHRFAM7A genotype and protein expression.
2. Determine pentamer stoichiometry and subunit orientation experimentally or evaluate a complete enumerated topology ensemble.
3. Confirm membrane assembly and surface expression before interpreting ligand docking.

### KIT, GABRA2, CACNA1D, and PDE9A

- KIT: obtain an experimentally constrained structure or ensemble for the 690–761 kinase-insert region before redocking KIT-223.
- GABRA2: reconcile transcript usage direction and verify whether the 73-aa product is translated.
- CACNA1D: use electrophysiology or conformational-state modeling to test C-terminal regulatory effects beyond the retained pocket.
- PDE9A: identify the exact alternate transcript/protein and verify that it changes the catalytic domain before comparative docking.

## 14. Presentation-ready evidence

The most defensible presentation sequence is:

1. FYN canonical validation and alternate kinase-pocket loss.
2. BACE1 canonical validation followed by the matched deletion-isoform comparison.
3. The all-candidate status figure showing why numerical results were not forced for every alternate.
4. CHRNA7/CHRFAM7A as a clearly labeled topology-sensitivity hypothesis.
5. GABRA2 and CACNA1D as structural exclusion/negative-control examples.
6. GNINA consensus as supporting evidence for BACE1 and the CHRFAM7A B-face hypothesis, with FYN disagreement shown rather than hidden.

Recommended figures:

- `../figures/expanded_campaign/matched_comparative_docking.png`
- `../figures/expanded_campaign/all_candidate_status.png`
- `../figures/expanded_campaign/canonical_crossdock_stability.png`
- `../figures/expanded_campaign/3D_BACE1_variant_pocket_overlay.png`
- `../figures/expanded_campaign/3D_CHRNA7_canonical_encenicline.png`
- `../figures/expanded_campaign/3D_CHRFAM7A_B_face_encenicline.png`
- `../figures/expanded_campaign/3D_CHRNA7_topology_site_overlay.png`
- `../figures/gnina_comparison/gnina_pose_selection_validation.png`
- `../figures/gnina_comparison/gnina_matched_comparison.png`
- `../figures/gnina_comparison/vina_gnina_rank_agreement.png`

Machine-readable values and per-run provenance are available in `../analysis/aggregate/`, `../analysis/gnina/`, and `../systems/`.

# Differential Transcript Usage, Isoform Structure, and Failed Neurodegenerative-Disease Drug Response

## 1. Executive summary

This project asks whether disease-associated changes in transcript usage can produce protein isoforms with altered drug response. The working model is that a drug may have been adequately designed for the canonical protein but may be less effective in a disease state because an alternative isoform becomes more abundant.

The analysis was organized around three structural mechanisms:

- **Mechanism A:** the canonical drug-binding region is absent from the alternative protein. The alternative cannot form the canonical drug site, so alternate docking is unnecessary.
- **Mechanism B:** the isoform difference intersects the canonical drug-binding region. The sequence change removes or changes drug-contact residues, potentially producing a different binding pose or weaker recognition.
- **Mechanism C:** the isoform difference is outside the direct binding region but may change allostery, gating, trafficking, conformational state, or binding mechanics. Static docking may show little difference even when physiological drug response changes.

The completed work identified:

1. **A: FYN–saracatinib** as a strong missing-pocket example.
2. **B: BACE1-202–verubecestat** as the strongest direct pocket-disruption example. The observed 401-aa isoform deletes canonical residues 21–120, including catalytic Asp93 and 12 of 28 experimentally defined drug-contact residues. Canonical redocking succeeds below 1 Å RMSD, whereas the alternate fails to recover a near-native pose in five independent runs.
3. **C: CACNA1D-214–isradipine** as the strongest distal-regulation example. The AD-associated isoform retains the mapped isradipine-contact region but lacks a distal C-terminal region associated with channel regulation. Published CaV1.3 studies provide a plausible link between distal splicing, gating, and dihydropyridine sensitivity.

The BACE1 structural result is strong but its transcript-level increase does not pass the available statistical significance gate. CACNA1D-214 has the stronger disease-associated transcript evidence, but its mechanism is necessarily a dynamic/state-dependent hypothesis rather than a simple docking-score difference.

These results are hypothesis-generating explanations for possible clinical failure mechanisms. They do not prove that an isoform caused a trial failure.

## 2. Scientific question

The central question is:

> Can differential transcript usage in a neurodegenerative disease produce protein isoforms that alter the binding, recognition, or physiological action of a drug that was developed against the canonical protein?

The more specific translational question is:

> For drugs that failed neurodegenerative-disease trials, can we identify at least one canonical/alternative protein pair where the canonical protein binds the drug plausibly, while the disease-associated alternative either lacks the binding site, changes the binding site, or changes the mechanics governing drug action?

The desired outcome was not simply a list of low docking scores. A useful example needed to satisfy three independent requirements:

1. **Clinical relevance:** the drug must have a documented failed or negative neurodegenerative-disease trial.
2. **Transcript relevance:** the alternative isoform must be present in the project’s differential-transcript-usage results, preferably with disease-associated usage gain.
3. **Structural/mechanistic relevance:** the canonical/alternative sequence difference must be mapped against an experimental binding site or a defensible regulatory region.

## 3. Biological rationale

Transcript usage changes can alter protein products in ways that are invisible to gene-level expression analysis. A gene may have unchanged total abundance while the relative abundance of its protein isoforms changes. Depending on the splice event, an alternative isoform can:

- delete a ligand-binding domain;
- remove individual contact residues;
- change the shape or electrostatics of a pocket;
- alter the catalytic machinery;
- modify an oligomeric interface;
- change trafficking or localization;
- change conformational-state occupancy;
- alter allosteric coupling between a distal regulatory region and the ligand site.

For this reason, binding affinity alone is not sufficient. A receptor can receive a favorable docking score in a non-native pose, particularly when a pocket is partially preserved or when the alternative structure is predicted. The workflow therefore prioritizes experimental pose recovery, sequence-to-contact mapping, structural confidence, and orthogonal rescoring.

## 4. Candidate definitions

### 4.1 Mechanism A: binding region absent

The canonical structure is used to define the drug-contact region. If the alternative sequence lacks the complete region or the receptor architecture required to form that region, the result is reported as structural absence. Docking an arbitrary ligand into an unrelated residual surface would not answer the biological question.

### 4.2 Mechanism B: binding-region change

The alternative sequence is aligned to the canonical sequence. Canonical drug-contact residues are defined from an experimental complex, generally using a heavy-atom distance cutoff. The candidate qualifies for B when the isoform difference removes or changes residues inside that contact shell. Comparative docking is then meaningful if the remaining pocket is modeled with adequate local confidence.

### 4.3 Mechanism C: distal or mechanical change

The direct binding region is retained, but a distal deletion, insertion, or alternative terminal region may alter receptor dynamics or functional coupling. A similar static pocket does not rule out altered physiological drug response. For ion channels and other state-dependent targets, electrophysiology, molecular dynamics, or conformational-state modeling is more informative than a single rigid-receptor score.

## 5. Input data and candidate screening

The initial candidate list was read from `outputs/DS_docking_candidate_pairs.md`. Candidate rows were screened against:

- disease-associated transcript usage;
- drug trial outcome;
- availability of a canonical experimental structure;
- availability or constructability of an alternative protein sequence;
- whether the direct drug-binding site could be mapped;
- whether the alternative receptor architecture was biologically defined.

The expanded screening considered additional combinations where the original row was structurally incomplete or not a clean B/C example. This was necessary because a candidate with an interesting score is not automatically a valid mechanistic example.

## 6. Computational methods

### 6.1 Canonical structural anchoring

Experimental protein–ligand complexes were downloaded and retained as immutable input files. Canonical binding sites were defined from the deposited ligand coordinates rather than from an arbitrary docking box. The principal structures were:

- FYN–saracatinib: PDB 10DJ;
- BACE1–verubecestat: PDB 5HU1;
- CaV1-family dihydropyridine structural framework: PDB 8E59 for the human channel coordinate frame and related published CaV structural evidence.

For BACE1, 5HU1 contains human BACE1 bound to verubecestat at high resolution and was used for both contact definition and coordinate alignment.

### 6.2 Sequence and contact mapping

Protein sequences were aligned with residue-level mapping. For each candidate, the analysis recorded:

- canonical and alternate lengths;
- inserted, deleted, and substituted regions;
- canonical residues contacting the experimental ligand;
- which contacts were retained, changed, or absent;
- whether catalytic or interface residues were affected;
- residue mapping into alternate numbering.

For BACE1, a 6 Å heavy-atom contact shell around verubecestat identified 28 canonical contact residues. For CACNA1D, 14 pocket seed residues were mapped into the alternative model.

### 6.3 Alternative-structure prediction

The BACE1-202 401-aa alternative was modeled using a five-member no-template ColabFold/AlphaFold2-ptm ensemble. Models were independently ranked and aligned to the canonical 5HU1 frame. Structural QC used:

- global aligned Cα RMSD;
- local retained-contact Cα displacement;
- local pLDDT;
- pTM;
- number of mapped contact residues;
- agreement across independent model seeds.

CACNA1D-214 had an existing five-model structural ensemble. The analysis used the same local pocket mapping and confidence checks.

### 6.4 AutoDock Vina

AutoDock Vina was used for conformational search and scoring. Matched BACE1 canonical/alternate runs used:

- five seeds: 1103, 2207, 3301, 4409, 5519;
- exhaustiveness: 32;
- 20 poses per seed;
- identical ligand preparation;
- identical search center and box;
- identical receptor-preparation convention;
- serial execution;
- maximum 16 CPU cores.

Canonical redocking was performed first. An alternate result was not interpreted as a binding change unless the canonical protocol recovered the known experimental pose.

### 6.5 RMSD and pose-recovery QC

For exact cognate redocking, heavy-atom RMSD was calculated in the fixed experimental receptor frame. A top-ranked RMSD below 2 Å was treated as successful pose recovery. RMSD was interpreted separately from affinity score because a favorable score in a displaced pose does not demonstrate preservation of the canonical binding mode.

### 6.6 GNINA orthogonal rescoring

GNINA 1.3.3 was used after Vina in `score_only` mode. It rescored the same Vina pose coordinates using CNN-based metrics without performing a second coordinate-optimization search. Coordinate-preservation QC confirmed that the rescoring step did not change ligand coordinates.

### 6.7 Structural visualization

Presentation figures were generated as high-resolution PNG, SVG, and PDF outputs where appropriate. PyMOL was used for 3D structural overlays and ligand-pose figures. Matplotlib was used for quantitative comparisons and mechanism diagrams. Figures show:

- canonical versus alternate structures;
- deleted or retained regions;
- experimentally defined contact residues;
- crystal ligand and alternate docked poses;
- per-seed docking variability;
- Vina pose recovery and GNINA confidence.

## 7. Results

### 7.1 A example: FYN–saracatinib

FYN–saracatinib is the strongest completed mechanism-A example. Saracatinib was tested in a Phase 2a Alzheimer’s disease study that did not provide the expected clinical benefit. The canonical FYN kinase structure contains a validated saracatinib-binding site, and exact redocking recovered the experimental ligand pose in nine independent runs, with median RMSD approximately 1.90 Å.

The AD-enriched alternative transcript `transcript160449.chr6.nic` showed approximately 58.0% usage in the AD group versus 2.15% in controls. The predicted alternative product is only 115 aa and lacks the entire kinase pocket. Therefore, alternate docking is not required: the canonical saracatinib-binding architecture is absent.

This example supports the strongest possible structural statement: the canonical drug target region is not present in the alternative product. It does not, by itself, establish that the alternative protein is translated, correctly localized, or responsible for the trial outcome.

### 7.2 B example: BACE1-202–verubecestat

#### Transcript and clinical context

Verubecestat failed the Phase 3 EPOCH Alzheimer’s trial for futility ([Egan et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29719179/)). The prodromal APECS trial was also negative ([Egan et al., 2019](https://www.nejm.org/doi/full/10.1056/NEJMoa1812840)).

The candidate transcript ENST00000392937 corresponds to UniProt isoform P56817-5, a 401-aa BACE1 product. In the available DTU table, usage increased in oligodendrocytes from 8.54% in controls to 18.50% in AD, a +9.97 percentage-point change. However, adjusted p = 0.655 and empirical FDR = 0.977. This is therefore a structural B example with limited transcript-level statistical support.

#### Sequence-to-pocket result

The alternative replaces canonical residues 1–20 and deletes canonical residues 21–120. The 5HU1 verubecestat contact shell contains 28 residues. The alternative lacks 12 of them:

- residues 70–75;
- residues 91–96.

Catalytic Asp93 is absent, while catalytic Asp289 remains and maps to alternate residue 189. Thus, the isoform change directly intersects the binding region and removes a catalytic residue.

#### Model QC

Across five predicted models:

- mean model pLDDT: approximately 85.9;
- pTM: 0.800;
- retained contacts mapped: 16/16 in every model;
- retained-contact mean pLDDT: approximately 87.7;
- mean retained-contact Cα displacement: 0.91 Å;
- mean global aligned Cα RMSD: 2.45 Å.

The remaining pocket is sufficiently modeled for hypothesis-generating comparative docking. The missing contacts are sequence-defined absences and cannot be corrected by model refinement.

#### Matched Vina and GNINA results

| Metric | Canonical BACE1 | BACE1-202 | Meaning |
|---|---:|---:|---|
| Vina affinity | −9.505 ± 0.019 kcal/mol | −9.374 ± 0.393 kcal/mol | Similar raw scores can be misleading |
| Top-pose RMSD | 0.957 ± 0.005 Å | 10.607 ± 0.630 Å | Alternate top poses are displaced |
| Best RMSD among 20 poses | 0.957 ± 0.005 Å | 5.245 ± 0.037 Å | No alternate run reaches near-native recovery |
| GNINA CNNscore | 0.920 ± 0.008 | 0.225 ± 0.220 | Alternate poses are less plausible |
| GNINA CNNaffinity | 7.367 ± 0.014 | 5.517 ± 0.194 | Orthogonal score favors canonical |

Canonical redocking succeeded in every seed below 1 Å RMSD. BACE1-202 did not produce a pose below 2 Å in any of the 100 sampled alternate poses. The main conclusion is a change in binding geometry, not simply a weaker numerical score.

#### B interpretation

The deletion removes a substantial portion of the experimentally observed verubecestat site, including Asp93. Vina can still score poses in the residual structure, but those poses are displaced from the canonical binding mode and are rejected by GNINA more strongly. This is direct computational support for mechanism B.

The correct wording is:

> BACE1-202 is a structurally compelling mechanism-B candidate: it removes catalytic Asp93 and 12/28 verubecestat-contact residues, and the alternative fails canonical pose recovery. Its transcript usage increase is observed in the AD dataset but does not pass the current statistical significance gate.

### 7.3 C example: CACNA1D-214–isradipine

#### Transcript and clinical context

Isradipine failed to slow progression in the STEADY-PD III Parkinson’s disease trial ([Parkinson Study Group, 2020](https://pubmed.ncbi.nlm.nih.gov/32227247/)). This is a neurodegenerative-disease trial failure, although not an Alzheimer’s disease trial.

CACNA1D-214 usage increased in AD inhibitory neurons from 1.81% to 28.11%, a +26.30 percentage-point change. The adjusted p-value was 4.12×10⁻³⁶ and empirical FDR was 0.0259, making this a statistically supported AD-associated isoform shift.

#### Sequence-to-pocket result

The canonical CACNA1D alpha-1 subunit is 2,161 aa, while CACNA1D-214 is 1,625 aa. The alternate contains a 20-aa insertion near canonical residue 492 and lacks canonical residues 1606–2161. The mapped isradipine-pocket residues—1078, 1081, 1082, 1085, 1154, 1156, 1194, 1198, 1205, 1209, 1212, 1489, 1492, and 1493—are retained and sequence-identical.

This is a clear distal-change pattern: the static binding region is retained while a distant regulatory C-terminal region is altered.

#### Structural QC

Across five alternate models:

- all 14 pocket residues mapped;
- mean pocket pLDDT: 83.96;
- minimum pocket pLDDT: approximately 72.3;
- local pocket Cα RMSD: 1.486 Å, range 1.437–1.522 Å;
- local shared-heavy-atom RMSD: 1.620 Å.

The pocket is structurally similar and sufficiently confident for a retained-pocket conclusion. The existing isradipine comparison that used duplicated receptor coordinates was excluded as evidence for an isoform-specific docking effect.

#### C interpretation

Dihydropyridine action depends on channel state and gating. Published CaV1.3 studies connect C-terminal splice variation to gating behavior and dihydropyridine sensitivity ([Huang et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23924992/); [Bock et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3234967/); [Ortner et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28592699/)).

The working mechanism is therefore:

1. CACNA1D-214 becomes more abundant in AD inhibitory neurons.
2. The direct isradipine-contact region remains present.
3. The distal C-terminal splice change alters gating, state occupancy, trafficking, or coupling.
4. The fraction of channels occupying drug-sensitive states may change.
5. A drug can retain a geometrically valid pocket while producing a different physiological response.

This is a mechanism-C hypothesis. Static docking cannot determine whether the channel spends more time in an isradipine-sensitive or insensitive state. Electrophysiology, state-resolved simulation, or targeted functional assays are needed for validation.

The correct wording is:

> CACNA1D-214 is a strong mechanism-C candidate: its AD-associated usage increase is statistically supported, the static isradipine pocket is retained, and published CaV1.3 biology provides a plausible route from distal splicing to altered state-dependent drug action.

### 7.4 Candidates not promoted

#### SORT1-209–latozinemab

SORT1-209 showed a strong AD excitatory-neuron usage gain and deletes canonical residues 1–137. However, the documented latozinemab epitope lies within residues 207–231 and is retained with identical sequence in all five alternate models. The local epitope RMSD was approximately 0.35 Å with high confidence. This is not a clean B example. It remains a possible trafficking or maturation-related C hypothesis because the deleted region includes the signal peptide/propeptide and part of the mature beta-propeller.

#### GFRA1-209–liatermin

The GFRA1 alternative deletes five residues outside the experimentally mapped GDNF interface. Interface displacement averaged only 0.26 Å, and confidence near the deletion was poor. It was not promoted as a convincing B or C example.

#### BACE1-476 and BACE1-457

These earlier BACE1 models are useful structural method-development comparisons and showed reproducible pose displacement. However, the exact isoforms were not confirmed as the AD-increased transcript in the available DTU table. They should not replace the evidence-qualified BACE1-202 result in the presentation.

## 8. Overall findings

### Finding 1: Docking score alone is insufficient

BACE1-202 and canonical BACE1 had similar mean Vina affinity values, but their pose geometries were radically different. This demonstrates why a successful canonical redocking control, pose RMSD, and orthogonal rescoring are essential.

### Finding 2: Mechanism B has the most direct structural link to altered binding

BACE1-202 removes experimentally defined contact residues and catalytic Asp93. Its alternate docking poses fail to recover the canonical geometry. This is the most direct computational explanation of how an isoform change could alter drug recognition.

### Finding 3: Mechanism C requires a functional-dynamics interpretation

CACNA1D-214 retains the static pocket, so a simple static docking comparison would miss the likely mechanism. The relevant variables are channel gating, conformational state, trafficking, and state-dependent ligand access.

### Finding 4: Transcript evidence and structural evidence must be reported separately

BACE1-202 has strong structural evidence but weak DTU significance. CACNA1D-214 has strong DTU evidence and plausible structural biology but a less direct docking endpoint. Combining these into one undifferentiated “best candidate” score would obscure important uncertainty.

### Finding 5: Candidate rejection is scientifically informative

SORT1 demonstrates that a large deletion does not automatically imply a binding-site deletion: the antibody epitope is retained. GFRA1 demonstrates that a small deletion does not automatically imply interface disruption. These exclusions prevent overinterpretation.

## 9. Limitations

1. Alternative protein structures are predictions, not experimentally determined structures.
2. Rigid-receptor docking does not reproduce protein dynamics, solvent rearrangement, membrane effects, or conformational-state populations.
3. Vina and GNINA scores are not experimental binding free energies or clinical efficacy measures.
4. A transcript usage change does not prove protein translation, correct folding, localization, or disease-specific functional abundance.
5. The BACE1-202 transcript increase is not statistically significant in the current table.
6. The CACNA1D-214 drug was tested in Parkinson disease, whereas its transcript enrichment evidence comes from AD; this is a cross-disease mechanistic hypothesis, not direct trial-subgroup evidence.
7. Trial failure is multifactorial and may reflect target biology, timing, safety, pharmacokinetics, patient selection, or inadequate target engagement rather than isoform switching.

## 10. Recommended validation experiments

### BACE1-202

- Confirm ENST00000392937 transcript usage with independent long-read or junction-specific RNA sequencing.
- Verify translation and cellular localization of the 401-aa product.
- Express canonical BACE1 and BACE1-202 in matched systems.
- Measure verubecestat binding and enzymatic inhibition directly.
- Test whether catalytic activity is lost because Asp93 is absent.
- Determine an experimental alternate structure or perform flexible/ensemble docking after biochemical confirmation.

### CACNA1D-214

- Confirm disease-associated isoform expression in relevant neuronal populations.
- Express canonical and CACNA1D-214 channels in the same cellular background.
- Compare voltage dependence, activation/inactivation, recovery, and trafficking.
- Measure isradipine inhibition across channel states rather than at one holding condition.
- Use molecular dynamics or enhanced-sampling methods to test coupling between the distal C terminus and the dihydropyridine pocket.

## 11. Presentation recommendations

The clearest presentation sequence is:

1. Introduce the hypothesis: disease-associated isoform switching can change drug response without changing total gene expression.
2. Show FYN–saracatinib as the simple A case where the canonical pocket is absent.
3. Show BACE1-202–verubecestat as the direct B case: contact residues and Asp93 are deleted, and the alternate fails pose recovery.
4. Show CACNA1D-214–isradipine as the C case: the pocket remains, but distal splicing may alter channel state dependence.
5. End with the limitations and proposed experimental tests.

Avoid saying that these analyses “explain why the drug failed” as a proven causal statement. Use “provide a mechanistic hypothesis for reduced or altered drug action in an isoform-shifted disease state.”

## 12. Output files

- Detailed candidate assessment: [BC_CANDIDATE_ASSESSMENT.md](BC_CANDIDATE_ASSESSMENT.md)
- Quantitative BACE1 metrics: `../analysis/bc_candidates/BACE1_202_401/matched_docking_metrics.csv`
- BACE1 docking plot: `../figures/bc_candidates/B_BACE1_matched_docking.png`
- BACE1 mechanism diagram: `../figures/bc_candidates/B_BACE1_isoform_mechanism.png`
- BACE1 3D overlay: `../figures/bc_candidates/B_BACE1_isoform_overlay.png`
- BACE1 pose displacement: `../figures/bc_candidates/B_BACE1_pose_displacement.png`
- CACNA1D mechanism diagram: `../figures/bc_candidates/C_CACNA1D_isoform_mechanism.png`
- CACNA1D full overlay: `../figures/bc_candidates/C_CACNA1D_full_structure_overlay.png`
- CACNA1D retained-pocket overlay: `../figures/bc_candidates/C_CACNA1D_retained_pocket_overlay.png`

All new docking and rescoring processes used a maximum of 16 CPU cores and completed successfully.

# Preliminary protein–drug docking analysis: detailed process

**Analysis date:** 24–25 August 2026  
**Primary engine:** AutoDock Vina 1.2.7  
**Status:** Preliminary, presentation-ready computational analysis  
**Main validated example:** canonical FYN–saracatinib redocking

> **Update:** This historical preliminary report is superseded by `EXPANDED_DOCKING_PROCESS.md` for the all-candidate campaign. Its CaV1.3 section has been corrected because the original staging used `3PE`, a lipid, rather than amiodarone (`BBI`).

## 1. Purpose and analysis strategy

The immediate objective was to obtain at least one scientifically defensible protein–drug docking example while BIOVIA Discovery Studio was unavailable. The candidate list contained several canonical-versus-alternate-protein comparisons, but not every pair was suitable for numerical docking. The workflow therefore began with structural triage rather than automatically assigning a docking score to every candidate.

The analysis followed this sequence:

1. Review each proposed protein, alternate transcript or isoform, and drug pair.
2. Determine whether the alternate product contains the relevant drug-binding domain.
3. Prefer an experimental protein–ligand structure when available.
4. Validate the docking protocol by redocking a crystallographic ligand.
5. Repeat docking with independent random seeds to measure reproducibility.
6. Retain only results supported by the structure and validation checks.
7. Report structural limitations instead of producing unsupported mutant or isoform scores.

This distinction is essential: a reproducible docking score is not automatically a validated binding prediction. A cognate redocking experiment additionally asks whether the method can recover a ligand pose already observed experimentally.

## 2. Candidate triage

The starting list was `outputs/DS_docking_candidate_pairs.md`. Candidates were evaluated for structural coverage, availability of an appropriate experimental template, and whether the proposed alternate protein retained the binding pocket.

| Candidate | Decision | Reason |
|---|---|---|
| FYN–saracatinib | Highest priority; perform canonical redocking | PDB 10DJ contains human FYN with crystallographic saracatinib, allowing direct pose-recovery validation. |
| Alternate FYN product | Do not calculate a pocket score | The predicted product ends at residue 115 and lacks complete SH3, SH2, and kinase domains. The validated kinase pocket is therefore absent. |
| KIT–masitinib | Retain as a canonical-pocket baseline | PDB 1T46 provides an experimental c-KIT ATP-pocket template, although its ligand is imatinib rather than masitinib. |
| KIT-223 | Do not report a comparative score | The deletion is at canonical residue 715, inside a kinase-insert segment unresolved in PDB 1T46 (residues 690–761). |
| CACNA1D control | Attempt a crystallographic-ligand control, then treat as exploratory | Corrected amiodarone docking recovered a near-native pose only at rank 7, not as the top-ranked pose. |
| GABRA2, BACE1, CHRNA7/CHRFAM7A, PDE9A | Not advanced in this handoff | The required alternate transcript, genotype, expression evidence, or structurally appropriate model was not available in the current analysis inputs. |

## 3. Software environment

Docking and analysis were run in the isolated `pocket_dock` Conda environment.

| Component | Version | Role |
|---|---:|---|
| Python | 3.10.20 | Workflow execution |
| AutoDock Vina | 1.2.7 | Docking search and Vina scoring |
| Meeko | 0.7.1 | Ligand/PDBQT preparation support |
| RDKit | 2025.09.5 | Ligand chemistry, 3D conformer generation, and SDF handling |
| Open Babel | 3.1.1 | Molecular-format and topology conversion |
| Biopython | 1.87 | PDB parsing, chain/ligand extraction, and structural alignment |
| NumPy | 2.2.6 | Coordinate and RMSD calculations |
| PyMOL open source | isolated `docking_viz` environment | Ray-traced molecular rendering |

GNINA was staged as a possible orthogonal rescoring method, but the available CUDA binary required `libcudnn.so.9`, which was not present. Consequently, all reported numerical results come from AutoDock Vina; no GNINA score is mixed into the analysis.

## 4. General docking configuration

The stored receptor and ligand inputs were converted to PDBQT format. The receptor was treated as rigid and the ligand as flexible according to the rotatable-bond definitions in its PDBQT file. AutoDock Vina used its standard `vina` scoring function.

For every reported run:

- The search box was centered on an experimentally observed ligand pocket.
- A fixed random seed made each run reproducible.
- Ten output poses were requested.
- Up to 16 Vina CPU cores were assigned to an individual docking job.
- The top-ranked pose was the pose with the most favorable Vina score within that run.
- Raw poses were stored as `docked_poses.pdbqt`.
- Run metadata, scores, and RMSD values were stored as `result.json`.

The later KIT replicate jobs were run serially so that no more than 16 docking cores were active in total.

## 5. FYN–saracatinib redocking validation

### 5.1 Experimental structure

Human FYN–saracatinib structure PDB **10DJ** was used. Chain A was selected as the receptor. The crystallographic ligand is residue **H8H A601**.

The staging script performed the following operations:

1. Parsed PDB 10DJ with Biopython.
2. Wrote chain A amino-acid coordinates as the protein receptor.
3. Removed waters, ligands, ions, and other non-protein residues from the receptor file.
4. Extracted saracatinib H8H A601 separately.
5. Recorded the 38 crystallographic ligand heavy atoms.
6. Calculated the ligand centroid and coordinate extent for grid definition.

The prepared receptor and ligand are:

- `FYN_saracatinib/prepared/fyn_chain_A.pdbqt`
- `FYN_saracatinib/prepared/saracatinib_H8H_A601.pdbqt`

### 5.2 Search space

The crystallographic ligand centroid defined the grid center:

```text
center = (-11.255, 14.853, -9.445) Å
box    = (20, 20, 26) Å
```

The underlying ligand extent was 8.155 × 8.711 × 14.715 Å. The final box was rounded and expanded from the extent-plus-padding calculation so that the complete ligand and surrounding pocket could be searched without clipping.

### 5.3 Independent docking runs

Nine independently seeded runs were included:

| Seed | Exhaustiveness | CPU limit | Output poses |
|---:|---:|---:|---:|
| 20260824 | 16 | 16 | 10 |
| 1103 | 32 | 16 | 10 |
| 2207 | 32 | 16 | 10 |
| 3301 | 32 | 16 | 10 |
| 4409 | 32 | 16 | 10 |
| 5519 | 64 | 16 | 10 |
| 6619 | 64 | 16 | 10 |
| 7723 | 64 | 16 | 10 |
| 8831 | 64 | 16 | 10 |

Using several seeds tests whether the result is robust to Vina’s stochastic search. Exhaustiveness values of 32 and 64 provide a more thorough search than Vina’s usual default of 8. The initial exhaustiveness-16 run was retained because it was completed with the same receptor, ligand, and grid and also passed the validation criterion.

### 5.4 Pose-recovery RMSD

The docking runner parsed heavy-atom coordinates from the prepared ligand and every docked pose. Because both remained in the same fixed experimental receptor coordinate frame, RMSD was calculated directly:

```text
RMSD = sqrt[(1/N) × Σ ||docked_atom_i − crystal_atom_i||²]
```

Hydrogen atoms were excluded. The ligand atom order preserved during PDBQT preparation was used to pair atoms. No additional receptor or ligand superposition was performed.

The validation criterion was:

```text
top-pose RMSD < 2.0 Å  →  successful pose recovery
```

Important limitation: this is an atom-order-based, fixed-frame heavy-atom RMSD and is **not symmetry-corrected**. Chemically equivalent symmetric atoms are not permuted to minimize RMSD. The value is therefore suitable as a transparent preliminary pose-recovery metric but should not be presented as a fully symmetry-aware benchmark.

### 5.5 Aggregation

All poses from all completed runs were combined into `FYN_saracatinib/analysis/all_poses.csv`. The top-scoring pose from each run was then used to calculate:

- number and proportion below 2 Å;
- mean and median top-pose RMSD;
- best-recovered top pose;
- the associated Vina scores.

The machine-readable summary is `FYN_saracatinib/analysis/summary.json`.

## 6. Structural interpretation of the alternate FYN product

The candidate transcript predicts a protein ending after residue 115. This was compared with the canonical 537-aa FYN domain architecture.

The key decision was structural, not numerical: a product ending at residue 115 does not contain a complete SH3 domain and entirely lacks the SH2 and kinase domains. Since saracatinib binds the kinase ATP pocket, there is no homologous pocket in the predicted short product into which the drug can be meaningfully docked.

For this reason, the alternate product was represented with a domain/pocket-loss schematic rather than an artificial docking score. This avoids treating an absent target site as if it were a weakly binding version of an intact site.

## 7. Canonical KIT–masitinib baseline

### 7.1 Experimental pocket template

PDB **1T46**, an imatinib-bound c-KIT kinase-domain structure, was selected. Chain A protein coordinates formed the receptor. The crystallographic imatinib ligand, STI A3, defined the ATP-pocket location.

The imatinib centroid and the final docking box were:

```text
center = (26.181, 26.061, 40.364) Å
box    = (22, 20, 28) Å
```

The final box was a rounded expansion of the ligand-derived 20.000 × 17.560 × 27.226 Å recommendation.

### 7.2 Masitinib preparation

Masitinib was generated from the stored SMILES string:

```text
Cc1ccc(cc1Nc2nc(cs2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C
```

RDKit ETKDGv3 generated a three-dimensional conformer with random seed 20260824. The conformer was optimized using the MMFF force field, saved as SDF, and converted to PDBQT for docking.

### 7.3 Replicate docking and summary

The same nine seeds and exhaustiveness schedule used for FYN were applied to canonical c-KIT. Ten poses were requested per run. The eight later replicates were executed serially under the aggregate 16-core limit.

Only the top Vina score from each run was used to summarize repeatability. Results were written to:

- `KIT_masitinib/analysis/masitinib_1T46_replicates.csv`
- `KIT_masitinib/analysis/masitinib_1T46_summary.json`

Unlike FYN–saracatinib, this calculation is **not** a cognate redocking validation because PDB 1T46 contains imatinib, not masitinib. There is no experimental masitinib pose in this structure against which to calculate a meaningful pose-recovery RMSD. Therefore, the KIT result measures numerical reproducibility within one prepared canonical pocket, not experimental accuracy.

### 7.4 Why the KIT-223 score was excluded

Early exploratory models of canonical KIT and KIT-223 were aligned to 1T46. The canonical and alternate models had alignment RMSDs of 1.728 Å and 1.332 Å, respectively, over 297 matched Cα atoms. However, those whole-model alignment values do not validate the local geometry around residue 715.

The decisive problem is that PDB 1T46 does not resolve residues 690–761. The KIT-223 deletion at residue 715 is located inside this missing segment. A comparison would therefore be determined by an unvalidated computationally rebuilt loop rather than experimental coordinates. Exploratory predicted-model docking outputs were excluded from the final quantitative interpretation, and no KIT-223 score is reported.

## 8. CACNA1D/CaV1.3 control — corrected

The human CaV1.3 structure PDB **8E59** contains amiodarone as residue **BBI A2201**. The original staging mistakenly extracted `3PE`, which is a phosphatidylethanolamine lipid. The staging script and control were corrected; the old `3PE` calculation remains only as provenance.

Corrected configuration:

```text
center = (151.334, 167.442, 149.793) Å
box    = (18.851, 27.189, 19.189) Å
seed   = 20260825
exhaustiveness = 32
```

The corrected top pose scored −1.041 kcal/mol with 14.95 Å RMSD. A near-native pose was found at rank 7 with 1.56 Å RMSD. Because the near-native pose was not ranked first, this is not a successful top-pose validation. CaV1.3/isradipine results were retained only as exploratory shared-pocket controls and not as validated affinity predictions.

## 9. Quantitative interpretation rules

The following rules were applied consistently:

- A FYN top pose below 2 Å was counted as successful crystallographic pose recovery.
- Vina scores were treated as model scores in kcal/mol, not measured binding free energies.
- Score repeatability across seeds was interpreted as numerical stability, not proof of accuracy.
- Different proteins were not ranked directly by raw Vina score.
- No score was assigned when the relevant pocket or altered residue lacked a defensible structural model.
- Failed controls were excluded from biological conclusions but retained in the audit trail.
- Docking alone was not used to claim cellular efficacy, disease causality, or an isoform-dependent drug response.

## 10. Visualization process

Docked PDBQT poses were converted to SDF while preserving proper ligand bond topology. PyMOL rendered the actual protein structures and docking coordinates as ray-traced 3000 × 2400 transparent PNG files. Protein backbones are shown as blue-gray cartoons; crystallographic saracatinib is teal, docked saracatinib is orange, and docked masitinib is purple.

Quantitative plots and domain/coverage schematics were generated directly from the CSV and JSON outputs with Matplotlib. Each plot is available as PNG, PDF, and SVG. No AI-generated protein artwork was used.

## 11. Reproducibility and audit trail

Key scripts:

- `02_SURVEYOR/master_surveyor/stage_fyn_saracatinib_redocking.py`
- `02_SURVEYOR/master_surveyor/run_vina_redock.py`
- `02_SURVEYOR/master_surveyor/launch_fyn_vina_replicates.sh`
- `02_SURVEYOR/master_surveyor/summarize_fyn_vina_redocking.py`
- `02_SURVEYOR/master_surveyor/stage_kit_masitinib_comparison.py`
- `02_SURVEYOR/master_surveyor/stage_kit_imatinib_redocking.py`
- `02_SURVEYOR/master_surveyor/launch_kit_masitinib_template_replicates.sh`
- `02_SURVEYOR/master_surveyor/summarize_kit_masitinib_template.py`
- `02_SURVEYOR/master_surveyor/create_standalone_docking_figures.py`
- `02_SURVEYOR/master_surveyor/render_standalone_docking_3d.pml`

Each run directory contains the docked PDBQT poses and a JSON record of the receptor, ligand, grid center, grid dimensions, random seed, CPU allocation, exhaustiveness, requested pose count, scores, and RMSD values. This makes the reported summaries traceable back to individual docking runs.

## 12. Recommended next steps

1. Confirm the alternate FYN transcript sequence, translation, and biological abundance before making a functional claim.
2. Build the KIT kinase-insert loop with a method that explicitly models residue 715, then validate local geometry and model uncertainty before redocking KIT-223.
3. Add a true experimental KIT–masitinib structure or an orthogonal validated protocol if one becomes available.
4. Repair the CUDA/cuDNN environment and use GNINA only as orthogonal rescoring; retain FYN cognate redocking as the primary pose-validation benchmark.
5. Consider symmetry-corrected ligand RMSD and interaction-fingerprint analysis for a publication-grade follow-up.

## 13. Method references

- Eberhardt J. et al. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *Journal of Chemical Information and Modeling* (2021). DOI: 10.1021/acs.jcim.1c00203.
- Forli S. et al. Computational protein–ligand docking and virtual drug screening with the AutoDock suite. *Nature Protocols* (2016). DOI: 10.1038/nprot.2016.051.

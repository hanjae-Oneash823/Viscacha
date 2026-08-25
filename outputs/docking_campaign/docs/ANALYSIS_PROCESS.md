# AutoDock Vina all-candidate campaign: detailed technical analysis process

**Analysis period:** 24–25 August 2026
**Primary docking engine:** AutoDock Vina 1.2.7
**Scope:** all eight rows in `outputs/DS_docking_candidate_pairs.md`
**Compute constraint:** maximum aggregate allocation of 16 CPU cores
**Document status:** current authoritative process report

## 1. Analytical objective

The objective was to test every protein–drug proposal in the candidate list while avoiding numerical comparisons that were not structurally or biologically defined. The candidate list did not represent one uniform docking problem. It contained exact cognate redocking controls, conventional canonical-versus-isoform comparisons, cross-docking into homologous pockets, truncations that remove the receptor or binding domain, and candidates for which the alternate protein sequence was not specified.

The campaign therefore used a decision-gated workflow. Docking was performed only after establishing that a receptor model contained a defensible binding site and that the comparison had a defined molecular interpretation. When an alternate protein lacked the relevant pocket, the endpoint was a structural exclusion rather than an artificial low-affinity score. When an alternate was undefined, the endpoint was an explicitly incomplete comparison rather than substitution of an arbitrary isoform.

## 2. Evidence classes and decision rules

Every result was assigned to one of five evidence classes before interpretation.

| Evidence class | Definition | Permitted interpretation |
|---|---|---|
| Exact cognate redocking | The deposited receptor and its co-crystallized ligand were separated, prepared, and docked back into the crystallographic site. | Tests whether the preparation, search box, scoring function, and search protocol recover a known pose. |
| Matched comparative docking | Canonical and alternate receptors used the same ligand preparation, coordinate frame, grid, Vina settings, and random seeds. | Supports a relative, model-dependent structural hypothesis within that protein–drug system. |
| Cross-docking | A drug without a cognate pose in the chosen structure was docked into a pocket defined by another ligand or a homologous receptor. | Exploratory pose/score result; no crystallographic pose-recovery claim. |
| Structural exclusion | The alternate product lacks the domain, interface, or assembly required to form the drug-binding site. | Report pocket/receptor absence; do not assign a Vina score. |
| Unresolved comparison | The alternate sequence, stoichiometry, genotype, expression evidence, or local structure was not sufficiently specified. | Report what is missing and do not infer an alternate affinity. |

The central validation criterion was a top-ranked heavy-atom pose RMSD below 2.0 Å for exact cognate redocking. Numerical repeatability across random seeds was treated separately from pose accuracy. A low score standard deviation indicates that Vina repeatedly converged on a similar score under the chosen setup; it does not demonstrate experimental accuracy.

## 3. Compute policy and software environment

### 3.1 CPU and GPU controls

- Every new Vina job was configured with `cpu=16`.
- Jobs were launched serially; two 16-thread Vina calculations were never allowed to overlap.
- ColabFold processes were pinned to logical CPUs 0–15 with `taskset -c 0-15`.
- `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `NUMEXPR_NUM_THREADS=1` prevented hidden BLAS thread pools from exceeding the aggregate limit.
- AlphaFold inference used one RTX 4090 GPU, while its host-side CPU affinity remained restricted to CPUs 0–15.

### 3.2 Recorded software versions

| Component | Version | Function |
|---|---:|---|
| Python | 3.10.20 | Workflow execution and aggregation |
| AutoDock Vina Python API | 1.2.7 | Docking search and Vina scoring |
| Meeko | 0.7.1 | Ligand/receptor PDBQT preparation support |
| RDKit | 2025.09.5 | SMILES parsing, 3D embedding, and force-field optimization |
| Open Babel | 3.1.x | Hydrogen addition, bond perception, and coordinate-format conversion |
| Biopython | 1.87 | PDB/mmCIF parsing, sequence mapping, and structural superposition |
| NumPy | 2.2.6 | Coordinate calculations and RMSD |
| Matplotlib | 3.10.9 | Quantitative figures |
| ColabFold | 1.6.2 | BACE1 deletion-isoform structure prediction |
| PyMOL open source | project `docking_viz` environment | Ray-traced structural figures |

GNINA was staged as a possible orthogonal CNN-based rescoring method. The available CUDA executable required `libcudnn.so.9`, which was not installed. No GNINA value was included in any table or conclusion; all reported numerical docking results are AutoDock Vina results.

## 4. Input provenance and structural triage

Public structures and reviewed sequences were downloaded into system-specific `inputs/` directories. For the expanded-stage downloads, source URLs, file sizes, SHA-256 hashes, and local paths are stored in `../analysis/metadata/expanded_inputs_manifest.json`. The earlier FYN, KIT, and CaV1.3 inputs predate that combined manifest and are traced through their retained source files and system-specific staging metadata. Raw structures were retained unchanged so every prepared receptor can be traced to its source.

### 4.1 Experimental structures

| System | Structure and ligand | Analytical role |
|---|---|---|
| FYN–saracatinib | PDB 10DJ, chain A, `H8H A601` | Exact cognate control and kinase-site definition |
| KIT–masitinib | PDB 1T46, chain A, imatinib `STI A3` | Experimental KIT ATP-site template for masitinib cross-docking |
| BACE1–verubecestat | PDB 5HU1, chain A, `66F A501` | Exact cognate control and model-alignment frame |
| CHRNA7–encenicline | PDB 7EKP, chains A–E, `I33 A601` | Exact cognate control and intersubunit-site definition |
| CHRFAM7A fusion | PDB 9QTO, chain A extracellular-domain construct | Experimental fusion-domain conformation for hybrid hypotheses |
| GABRA2 receptor | PDB 9CSB, native human β3–α1–β2–α2–γ2 pentamer | Target assembly containing adjacent α2 and γ2 subunits |
| GABA-A site donor | PDB 6X3X, diazepam `DZP D404` | Donor geometry for the α/γ benzodiazepine interface |
| CaV1.3 | PDB 8E59, chain A, amiodarone `BBI A2201` | Human channel frame and cognate control |
| PDE9A | PDB 4GH6, chain A, `LUO A601` | Catalytic-pocket definition for BI 409306 cross-docking |

### 4.2 Row-by-row triage decisions

| Candidate row | Structural finding before docking | Final analysis plan |
|---:|---|---|
| 1. FYN | The proposed 115-aa product has no complete SH3 domain and lacks the SH2 and kinase domains. | Validate canonical FYN; classify the alternate kinase pocket as absent. |
| 2. KIT | Ser715 lies within residues 690–761, which are unresolved in PDB 1T46. | Retain an experimental canonical baseline; attempt local model refitting, then reject the alternate result if steric quality fails. |
| 3. GABRA2 | The proposed 73-aa GABRA2-206 product cannot provide the extracellular ligand-binding domain or transmembrane architecture needed for a pentamer. | Build a native canonical α2/γ2 site; classify the alternate receptor as absent. |
| 4. CACNA1D | The 8E59 model ends at residue 1589, before the CACNA1D-214 product ends at residue 1625. | Use identical resolved coordinates as a shared-pocket negative control. |
| 5. BACE1-476 | A reviewed 476-aa deletion isoform alters the catalytic ectodomain. | Predict the alternate structure, validate its local pocket, and perform matched docking. |
| 6. BACE1-457 | A reviewed 457-aa deletion isoform removes a larger segment, including five residues in the experimental contact shell. | Predict and locally validate the model, then perform matched docking. |
| 7. CHRNA7/CHRFAM7A | Fusion-subunit stoichiometry and orientation were not specified. | Validate canonical α7 and evaluate two explicit one-fusion-subunit interface hypotheses. |
| 8. PDE9A | No coding-altered alternate transcript or sequence was named. | Perform canonical cross-docking only; do not invent an alternate. |

## 5. Receptor preparation

Protein chains were extracted from the experimental structures with Biopython. Water, ions, glycans, crystallographic ligands, lipids, and other heteroatoms were excluded from the rigid receptor unless explicitly required to define a reference ligand. Protein coordinates were then converted to PDBQT with polar hydrogens and AutoDock-compatible atom types/partial charges. This protein-only policy also removed the crystallographic Zn²⁺ and Mg²⁺ ions from 4GH6; the resulting PDE9A calculation was therefore classified as exploratory and not as a fully reconstructed catalytic-site model.

Meeko-prepared receptors were used where the conversion completed consistently. For the BACE1 and CHRNA7 matched comparisons, Open Babel-prepared PDBQT receptors were used for all members of a comparison. This ensured that canonical and alternate receptors shared one preparation convention. Crucially, the corresponding experimental canonical receptor independently passed exact redocking after that preparation choice, so the comparative protocol was not accepted solely because an alternate score was obtained.

All receptors were rigid during Vina docking. This approximation means side-chain, loop, and backbone relaxation induced by ligand binding was not modeled. Its consequences are most important for predicted alternate structures and for KIT, where the altered segment is intrinsically unresolved.

## 6. Ligand preparation

### 6.1 Cognate ligands

Saracatinib, verubecestat, encenicline, and corrected amiodarone were extracted from their deposited coordinates. Hydrogens and bond orders were assigned during SDF/PDBQT conversion, and Meeko generated flexible-ligand torsion trees and Gasteiger-type charges. The heavy-atom order and experimental coordinate frame were preserved so docked poses could be compared directly with the input crystallographic ligand.

### 6.2 Cross-docked ligands

Masitinib, AZD7325, isradipine, and BI 409306 were generated from explicit SMILES. RDKit ETKDGv3 embedded a three-dimensional conformer with random seed 20260825. MMFF94 optimization was used where parameters were available, followed by PDBQT conversion. Because these ligands were not initialized from a cognate pose in the selected receptor, their output poses do not have a meaningful fixed-frame redocking RMSD.

The exact stored SMILES and conformer metadata for AZD7325, isradipine, and BI 409306 are recorded in `../analysis/metadata/cross_docking_stage_metadata.json`; the masitinib SMILES is recorded in the KIT staging metadata.

## 7. Binding-site and search-box definition

For a cognate structure, the search center was the arithmetic centroid of crystallographic ligand heavy atoms. The initial box dimensions were the ligand coordinate extent plus 12 Å total padding for the expanded controls; the earlier FYN/KIT boxes were rounded expansions of the same ligand-centered principle. For homologous-site transfer, the donor ligand was transformed into the target receptor frame and its heavy-atom centroid/extent defined the box.

| System | Grid center (Å) | Box dimensions (Å) | Definition |
|---|---|---|---|
| FYN | (−11.255, 14.853, −9.445) | 20.000 × 20.000 × 26.000 | 10DJ saracatinib `H8H` |
| KIT | (26.181, 26.061, 40.364) | 22.000 × 20.000 × 28.000 for reported baseline | 1T46 imatinib `STI` |
| BACE1 | (24.887, 10.422, 21.858) | 23.129 × 16.370 × 18.460 | 5HU1 verubecestat `66F` |
| CHRNA7 | (142.818, 138.419, 90.348) | 22.027 × 19.752 × 16.264 | 7EKP encenicline `I33` |
| GABRA2 | (143.009, 97.405, 140.440) | 19.879 × 15.894 × 19.575 | 6X3X diazepam transferred to 9CSB |
| CACNA1D | (151.334, 167.442, 149.793) | 18.851 × 27.189 × 19.189 | corrected 8E59 amiodarone `BBI` |
| PDE9A | (77.260, 51.202, 39.870) | 19.790 × 18.581 × 23.684 | 4GH6 ligand `LUO` |

Each canonical/alternate matched comparison used exactly the same center and dimensions. Therefore, a difference was not caused by shifting or resizing the search region between receptor variants.

## 8. Alternate-structure construction and validation

### 8.1 BACE1 deletion isoforms

Reviewed UniProt sequences were used: P56817-1 (501 aa), P56817-2 (476 aa; deletion corresponding to canonical residues 190–214), and P56817-3 (457 aa; deletion corresponding to residues 145–188).

Each alternate was modeled independently with ColabFold 1.6.2 using MMseqs2 UniRef plus environmental MSAs, `alphafold2_ptm` model 1, one model, one seed (20260825), three recycles, no dropout, no template input, and no Amber relaxation. Model-level confidence was not treated as sufficient by itself. Models were sequence-aligned to experimental 5HU1 chain A using a local pairwise alignment with match score 2, mismatch −1, gap-open −8, and gap-extension −0.5. Matched Cα atoms defined the rigid superposition into the experimental coordinate frame.

The experimental pocket-contact set comprised protein residues with any atom within 6.0 Å of a verubecestat heavy atom. Local Cα RMSD and pLDDT were calculated only for mapped contact residues.

| Model | Matched Cα globally | Global fit RMSD | Mean model pLDDT | pTM | Mapped contact Cα | Pocket RMSD | Mean/minimum pocket pLDDT |
|---|---:|---:|---:|---:|---:|---:|---:|
| BACE1-476 | 363 | 5.499 Å | 77.75 | 0.722 | 28/28 | 1.484 Å | 86.87 / 73.19 |
| BACE1-457 | 344 | 2.816 Å | 78.70 | 0.750 | 23/28 | 1.937 Å | 85.27 / 69.81 |

BACE1-457 mapped only 23 of 28 contacts because its deletion removes five residues from the 5HU1 contact shell. The local pocket metrics were considered adequate for hypothesis-generating matched docking, but the receptors remain predictions and were not promoted to experimental structures.

### 8.2 KIT-223 local-refit attempt

Initial full-domain predicted models were aligned to 1T46, but whole-domain RMSD could not validate the local Ser715 environment because residues 690–761 are absent from the experimental template. A second attempt selected template residues with any atom within 12 Å of the imatinib centroid, fitted the predicted models using 66–67 local Cα pairs, and retained protein residues within 22 Å of the center for docking.

The canonical and KIT-223 refits gave local Cα RMSDs of 0.692 and 0.690 Å, respectively, with mean local pLDDT values of 86.31 and 85.34. The diagnostic refit runs used the same center as the reported baseline and the unrounded 20.000 × 17.560 × 27.226 Å imatinib-derived box. However, minimum local pLDDT values were only 43.62 and 41.41, and docking produced strongly positive top scores (+0.183 and +16.529 kcal/mol). Visual/energetic inspection indicated severe steric incompatibility rather than a credible binding solution. These runs were retained for audit but rejected from biological interpretation.

### 8.3 GABRA2 interface transfer

The first candidate pentamer, 9CXC, was rejected because its α2 and γ2 chains were not adjacent at the intended benzodiazepine interface. No result from that geometry was retained.

For the accepted geometry, 6X3X α1 chain D was sequence-aligned and superposed onto 9CSB α2 chain D. Interface transfers used local pairwise alignment with match score 2, mismatch −1, gap-open −6, and gap-extension −0.5. The fit used 331 Cα pairs, 87.01% identity over the aligned region, and produced 0.985 Å Cα RMSD. After applying the same transformation, the neighboring γ2 chains overlaid at 0.740 Å Cα RMSD. This two-subunit check confirmed that transferred diazepam occupied a geometrically corresponding α2/γ2 interface rather than a rotationally mismatched site.

### 8.4 CHRNA7/CHRFAM7A hypotheses

The canonical encenicline site in 7EKP is between chains A and B. The experimental 9QTO fusion extracellular-domain chain was aligned separately to each site-forming chain. Both alignments used 206 identical matched Cα positions and produced a 1.147 Å fit RMSD. Two receptors were constructed: fusion replacing chain A and fusion replacing chain B. Interface checks found no atom pairs closer than 1.5 Å after replacement.

These models are intentionally topology hypotheses. They do not establish the number of fusion subunits in a biological pentamer, which face of a ligand site is altered, membrane assembly, surface expression, or sample genotype.

## 9. Vina execution and RMSD calculation

The common runner instantiated `Vina(sf_name="vina", cpu=16, seed=<seed>)`, loaded one rigid receptor and one flexible ligand, computed maps for the recorded box, and called `dock(exhaustiveness=32, n_poses=20)` for expanded-campaign jobs. Poses and energies were written to `docked_poses.pdbqt` and `result.json` in each run directory.

For cognate ligands, hydrogen atoms were excluded and fixed-frame RMSD was calculated from preserved PDBQT atom order:

```text
RMSD = sqrt[(1/N) × Σ ||x_docked,i − x_crystal,i||²]
```

No receptor or ligand superposition was performed after docking because the input and output already occupied the same receptor coordinate frame. This metric is atom-order based and not symmetry corrected; chemically equivalent symmetric atoms were not permuted to minimize RMSD. It is therefore transparent and reproducible but less rigorous than a symmetry-aware benchmark.

For cross-docked ligands generated from SMILES, `--skip-rmsd` was used. Any historical RMSD field generated against an arbitrary non-cognate starting conformer was ignored.

## 10. Replication and statistical aggregation

Expanded comparisons used pilot seed 20260825 plus five independent seeds: 1103, 2207, 3301, 4409, and 5519. Every accepted group therefore contains six top-ranked observations. The BACE1 and CHRNA7 comparisons were paired by seed. For each seed, the alternate-minus-canonical score difference was calculated before computing the mean and sample standard deviation.

The earlier FYN and canonical KIT campaigns used nine reported runs. FYN included one exhaustiveness-16 run, four exhaustiveness-32 runs, and four exhaustiveness-64 runs. KIT used the same general nine-seed schedule for the experimental canonical pocket. These established results remain in their system-specific analyses rather than being duplicated in the 54-row expanded aggregate table.

For every group, the analysis script calculated the arithmetic mean, sample standard deviation (`n−1` denominator), minimum, maximum, and—where meaningful—mean top-pose RMSD and count below 2.0 Å. No null-hypothesis tests or p-values were applied: the six seeds measure stochastic search reproducibility, not independent biological replicates. The aggregate contains 54 accepted replicate records: nine six-run groups. Raw system directories contain 82 `result.json` files because they additionally retain pilots, alternate preparation checks, historical controls, and rejected KIT calculations.

## 11. System-specific execution record

### 11.1 FYN–saracatinib

Chain A and saracatinib `H8H A601` were extracted from 10DJ. The 38-heavy-atom ligand defined the ATP-site box. Nine seeded cognate redocking runs were performed. The alternate product was not docked because a 115-aa product cannot contain the kinase ATP pocket.

### 11.2 KIT–masitinib

The imatinib site in 1T46 defined the canonical cross-docking box. Masitinib was embedded and optimized from SMILES, then docked in nine independent canonical runs. Predicted canonical/KIT-223 whole-domain and local-refit comparisons were attempted. Because Ser715 is in an unresolved experimental segment and the refit docking showed severe clashes/positive scores, the comparative values were excluded.

### 11.3 GABRA2–AZD7325

A native 9CSB pentamer was constructed after rejecting 9CXC site geometry. The 6X3X diazepam site was transferred and validated across both α and γ subunits. AZD7325 was docked in six canonical cross-docking runs. No 73-aa alternate receptor was constructed because it cannot form the binding domain or pentamer.

### 11.4 CACNA1D–isradipine

The original control mistakenly extracted 8E59 residue `3PE`, a phosphatidylethanolamine lipid. PDB chemical annotation identified amiodarone as `BBI A2201`; staging, coordinates, ligand PDBQT, grid, wrapper script, and documentation were corrected. The old `3PE` result remains provenance only.

Corrected amiodarone redocking was evaluated over 20 poses. Isradipine was then cross-docked in the shared 8E59 pocket. Canonical and CACNA1D-214 pilot commands used the identical receptor file, ligand, grid, seed, and settings and therefore produced identical output, which was the expected negative-control outcome.

### 11.5 BACE1–verubecestat

Canonical Meeko and Open Babel receptor preparations both passed pilot redocking. The Open Babel preparation was selected for the matched canonical/alternate series because both predicted alternate receptors could be prepared consistently by the same route. Six paired runs were completed for canonical BACE1, BACE1-476, and BACE1-457.

### 11.6 CHRNA7/CHRFAM7A–encenicline

Canonical Meeko and Open Babel preparations both passed pilot redocking. The Open Babel series was selected for consistent preparation across canonical and hybrid pentamers. Six paired runs were performed for canonical α7, fusion-at-A-face, and fusion-at-B-face receptors.

### 11.7 PDE9A–BI 409306

PDE9A chain A was extracted from 4GH6. The `LUO A601` pocket defined the search box, and BI 409306 was cross-docked in six canonical runs. The protein-only receptor excluded 4GH6 Zn²⁺ and Mg²⁺, further limiting the score to exploratory use. No alternate model was created because the candidate row did not identify a coding-altered sequence.

## 12. Quality-control and exclusion logic

Results entered the main interpretation only if all applicable conditions were satisfied:

1. source structure or sequence was explicitly recorded;
2. the target site existed in the receptor product;
3. the search box was defined from an experimental ligand or validated structural transfer;
4. canonical cognate redocking passed when a cognate complex was available;
5. alternate models had interpretable local pocket confidence/geometry;
6. matched comparisons used the same ligand, grid, preparation convention, seeds, and Vina settings;
7. cross-docking was labeled exploratory and did not receive a fabricated RMSD;
8. steric clashes, undefined alternates, and absent pockets were preserved as exclusions;
9. scores were compared only within the same receptor–drug system;
10. docking was not used to claim binding free energy, efficacy, disease causality, brain exposure, or clinical benefit.

## 13. Visualization and reporting

Quantitative plots were generated directly from the aggregate CSV/JSON files. Points represent individual seeded runs; diamonds/error bars represent mean ± sample SD. The 2 Å line appears only where pose-recovery RMSD is meaningful.

PyMOL rendered molecular views from the actual prepared proteins and stored docking poses. Canonical and alternate pockets were superposed in the experimental coordinate frame before rendering. The images are structural renders, not generated illustrations.

## 14. Reproducibility map

### 14.1 Primary data and reports

- Aggregate summary: `../analysis/aggregate/expanded_summary.json`
- Accepted replicate table: `../analysis/aggregate/replicate_results.csv`
- Candidate disposition: `../analysis/aggregate/candidate_status.csv`
- Download provenance: `../analysis/metadata/expanded_inputs_manifest.json`
- Interface-transfer metadata: `../analysis/metadata/interface_model_metadata.json`
- Cross-docking ligand metadata: `../analysis/metadata/cross_docking_stage_metadata.json`
- System inputs, models, prepared files, logs, poses, and results: `../systems/`
- Presentation figures: `../figures/expanded_campaign/`

### 14.2 Core scripts

- `02_SURVEYOR/master_surveyor/stage_expanded_candidate_inputs.py`
- `02_SURVEYOR/master_surveyor/stage_exact_cognate_controls.py`
- `02_SURVEYOR/master_surveyor/stage_cross_docking_candidates.py`
- `02_SURVEYOR/master_surveyor/stage_interface_models.py`
- `02_SURVEYOR/master_surveyor/stage_bace1_alternate_models.py`
- `02_SURVEYOR/master_surveyor/refit_kit_models_to_pocket.py`
- `02_SURVEYOR/master_surveyor/run_vina_redock.py`
- `02_SURVEYOR/master_surveyor/run_expanded_docking_replicates.py`
- `02_SURVEYOR/master_surveyor/analyze_expanded_docking_campaign.py`
- `02_SURVEYOR/master_surveyor/render_expanded_docking_3d.pml`

## 15. Methodological limitations

Vina uses a simplified empirical scoring function and a rigid receptor. It omits explicit bulk-solvent thermodynamics, long-timescale conformational changes, membrane dynamics, protonation-state ensembles, induced fit, and cellular context. Membrane proteins were docked without an explicit lipid bilayer, and the protein-only PDE9A receptor omitted catalytic-site metals. Score differences of approximately 1 kcal/mol are therefore hypotheses, not measured free-energy changes.

The BACE1 alternates are AlphaFold-derived single models without relaxation or ensemble sampling. The CHRFAM7A calculations depend on assumed subunit topology. The GABRA2 and PDE9A calculations are cross-docking rather than cognate validation. The CaV1.3 protocol did not rank its near-native amiodarone pose first. KIT-223 lacks an experimentally anchored local structure. Finally, transcript detection does not by itself prove translation, stability, assembly, localization, or pharmacological relevance of the alternate protein.

These limitations define the appropriate use of the campaign: prioritization of structural hypotheses and presentation of one or more validated canonical docking examples, followed by targeted experimental validation.

## 16. Method references

- Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry* (2010). DOI: 10.1002/jcc.21334.
- Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. *Journal of Chemical Information and Modeling* (2021). DOI: 10.1021/acs.jcim.1c00203.
- Forli S et al. Computational protein–ligand docking and virtual drug screening with the AutoDock suite. *Nature Protocols* (2016). DOI: 10.1038/nprot.2016.051.
- Jumper J et al. Highly accurate protein structure prediction with AlphaFold. *Nature* (2021). DOI: 10.1038/s41586-021-03819-2.
- Mirdita M et al. ColabFold: making protein folding accessible to all. *Nature Methods* (2022). DOI: 10.1038/s41592-022-01488-1.

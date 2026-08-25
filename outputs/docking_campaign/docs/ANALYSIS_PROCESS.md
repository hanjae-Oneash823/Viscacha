# Expanded docking campaign: detailed analysis process

## Purpose and scope

This campaign attempted every row in `outputs/DS_docking_candidate_pairs.md`. The list contains three different scientific situations:

1. a conventional canonical-versus-alternate docking comparison;
2. an alternate that physically lacks the drug-binding pocket; or
3. an alternate that is not sufficiently defined by the candidate file.

These cases were not forced into the same numerical format. A pocket-absent result is a structural result, not a failed docking run. An unspecified alternate remains explicitly incomplete rather than being replaced with an arbitrary isoform.

## Compute policy

- Every new Vina calculation used `cpu=16`.
- Vina calculations were run serially, so two 16-core jobs never overlapped.
- AlphaFold/ColabFold jobs were pinned with `taskset -c 0-15`.
- BLAS and related thread pools were restricted with `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `NUMEXPR_NUM_THREADS=1`.
- GPU calculations used one RTX 4090 and remained restricted to the same 16 logical CPU IDs.

## Candidate triage before docking

| Row | Structural question | Numerical plan |
|---:|---|---|
| 1, FYN | The 115-aa alternate loses the kinase domain. | Validate canonical redocking; classify alternate as pocket absent. |
| 2, KIT | KIT-223 deletes Ser715, but this segment is unresolved in 1T46. | Retain experimental canonical baseline; attempt predicted models and reject them if the pocket is invalid. |
| 3, GABRA2 | GABRA2-206 is only 73 aa and cannot form a pentameric receptor. | Cross-dock AZD7325 to a native α2/γ2 pocket; classify alternate as receptor absent. |
| 4, CACNA1D | The alternate truncation occurs after the experimentally resolved drug pocket. | Run a shared-pocket negative control with identical canonical/alternate coordinates. |
| 5–6, BACE1 | Reviewed 476- and 457-aa deletion isoforms alter the catalytic ectodomain. | Predict both alternate structures, quality-check the local pocket, then compare with validated canonical redocking. |
| 7, CHRNA7/CHRFAM7A | Mixed-receptor stoichiometry and orientation are not specified. | Validate canonical redocking; test two explicit one-fusion-subunit interface hypotheses. |
| 8, PDE9A | The candidate file does not name the coding-altered transcript. | Run canonical BI 409306 cross-docking only; do not invent the alternate. |

## Experimental structures and pocket definitions

| System | Experimental structure | Role |
|---|---|---|
| FYN–saracatinib | PDB 10DJ | Exact cognate redocking control and ATP-site definition. |
| KIT–masitinib | PDB 1T46 with imatinib | Experimental kinase-pocket template for masitinib cross-docking. |
| BACE1–verubecestat | PDB 5HU1, ligand `66F` | Exact cognate redocking control and alternate-model alignment frame. |
| α7–encenicline | PDB 7EKP, ligand `I33` | Exact cognate redocking control and mixed-receptor site frame. |
| CHRFAM7A fusion | PDB 9QTO | Experimental fusion-subunit extracellular-domain conformation. |
| GABRA2 | PDB 9CSB | Native human β3–α1–β2–α2–γ2 pentamer. |
| GABA-A benzodiazepine site | PDB 6X3X, diazepam `DZP` | α/γ interface pocket transferred onto the 9CSB α2/γ2 interface. |
| CaV1.3 | PDB 8E59, amiodarone `BBI` | Human channel pocket and corrected redocking control. |
| PDE9A | PDB 4GH6, ligand `LUO` | Human catalytic-domain pocket for BI 409306 cross-docking. |

The downloaded file URLs, SHA-256 hashes, and local paths are recorded in `../analysis/metadata/expanded_inputs_manifest.json`.

## Important CaV1.3 input correction

The original CaV1.3 staging script incorrectly treated residue `3PE` in 8E59 as amiodarone. The PDB chemical annotation shows that:

- `BBI` is amiodarone;
- `3PE` is a phosphatidylethanolamine lipid.

The staging script was corrected to extract `BBI A2201`. The old failed `3PE` run is retained only as provenance and must not be presented as an amiodarone result.

## Ligand preparation

- Cognate redocking ligands were extracted in the experimental receptor frame.
- Open Babel added hydrogens and inferred ligand bonding from the deposited coordinates.
- Meeko 0.7.1 generated ligand PDBQT files with Gasteiger charges.
- AZD7325, isradipine, and BI 409306 were generated from explicit SMILES with RDKit ETKDGv3, seed 20260825, followed by MMFF94 optimization.
- Cross-docked ligands did not begin in the experimental frame, so fixed-frame RMSD was deliberately disabled for those runs.

## Receptor preparation

- Protein chains were extracted without crystallographic ligands, waters, glycans, or other heteroatoms.
- Meeko preparation was used where it produced an internally consistent receptor.
- For BACE1 and α7 matched comparisons, Open Babel preparation was used for every member of a comparison after its experimental canonical control independently passed redocking.
- Receptor preparation warnings and failed model attempts were retained. A method was accepted for comparative use only after its matched canonical redocking control passed.

## BACE1 alternate modeling and validation

Reviewed UniProt sequences were used:

- P56817-1: 501 aa canonical;
- P56817-2: 476 aa isoform, deleting canonical residues 190–214;
- P56817-3: 457 aa isoform, deleting canonical residues 145–188.

Each alternate was predicted independently with local ColabFold 1.6.2 using:

- MMseqs2 UniRef plus environmental MSA;
- AlphaFold2-ptm model 1;
- three recycles;
- one model and one seed;
- no Amber relaxation.

Model-level scores were pLDDT 77.8/pTM 0.722 for BACE1-476 and pLDDT 78.7/pTM 0.750 for BACE1-457. More importantly, local experimental-pocket checks gave:

| Alternate | Mapped 5HU1 contacts | Pocket Cα RMSD | Mean pocket pLDDT | Minimum pocket pLDDT |
|---|---:|---:|---:|---:|
| BACE1-476 | 28/28 | 1.484 Å | 86.87 | 73.19 |
| BACE1-457 | 23/28 | 1.937 Å | 85.27 | 69.81 |

BACE1-457 maps only 23 contacts because its deletion removes five residues that contact verubecestat in 5HU1.

## Interface-transfer checks

### GABRA2

The first native template considered, 9CXC, failed because its α2 and γ2 chains are not adjacent. No docking was retained from that geometry.

The accepted 9CSB transfer aligned 6X3X α1 onto 9CSB α2 with 0.985 Å Cα RMSD. After the same transform, the neighboring γ2 chains overlaid at 0.740 Å. This confirmed that the transferred diazepam box represented a real α2/γ2 interface rather than a rotationally mismatched pentamer site.

### CHRNA7/CHRFAM7A

The exact 7EKP encenicline site lies between chains A and B. Two one-fusion-subunit hypotheses were built:

- replace chain A with the aligned 9QTO fusion extracellular domain;
- replace chain B with the aligned 9QTO fusion extracellular domain.

The experimental extracellular-domain overlay was 1.147 Å Cα RMSD. Interface clash checks found no atom pairs below 1.5 Å. These are explicit topology hypotheses; they do not establish the biological stoichiometry of a sample.

## Vina settings and replication

- AutoDock Vina Python API 1.2.7;
- rigid receptor and flexible ligand;
- 16 CPU threads;
- exhaustiveness 32;
- up to 20 reported poses;
- pilot seed 20260825 plus replicate seeds 1103, 2207, 3301, 4409, and 5519;
- one common grid for every member of a matched comparison.

Fixed-frame heavy-atom RMSD was calculated only when the input ligand coordinates came from a cognate structure in the same receptor frame. It is atom-order based and not symmetry corrected.

## Quality criteria

1. Exact cognate redocking: top-pose RMSD below 2 Å.
2. Alternate model: local pocket confidence and structural overlay inspected before docking.
3. Matched comparison: identical ligand preparation, grid, Vina version, exhaustiveness, and seeds.
4. Cross-docking: score reported as exploratory; no invented RMSD.
5. Absolute scores were not compared across different proteins or drugs.
6. Model clashes, absent pockets, and undefined alternates were recorded as exclusions rather than converted into biological claims.

## Reproducible outputs

- Aggregate machine-readable summary: `../analysis/aggregate/expanded_summary.json`
- All replicate-level values: `../analysis/aggregate/replicate_results.csv`
- Per-row disposition: `../analysis/aggregate/candidate_status.csv`
- Quantitative figures: `../figures/expanded_campaign/`
- Per-run PDBQT poses and JSON: each candidate's `runs/` directory
- BACE1 model-quality report: `../systems/BACE1_verubecestat/prepared/alternate_model_quality.json`
- Interface geometry report: `../analysis/metadata/interface_model_metadata.json`

## Interpretation limits

Docking ranks poses under a simplified scoring function. It does not directly measure binding affinity, efficacy, brain exposure, or clinical benefit. Predicted alternate structures add uncertainty even when their local pockets have high pLDDT. The CHRFAM7A result additionally depends on an assumed interface topology. Finally, the local expression data must independently verify that BACE1 alternatives, GABRA2-206, and CHRFAM7A are relevant to the AD samples before a disease-specific conclusion is made.

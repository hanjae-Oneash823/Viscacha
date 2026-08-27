# Why these candidate pairs were selected for Discovery Studio

## Executive summary

This shortlist was designed to test a specific, falsifiable explanation for an unsuccessful neurotherapeutic trial: a drug can bind the reference (canonical) protein but bind an AD-enriched protein isoform less well, or the disease-associated isoform may no longer contain the drug-binding pocket. A useful result therefore requires three linked evidence streams:

1. A disease-associated isoform switch in the supplied AD/control transcript data.
2. A drug with a documented target and an unsuccessful or inconclusive human neuro-disease clinical program.
3. A structural change that could alter the relevant binding site, and a defensible model for testing it.

The candidates are deliberately not all equivalent. **FYN–saracatinib** is the strongest presentation candidate: the alternate transcript is strongly AD-enriched and loses the complete kinase domain used by saracatinib. **KIT–masitinib** is the strongest conventional two-structure docking comparison: the alternate is a known, very small splice change and can be compared in the same kinase model, although a large docking difference is not expected. **BACE1–verubecestat** is an established AD drug-failure reference experiment, but it must be retained only if the project's own data demonstrate an AD-enriched BACE1 coding isoform. CACNA1D, GABRA2, PDE9A, and CHRNA7/CHRFAM7A provide boundary cases and controls that prevent over-interpreting a single positive-looking docking score.

**Important interpretation rule:** docking is hypothesis-generating. It cannot prove why a clinical trial failed, because exposure, blood–brain-barrier delivery, target engagement, cell type, disease stage, and downstream biology also matter. The appropriate conclusion is “consistent with an isoform-dependent loss of target engagement,” not “explains the trial failure.”

## Selection method

### Dataset filters

The analysis first prioritized coding-alternate transcripts in the project candidate table, then required a large and directionally consistent AD/control usage difference, structural annotation, and a biologically relevant cell type. The most useful events have either a change in the drug pocket itself or loss of the domain containing that pocket. A strong transcript-level signal alone was not enough.

### Drug and structural filters

For each gene, the drug had to meet most of the following conditions:

- a clear, direct protein target;
- an available experimental structure, or a credible structure/model of the relevant domain;
- a reported unsuccessful/inconclusive neuro-disease clinical study or regulatory outcome;
- a small-molecule or defined ligand amenable to Discovery Studio docking; and
- an alternate form that is meaningful to compare structurally.

This explains why some otherwise strong AD isoform switches were not selected: they lacked a suitable neurotherapeutic trial, a direct ligand, or a change near a tractable pocket.

## Candidate-by-candidate rationale

| Tier | Pair | Why it was selected | What the structural experiment can answer |
|---|---|---|---|
| Primary | FYN canonical vs `transcript160449.chr6.nic`; saracatinib | Strong AD-specific switch plus complete loss of the Fyn kinase domain. | Whether saracatinib re-docks to canonical Fyn; the alternate should be classified as **pocket absent**, not assigned a misleading docking score. |
| Primary | KIT canonical vs KIT-223; masitinib | Strong switch and a known one-residue splice variant that preserves a directly comparable kinase structure. | Whether a Ser715 deletion produces any reproducible change in masitinib pose/score or local kinase-insert geometry. |
| Conditional reference | BACE1-501 vs BACE1-476/457; verubecestat | Very strong AD clinical-failure and structural precedent. | Whether an AD-enriched BACE1 alternate in this dataset changes the inhibitor pocket. |
| Structural boundary case | GABRA2 canonical vs GABRA2-206; AZD7325 | A strong switch that removes receptor architecture. | A clear canonical-only receptor-pocket figure; alternate docking is invalid because the functional receptor/pocket is absent. |
| Negative-control style comparison | CACNA1D canonical vs CACNA1D-214; isradipine | Strong switch but the long C-terminal truncation is expected to leave the dihydropyridine pocket intact. | Whether docking is robust to a distal isoform change; a null result is informative. |
| Conditional follow-up | PDE9A canonical vs a coding-altered AD-enriched PDE9A isoform; BI 409306 | A relevant clinical program and tractable catalytic pocket. | Only useful when the selected isoform modifies/truncates the catalytic domain; N-terminal-only variants should be excluded. |
| Conditional translational model | CHRNA7 homopentamer vs CHRNA7/CHRFAM7A mixed pentamer; encenicline | Directly models a human-specific altered receptor context linked to variable pharmacology. | Whether mixed receptor composition changes the ligand site; this is not a standard splice-isoform experiment. |

### 1. FYN — saracatinib (AZD0530): highest-priority presentation pair

**Dataset evidence.** In excitatory neurons, `transcript160449.chr6.nic` showed usage of **58.0% in AD versus 2.15% in controls** (29 versus 2 supporting counts; adjusted chi-square *P* = 3.56 × 10⁻¹²; permutation *P* = 0.0379 in the candidate table). It is annotated as novel and predicts a C-terminal truncation after residue 115. The canonical Fyn SH3, SH2, and kinase regions are consequently lost.

**Why the drug pairing is unusually clean.** Saracatinib is an ATP-competitive Src-family kinase inhibitor with Fyn as a key target. The human Fyn–saracatinib co-crystal structure (PDB **10DJ**) covers the kinase domain (residues 260–537), giving an unusually strong positive-control structure for re-docking. A phase II Alzheimer’s study of saracatinib did not show significant effects on its clinical/imaging outcomes. [Fyn–saracatinib structure](https://www.rcsb.org/structure/10DJ); [AD phase II trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC6646979/).

**Predicted outcome and proper framing.** Saracatinib should recover a credible pose in the canonical kinase pocket. The alternate predicted protein lacks the entire kinase domain; no homologous ATP pocket exists to dock. The result should be presented as “AD-enriched transcript associated with loss of the druggable Fyn kinase domain,” rather than as a numerical canonical-versus-alternate docking-score difference.

**Main risk.** This is a novel transcript with limited control counts. Confirm junction support, open reading frame, and ideally protein translation before treating it as a biological isoform rather than an RNA-level event. This is the strongest mechanistic visual, but not the strongest proof of a translated alternate protein.

### 2. KIT — masitinib: best direct two-isoform docking test

**Dataset evidence.** In inhibitory neurons, KIT-223 was used in **37.7% of AD versus 3.5% of control** cells (106 versus 26 counts; adjusted chi-square *P* = 3.87 × 10⁻⁴⁹; permutation *P* reported as 0). Its annotated protein change is a one-amino-acid deletion at Ser715 in the kinase-insert region.

**Why it is structurally practical.** KIT-223 corresponds to a naturally occurring Ser715-minus splice isoform, so it is a more conventional and credible protein comparison than a novel truncation. Masitinib directly inhibits KIT, and KIT kinase structures provide a well-defined ATP-pocket template (for example, imatinib-bound KIT, PDB **1T46**). Masitinib’s ALS regulatory submission was refused after the evidence did not establish a favorable benefit–risk balance; this supplies the required failed-neurotherapeutic context, while not proving the same mechanism applies in AD. [KIT Ser715 isoform](https://pmc.ncbi.nlm.nih.gov/articles/PMC1850714/); [KIT structure](https://www.rcsb.org/structure/1T46); [EMA assessment](https://www.ema.europa.eu/en/medicines/human/EPAR/masitinib-ab-science).

**Predicted outcome and proper framing.** The deletion is not a known core masitinib-contact residue, so a large score difference is not expected. That is a strength, not a weakness: this pair tests whether the docking workflow can identify a subtle local structural change without claiming that every AD-enriched splice change must disrupt binding. Report pose reproducibility, local conformational changes, and interaction fingerprints, not only a single score.

### 3. BACE1 — verubecestat: established AD reference, conditional on the dataset

**Why it remains on the list.** BACE1 is directly tied to amyloid biology, verubecestat had major negative Alzheimer’s trials, and canonical BACE1 has experimentally solved inhibitor-bound structures. BACE1-476 and BACE1-457 are shorter alternatives whose changes can affect regions relevant to catalytic function and inhibitor-pocket architecture. [BACE1 isoform review](https://pmc.ncbi.nlm.nih.gov/articles/PMC9286785/); [verubecestat trial](https://www.nejm.org/doi/full/10.1056/NEJMoa1812840).

**Why it is conditional.** BACE1 was chosen as a high-value external reference, not as a confirmed hit from the supplied Viscacha candidate table. It should move forward only after the project’s AD/control data show that one of these coding alternatives is genuinely AD-enriched in the relevant cells. Without that evidence, it is a useful validation benchmark but does not support the presentation’s disease-associated isoform-switch narrative.

**Experiment.** Re-dock verubecestat into an experimental canonical BACE1 structure, build each alternate with an explicit quality check, then compare whether the active-site geometry and key contacts are retained. Exclude a model if the alternate cannot support a credible folded catalytic domain.

### 4. GABRA2 — AZD7325 (BAER-101): receptor-loss boundary case

**Dataset evidence.** GABRA2-206 was detected in astrocytes at **9.3% in AD and 0% in controls**. The predicted alternate is extensively C-terminally truncated and loses much of the ligand-binding and transmembrane architecture required for a functional GABA-A receptor.

**Why it was kept.** AZD7325/BAER-101 is an alpha2/alpha3-selective GABA-A receptor modulator with an unsuccessful/limited neurodevelopmental clinical signal. It is a useful companion to FYN because it makes the same logical point—loss of a druggable architecture—at a receptor rather than a kinase.

**Critical modeling limitation.** GABA-A receptors are pentamers. Docking an isolated alpha2 subunit is not biologically adequate, and docking the truncated alternative as though it were a receptor would be invalid. Use a canonical pentamer model for the positive control and a domain/topology figure for the alternate. [BAER-101 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13183936/).

### 5. CACNA1D — isradipine: informative negative-control comparison

**Dataset evidence.** In inhibitory neurons, CACNA1D-214 was **28.1% in AD versus 1.8% in controls** (104 versus 5 counts; adjusted chi-square *P* = 4.1 × 10⁻³⁶; permutation *P* = 0.0259). The annotation predicts a 536-amino-acid C-terminal truncation.

**Why it was selected despite lower expectation.** Isradipine’s phase III Parkinson’s trial was negative. The channel’s dihydropyridine-binding site lies in the transmembrane pore-forming region, whereas this candidate’s change is C-terminal; the pocket is therefore expected to remain. This makes CACNA1D a valuable specificity test: an isoform can be strongly AD-enriched without necessarily changing binding. [STEADY-PD III trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC7465126/).

**Modeling caution.** CaV1.3 is a large membrane channel. Comparative structural inspection is preferable to overconfident rigid docking. If docking is attempted, use a membrane-aware channel structure/model, define the known dihydropyridine site, and phrase the result as a control rather than a trial-failure explanation.

### 6. PDE9A — BI 409306 (osoresnontrine): conditional pocket-change candidate

**Why it was selected.** PDE9A has a well-defined catalytic pocket suitable for small-molecule docking, and BI 409306 was studied in Alzheimer’s disease without a successful efficacy outcome. This is a structurally tractable way to ask whether isoform selection could change inhibitor engagement. [BI 409306 clinical study](https://pubmed.ncbi.nlm.nih.gov/30755255/).

**Why it is conditional.** PDE9A alternatives frequently differ in N-terminal regulatory/localization segments. Such variants should *not* be docked as a claimed binding-difference pair if the catalytic domain is unchanged. Select a PDE9A transcript only when the project data identify AD enrichment and its coding change overlaps, deletes, or destabilizes the catalytic domain.

### 7. CHRNA7/CHRFAM7A — encenicline: altered receptor-composition model

**Why it was selected.** CHRFAM7A is a human-specific fusion gene that can incorporate with CHRNA7-containing alpha7 nicotinic receptor complexes and alter pharmacology. It offers a clinically relevant example in which the drug target in humans is not a simple canonical homopentamer. Encenicline provides the failed cognitive-therapy context.

**Why it is conditional and separate.** This is a mixed-pentamer/composition question, not a canonical-versus-splice-isoform question. Include it only with sample-specific CHRFAM7A expression or genotype support, and model an explicitly defined mixed pentamer rather than an isolated chain. It is best used as an advanced follow-up, not the central figure. [CHRFAM7A pharmacology evidence](https://www.nature.com/articles/s41380-023-02389-1).

## Recommended Discovery Studio workflow and decision rules

1. **Verify the biological comparison first.** Confirm transcript identity, coding potential, and AD enrichment in the stated cell type. Use the exact canonical and alternate amino-acid sequences; do not substitute a similarly named isoform.
2. **Validate canonical docking.** Use a ligand-bound experimental structure when available. Re-dock the cognate ligand, measure pose recovery (RMSD), and retain the method only if it reproduces the known binding mode.
3. **Build and quality-control alternate models.** Perform sequence alignment, inspect the changed region, assess template coverage/confidence, and reject models whose pocket-bearing domain is absent or implausibly folded.
4. **Use the same protocol for valid pairs.** Keep protonation, binding-site definition, flexibility choices, and scoring settings matched between canonical and alternate structures. Generate multiple poses and inspect key interactions.
5. **Report the appropriate endpoint.** For a retained pocket, report score distribution, pose cluster, and contact changes. For a lost domain/pocket, report “not dockable: pocket absent” with a structural-domain map. Never treat an artificial score from a truncated chain as evidence of weak binding.
6. **Use controls.** FYN and GABRA2 are expected pocket-loss cases; KIT is a subtle direct comparison; CACNA1D is expected to retain its primary pocket. A defensible campaign needs all three behaviors.

## What can be claimed after the docking campaign

If canonical re-docking validates and a high-confidence AD-enriched alternate loses or materially changes the pocket, the conclusion can be: “The data support an isoform-dependent target-engagement hypothesis that is consistent with reduced drug binding in AD-associated cells.” It remains necessary to test protein expression, ligand binding/functional inhibition, and patient-level association before attributing a human trial outcome to the isoform.

## Final prioritization

For a Thursday presentation, use **FYN–saracatinib** as the main mechanistic slide and **KIT–masitinib** as the main side-by-side docking result. Add **CACNA1D–isradipine** as a planned negative-control test. Add BACE1 only if the project data confirm its alternate is AD-enriched; otherwise label it clearly as an external AD reference. Keep GABRA2 as a domain-loss visual, not a two-score docking comparison.

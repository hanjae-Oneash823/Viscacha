# Mechanism-B and Mechanism-C Candidate Assessment

## Executive decision

Two candidates are suitable for the presentation, but they have different evidence profiles.

| Class | Recommended pair | Decision | Main strength | Main limitation |
|---|---|---|---|---|
| B | BACE1-202 (ENST00000392937; P56817-5; 401 aa)–verubecestat | Use as a **structurally strong, transcript-evidence-limited B example** | The isoform deletes catalytic Asp93 and 12/28 crystallographic drug-contact residues; matched docking loses the experimental pose | Its increase in the available AD DTU table does not pass the statistical gate |
| C | CACNA1D-214–isradipine | Use as the **primary C example** | The pocket is retained, the distal splice change is statistically strong in AD, and published CaV1.3 work connects distal splice regulation to dihydropyridine sensitivity | No exact human CaV1.3–isradipine co-crystal is available; the gating mechanism remains a hypothesis |

These examples support mechanistic hypotheses. They do not establish that isoform switching caused either clinical-trial failure.

## B: BACE1-202–verubecestat

### Biological and trial context

Verubecestat is a small-molecule BACE1 inhibitor. The EPOCH randomized Phase 3 trial enrolled 1,958 participants with mild-to-moderate Alzheimer disease and was terminated for futility; the drug did not improve cognitive or functional outcomes and produced adverse events. The primary trial report is [Egan et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29719179/). The prodromal AD APECS trial was also negative ([Egan et al., 2019](https://www.nejm.org/doi/full/10.1056/NEJMoa1812840)).

The canonical structural anchor is human BACE1 bound to verubecestat in PDB [5HU1](https://www.rcsb.org/structure/5HU1), resolved at 1.50 Å. This eliminates ambiguity about the canonical binding geometry.

### Transcript evidence

The observed transcript is ENST00000392937, corresponding to UniProt isoform P56817-5 and a 401-aa protein. In oligodendrocytes, its mean usage increased from 8.54% in controls to 18.50% in AD, a change of +9.97 percentage points.

This row does **not** pass the available statistical evidence gate: gene-adjusted p = 0.655 and empirical FDR = 0.977. Therefore the presentation must call it “AD-observed” or “nominally increased,” not “significantly gained in AD.” This is the principal weakness of the B example.

### Sequence-to-pocket mapping

Pairwise sequence mapping shows that the alternate N terminus replaces canonical residues 1–20 and deletes canonical residues 21–120. Experimental contacts were defined directly from 5HU1 as canonical protein residues with any atom within 6.0 Å of verubecestat.

- 28 canonical contact residues were identified.
- 12/28 contacts are absent in the 401-aa isoform: residues 70–75 and 91–96.
- Catalytic Asp93 is absent.
- Catalytic Asp289 is retained and maps to alternate residue 189.
- 16/28 contact residues remain.

This satisfies the operational definition of mechanism B: the isoform difference intersects and removes part of the experimentally observed drug-binding region. It also removes half of the catalytic aspartate dyad, so the inferred biochemical consequences likely extend beyond altered inhibitor pose to loss or severe alteration of normal BACE1 catalytic function.

### Model quality and ensemble agreement

The 401-aa sequence was modeled independently in a five-member no-template ColabFold/AlphaFold2-ptm ensemble. Across the five ranked models:

- mean model pLDDT was approximately 85.9;
- pTM was 0.800;
- all 16 retained crystallographic contacts mapped in every model;
- mean pLDDT across retained contacts was 87.7;
- retained-contact Cα displacement from 5HU1 averaged 0.91 Å;
- global aligned Cα RMSD averaged 2.45 Å.

These values pass the local comparative-model gate for docking the remaining partial pocket. They cannot restore the 12 physically missing contacts; that absence comes from sequence mapping, not from uncertain model coordinates.

### Matched docking protocol

Canonical and alternate receptors used the same Open Babel receptor-preparation convention. Verubecestat coordinates and topology came from 5HU1. AutoDock Vina was run serially with five seeds (1103, 2207, 3301, 4409, and 5519), exhaustiveness 32, 20 poses per seed, and a hard maximum of 16 CPU cores. The search box was identical for canonical and alternate receptors.

The canonical protocol passed exact redocking before comparative interpretation: every canonical seed recovered the crystallographic pose below 1 Å. GNINA 1.3.3 then rescored the Vina poses in `score_only` mode, preserving pose coordinates and providing an orthogonal pose-quality assessment.

### Results

| Metric, top-ranked pose | Canonical mean ± SD | BACE1-202 mean ± SD | Interpretation |
|---|---:|---:|---|
| Vina affinity (kcal/mol) | −9.505 ± 0.019 | −9.374 ± 0.393 | Raw Vina score alone is misleadingly similar |
| RMSD to crystal pose (Å) | 0.957 ± 0.005 | 10.607 ± 0.630 | Alternate top poses occupy a different geometry |
| Best RMSD among 20 poses (Å) | 0.957 ± 0.005 | 5.245 ± 0.037 | No alternate run samples a near-native pose |
| GNINA CNNscore | 0.920 ± 0.008 | 0.225 ± 0.220 | Orthogonal model strongly disfavors most alternate poses |
| GNINA CNNaffinity | 7.367 ± 0.014 | 5.517 ± 0.194 | GNINA predicts weaker alternate-pose affinity |

The most important result is pose recovery, not Vina affinity. The canonical complex reproducibly returns to the crystallographic geometry, while the alternate never approaches it within the conventional 2 Å redocking threshold. The high-looking alternate Vina scores are therefore not evidence of preserved binding.

### Defensible interpretation

The 21–120 deletion removes a substantial part of the experimentally defined pocket, including Asp93. In the predicted alternate structure, Vina finds energetically scored poses elsewhere in the remaining cavity, but these are displaced from the canonical binding mode and receive much poorer GNINA scores. This is direct computational support for **different binding geometry caused by an isoform change within the binding region**.

The defensible presentation language is: “BACE1-202 is a structurally compelling mechanism-B example observed at higher mean usage in AD oligodendrocytes, but its DTU increase is not statistically supported in the present dataset.”

## C: CACNA1D-214–isradipine

### Biological and trial context

Isradipine is a dihydropyridine L-type calcium-channel blocker. In STEADY-PD III, 336 participants with early Parkinson disease were randomized; the adjusted treatment effect on the primary UPDRS outcome was −0.27 points with P = 0.85, providing no evidence of slowed progression ([Parkinson Study Group STEADY-PD III Investigators, 2020](https://pubmed.ncbi.nlm.nih.gov/32227247/)). This is a failed neurodegenerative-disease trial, although it is Parkinson disease rather than AD.

### Transcript evidence

CACNA1D-214 increased in AD inhibitory neurons from 1.81% to 28.11%, a change of +26.30 percentage points. The gene-adjusted p value was 4.12×10⁻³⁶ and the empirical FDR was 0.0259. This is a strong, statistically supported AD-associated isoform shift.

### Sequence and pocket mapping

The canonical CaV1.3 alpha-1 subunit is 2,161 aa; CACNA1D-214 is 1,625 aa. Exact sequence mapping identifies:

- a 20-aa insertion around canonical residue 492;
- loss of canonical residues 1606–2161 from the alternate C terminus;
- complete retention and sequence identity of all 14 mapped dihydropyridine-pocket seed residues used in this analysis: 1078, 1081, 1082, 1085, 1154, 1156, 1194, 1198, 1205, 1209, 1212, 1489, 1492, and 1493.

Thus, the dominant isoform difference is outside the static drug-contact region, satisfying the operational definition of mechanism C.

### Structural ensemble QC

Across five independently ranked alternate models:

- all 14 pocket seed residues mapped in every model;
- mean pocket pLDDT was 83.96;
- minimum pocket pLDDT was approximately 72.3;
- the 59-residue local pocket Cα RMSD averaged 1.486 Å (range 1.437–1.522 Å);
- local shared-heavy-atom RMSD averaged 1.620 Å.

The static pocket is therefore well supported and structurally similar. A large score difference obtained by docking to duplicated or effectively identical receptor coordinates would not constitute isoform evidence, so the earlier identical-coordinate CACNA1D comparison is explicitly excluded.

### Why a distal change can still matter

Dihydropyridines are state-dependent channel ligands: apparent potency depends on voltage-dependent gating and occupancy of favored channel states. Published CaV1.3 work shows that C-terminal splice variation changes channel gating and dihydropyridine sensitivity ([Huang et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23924992/); [Bock et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3234967/); [Ortner et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28592699/)). The structural drug-contact region itself is supported by CaV1-family dihydropyridine studies ([Cooper et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7502546/)).

The resulting mechanism-C hypothesis is that CACNA1D-214 preserves the local isradipine pocket but changes distal regulatory architecture, which may change channel gating, state occupancy, trafficking, or coupling and therefore alter effective drug action. Static rigid-receptor docking cannot test those dynamic quantities. Electrophysiology or state-resolved molecular simulation is the appropriate next validation layer.

### Defensible interpretation

The defensible presentation language is: “CACNA1D-214 is a strong mechanism-C candidate: its AD-associated usage increase is robust, the local dihydropyridine pocket remains structurally intact, and the distal splice change has a literature-supported route to altered state-dependent pharmacology.”

It should not be presented as proof that CACNA1D-214 caused isradipine’s Parkinson trial failure or that the trial enrolled an AD-defined isoform subgroup.

## Candidates investigated but not promoted

### SORT1-209–latozinemab

SORT1-209 has a strong AD excitatory-neuron gain and deletes canonical residues 1–137. However, patent-defined latozinemab epitope residues 207–231, including T218, Y222, S223, and S227, are retained and sequence-identical. Five-model epitope RMSD is only 0.35 Å. It is therefore not a B example. It remains a secondary C/trafficking hypothesis because the deletion removes the signal peptide/propeptide and part of the mature beta-propeller. Latozinemab’s Phase 3 FTD-GRN study failed its primary endpoint according to the [sponsor’s topline report](https://investors.alector.com/news-releases/news-release-details/alector-announces-topline-results-latozinemab-phase-3-trial).

### GFRA1-209–liatermin

The five-residue deletion lies outside the experimentally mapped GDNF interface; interface displacement averages only 0.26 Å, and the deletion-neighbor confidence is low. This is not a persuasive B candidate and is weaker than CACNA1D as C.

## Presentation recommendation

1. Lead with the established A example, FYN–saracatinib.
2. Present BACE1-202–verubecestat as the cleanest structural B mechanism, with the DTU non-significance visible on the figure rather than hidden.
3. Present CACNA1D-214–isradipine as the strongest C mechanism and explicitly explain why rigid docking is not the correct decisive test for a state-dependent ion channel.
4. Use SORT1 only as a backup C/trafficking example.

The presentation-ready vector figures are in `outputs/docking_campaign/figures/bc_candidates/`; the underlying per-seed measurements are in `outputs/docking_campaign/analysis/bc_candidates/BACE1_202_401/matched_docking_metrics.csv`.

## Relaxed SURVEYOR follow-up screen (completed)

To search beyond the trial-failure-only subset, we screened all SURVEYOR candidate groups for AD-enriched, protein-altered isoforms with drug annotations. The strongest newly testable C hypotheses were ERBB4-204–neratinib and PARP2-201–niraparib. Each was docked with AutoDock Vina using five matched seeds, exhaustiveness 32, and 16 CPU threads per process. The experimentally resolved kinase/catalytic-domain pocket was used for both isoforms.

| Candidate | Canonical mean ± SD (kcal/mol) | Alternate mean ± SD (kcal/mol) | Paired Δ | Interpretation |
|---|---:|---:|---:|---|
| ERBB4-204–neratinib | −10.374 ± 0.039 | −10.374 ± 0.039 | 0.000 | Distal 16-aa deletion; pocket retained in the modeled domain. Identical receptor coordinates mean docking is a pocket-retention control, not proof of equal whole-protein pharmacology. |
| PARP2-201–niraparib | −10.548 ± 0.018 | −10.548 ± 0.018 | 0.000 | N-terminal 13-aa insertion outside the catalytic pocket; same limitation as ERBB4. |

KIT-223–masitinib was also tested as a potential B case. The pocket-refit alternate gave a nominal score penalty of +16.30 ± 0.07 kcal/mol relative to the canonical receptor. This result fails structural QC: the altered residue is in an unresolved/poorly supported region of the 1T46 template and the alternate receptor contains preparation-induced clashes. It must not be presented as a validated affinity difference; at most it is a warning that an underdetermined isoform model can generate an artifact.

These runs therefore add biologically motivated C hypotheses but do not replace the stronger mechanistic examples already identified (FYN for pocket loss/B-like behavior and CACNA1D for a distal, state-dependent C mechanism). Structural references for the follow-up pockets are [ErbB4 PDB 3BBT](https://www.rcsb.org/structure/3BBT) and [PARP2 PDB 8HLQ](https://www.rcsb.org/structure/8HLQ).

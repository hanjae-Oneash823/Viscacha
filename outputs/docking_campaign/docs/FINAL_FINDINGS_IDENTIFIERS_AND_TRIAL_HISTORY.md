# Final A/B/C Findings: Isoform Identifiers, DTU Evidence, and Drug-Trial History

This document gives the exact identifiers and project-level DTU values for the three final examples, followed by the relevant clinical-development history. DTU values are from the project outputs; trial facts are linked to primary trial registries, trial publications, FDA records, or sponsor announcements.

## Summary table

| Class | Canonical protein / transcript | Alternate protein / transcript | Lengths | Isoform change | Drug and disease | Project DTU result |
|---|---|---|---:|---|---|---|
| A | FYN; UniProt **P06241-1**, Ensembl **ENST00000354650** | Novel `transcript160449.chr6.nic`; no stable UniProt accession | **537 aa → 115 aa** | C-terminal truncation / premature termination; SH3, SH2, and kinase regions absent | Saracatinib (AZD0530), Alzheimer disease | Control **2.15%**, AD **58.00%**, Δ **+55.85 percentage points**; adjusted p **3.56×10⁻¹²**, raw p **4.38×10⁻¹³**, permutation p **0.0379** |
| B | BACE1; UniProt **P56817-1**, Ensembl **ENST00000313005** | BACE1 isoform 5; UniProt **P56817-5**, Ensembl **ENST00000392937**; project label **BACE1-202** | **501 aa → 401 aa** | N-terminal alternative/replacement plus internal deletion of canonical residues 21–120; net 100-aa shortening, not a simple C-terminal truncation | Verubecestat (MK-8931/SCH 900931), Alzheimer disease | Control **8.54%**, AD **18.50%**, Δ **+9.97 percentage points**; transcript p **0.7508**, gene-adjusted p **0.6554**, empirical FDR **0.9771**; **not statistically significant** |
| C | CACNA1D/CaV1.3; UniProt **Q01668**, Ensembl **ENST00000350061** | Project label **CACNA1D-214**; Ensembl **ENST00000636627**; no curated UniProt accession assigned in the project | **2161 aa → 1625 aa** | Distal C-terminal truncation of canonical residues 1606–2161 plus a 20-aa insertion near residue 492 | Isradipine, early Parkinson disease | Control **1.805%**, AD **28.108%**, Δ **+26.303 percentage points**; adjusted p **4.115×10⁻³⁶**, raw p **7.737×10⁻³⁸**, permutation p **0.0259** |

## A — FYN / saracatinib

### Protein identifiers and structural change

The canonical reference is UniProt P06241-1, Ensembl ENST00000354650, a 537-aa Fyn tyrosine kinase. UniProt identifies P06241-1 as the canonical sequence and gives a length of 537 aa ([UniProt P06241](https://www.uniprot.org/uniprotkb/P06241/entry)).

The alternative is the project-specific novel transcript `transcript160449.chr6.nic`. It has no stable UniProt isoform accession. The translated product used in the project is predicted to be 115 aa long. This is a C-terminal truncation/premature-termination product: it does not contain the canonical SH3, SH2, or protein-kinase domains. The canonical kinase domain occupies the C-terminal part of the full protein, so the saracatinib-binding architecture is absent rather than merely altered.

### Project DTU evidence

The transcript was detected in excitatory neurons at:

- control usage: **2.1505%**;
- AD usage: **58.0000%**;
- absolute change: **+55.8495 percentage points**;
- raw DTU p-value: **4.3769×10⁻¹³**;
- adjusted DTU p-value: **3.5603×10⁻¹²**;
- permutation p-value: **0.0379**.

The alternate is therefore strongly AD-enriched in the project dataset, subject to the usual caveat that transcript detection does not prove stable translation or correct localization.

### Drug and development history

Saracatinib, also called AZD0530, is an oral Src-family kinase inhibitor with Fyn activity. It was originally developed by **AstraZeneca** for oncology; AstraZeneca supplied the compound for the Alzheimer studies. The Alzheimer program was an academic-industry repurposing effort rather than a conventional AstraZeneca-sponsored late-stage AD program.

The relevant Alzheimer development sequence was:

1. **Phase Ib, NCT01864655:** a 4-week safety, tolerability, pharmacokinetic, and CNS-availability study. The study was conducted at Yale and ran in 2013–2014; the publication reports the Yale study period as July 2013–March 2014 ([Phase Ib report](https://pmc.ncbi.nlm.nih.gov/articles/PMC4396171/)).
2. **Phase IIa, NCT02167256:** a randomized 52-week study in mild AD, with 159 participants, initiated in **December 2014**. ClinicalTrials.gov lists Yale University as sponsor and the study as Phase 2 ([registry](https://clinicaltrials.gov/study/NCT02167256)); the publication reports enrollment from December 23, 2014 through November 30, 2016 ([JAMA Neurology trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC6646979/)).

The maximum Alzheimer disease phase reached was **Phase IIa**. The trial did not demonstrate a statistically significant treatment effect on the primary cerebral glucose-metabolism endpoint or on clinical/biomarker secondary endpoints. The primary outcome difference in relative CMRgl decline was −0.006 units/year, p = 0.34. Gastrointestinal adverse events, especially diarrhea, were more frequent with saracatinib, and more participants discontinued treatment in the saracatinib group ([trial results](https://pmc.ncbi.nlm.nih.gov/articles/PMC6646979/)).

### Why it failed

The reported clinical reason was **lack of efficacy**, not a proven isoform-related mechanism: saracatinib did not slow cerebral metabolic decline or improve cognition/function in this Phase IIa study. The project’s A hypothesis is that, in an AD cell population with strong gain of the truncated FYN transcript, the canonical Fyn kinase drug target may not be present in the corresponding alternative protein product.

## B — BACE1-202 / verubecestat

### Protein identifiers and structural change

The canonical active beta-secretase is UniProt P56817-1, Ensembl ENST00000313005, 501 aa. UniProt lists P56817-1 as the canonical BACE1 product and maps it to ENST00000313005 ([UniProt BACE1](https://www.uniprot.org/uniprotkb/P56817-1/entry)).

The alternate used in the final analysis is BACE1 isoform 5, UniProt P56817-5, Ensembl ENST00000392937, project label **BACE1-202**, 401 aa. UniProt maps ENST00000392937 to P56817-5 in the same entry.

Sequence mapping shows an alternative N-terminal segment replacing canonical residues 1–20 and deletion of canonical residues 21–120. The remaining canonical sequence from approximately 121–501 is retained and renumbered in the alternate. This is best described as an **N-terminal alternative plus internal deletion**, with a net length reduction of 100 aa; it is not a simple terminal truncation.

The deletion removes canonical catalytic Asp93 and 12 of 28 experimentally defined verubecestat-contact residues. Catalytic Asp289 remains. This is why the candidate qualifies as mechanism B.

### Project DTU evidence

The transcript was measured in oligodendrocytes at:

- control usage: **8.5371%**;
- AD usage: **18.5023%**;
- absolute change: **+9.9653 percentage points**;
- transcript-level p-value: **0.7508**;
- gene-level adjusted p-value: **0.6554**;
- empirical FDR: **0.9771**.

This candidate must be described as **AD-observed or nominally increased**, not statistically significant. Its structural evidence is strong, but its disease-transcript evidence is weaker than the A and C candidates.

### Drug and development history

Verubecestat, also called **MK-8931** and **SCH 900931**, is an oral BACE1 inhibitor developed by **Merck & Co.** (Merck Sharp & Dohme; MSD outside the United States and Canada). Merck initiated the pivotal EPOCH program in late 2012 ([Merck announcement](https://www.merck.com/news/merck-initiates-phase-ii-iii-study-of-investigational-bace-inhibitor-mk-8931-for-treatment-of-alzheimers-disease/)).

The maximum Alzheimer disease phase reached was **Phase III**. The two key trials were:

1. **EPOCH, NCT01739348:** Phase II/III-to-Phase III study in mild-to-moderate AD, actual start **November 30, 2012**, Merck Sharp & Dohme sponsor, planned enrollment 1,958 in the pivotal Part I and 2,211 randomized across the recorded study. The trial was terminated early for futility ([ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT01739348); [NEJM report](https://www.nejm.org/doi/full/10.1056/NEJMoa1706441)).
2. **APECS, NCT01953601:** Phase III study in prodromal AD/amnesic mild cognitive impairment due to AD, actual start **November 5, 2013**, 1,454 participants, Merck sponsor. It was stopped in February 2018 after the monitoring committee concluded that superiority over placebo was unlikely ([ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT01953601); [Merck discontinuation announcement](https://www.merck.com/news/merck-announces-discontinuation-of-apecs-study-evaluating-verubecestat-mk-8931-for-the-treatment-of-people-with-prodromal-alzheimers-disease/)).

### Why it failed

EPOCH failed because verubecestat did not slow cognitive or functional decline. At week 78, neither the 12-mg nor 40-mg group significantly improved ADAS-Cog or ADCS-ADL relative to placebo; adverse events were more common in active-treatment groups. APECS likewise failed to prevent clinical progression, with higher-dose results favoring placebo in the published report ([APECS NEJM report](https://www.nejm.org/doi/full/10.1056/NEJMoa1812840)).

The project’s B hypothesis is not that this explains the entire Merck trial failure. Rather, it demonstrates a plausible isoform-specific failure mode: an AD-observed BACE1 product deletes part of the verubecestat pocket and catalytic Asp93, and docking shows that the alternate does not recover the canonical binding geometry.

## C — CACNA1D-214 / isradipine

### Protein identifiers and structural change

The canonical CaV1.3 alpha-1 subunit is encoded by CACNA1D. The project’s 2,161-aa canonical sequence corresponds to UniProt Q01668 and Ensembl ENST00000350061; Ensembl reports ENST00000350061 as a 2,161-aa CACNA1D transcript ([Ensembl transcript record](https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core%3Bg%3DENSG00000157388%3Br%3D3%3A53528683-53847760%3Bt%3DENST00000288139)).

The project alternative is labeled **CACNA1D-214** and is represented by Ensembl **ENST00000636627**. No curated UniProt isoform accession was assigned to this project-specific product. Its modeled protein length is 1,625 aa.

The alternative contains a 20-aa insertion near canonical residue 492 and lacks canonical residues 1606–2161, a **536-aa distal C-terminal truncation**. The mapped isradipine-contact region is retained and sequence-identical. This is mechanism C rather than B.

### Project DTU evidence

The transcript was measured in inhibitory neurons at:

- control usage: **1.8051%**;
- AD usage: **28.1081%**;
- absolute change: **+26.3031 percentage points**;
- raw DTU p-value: **7.7371×10⁻³⁸**;
- adjusted DTU p-value: **4.1153×10⁻³⁶**;
- permutation p-value: **0.0259**.

This is a statistically strong AD-associated isoform shift.

### Drug and development history

Isradipine is an L-type calcium-channel dihydropyridine blocker. It was originally developed by **Sandoz Pharmaceuticals** and approved for hypertension; Sandoz is now part of the Novartis group. The Parkinson study itself was not a conventional industry-sponsored registration program: it was investigator-led, sponsored by the **University of Rochester**, and funded/collaborated through NINDS, the Michael J. Fox Foundation, and the Parkinson Study Group ([FDA historical review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/pre96/19-546_Isradipine_statr.pdf); [ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT02168842)).

The maximum Parkinson disease phase reached was **Phase III**:

1. **STEADY-PD II:** an earlier Phase II tolerability/feasibility program that informed dose selection.
2. **STEADY-PD III, NCT02168842:** Phase III, randomized, double-blind, placebo-controlled study in early Parkinson disease. Enrollment began in **November 2014**, and 336 participants were enrolled ([trial-design paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5454402/)). Participants received immediate-release isradipine up to 5 mg twice daily versus placebo for 36 months.

### Why it failed

STEADY-PD III failed because isradipine did not slow Parkinson disease progression. The adjusted treatment effect on the primary UPDRS outcome was −0.27 points with p = 0.85; the study showed no meaningful difference from placebo ([randomized trial report](https://pmc.ncbi.nlm.nih.gov/articles/PMC7465126/)). The study had approximately 95% retention, so the negative result was not simply due to widespread dropout.

The project’s C hypothesis is different from the BACE1 mechanism. The direct pocket is retained, so a static docking score is not expected to reveal the key effect. Published CaV1.3 studies support the possibility that C-terminal splice variation changes gating and dihydropyridine sensitivity ([Huang et al.](https://pubmed.ncbi.nlm.nih.gov/23924992/); [Bock et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3234967/); [Ortner et al.](https://pubmed.ncbi.nlm.nih.gov/28592699/)). The proposed mechanism is altered channel-state occupancy, gating, trafficking, or coupling rather than simple loss of ligand contacts.

## Recommended wording for the presentation

| Example | Safe wording |
|---|---|
| A / FYN | “An AD-enriched novel FYN transcript encodes a predicted 115-aa product lacking the complete canonical saracatinib-binding kinase domain.” |
| B / BACE1 | “The AD-observed BACE1-202 product deletes catalytic Asp93 and 12/28 verubecestat-contact residues; canonical docking succeeds, but the alternate fails pose recovery. Its DTU increase is not statistically significant in the current dataset.” |
| C / CACNA1D | “The statistically supported AD-associated CACNA1D-214 isoform retains the isradipine pocket but deletes a distal C-terminal regulatory region, providing a hypothesis for altered state-dependent drug action.” |

Do not state that any isoform has been proven to cause a clinical-trial failure. The evidence supports candidate mechanisms that require independent validation of translation, localization, target engagement, and functional drug response.

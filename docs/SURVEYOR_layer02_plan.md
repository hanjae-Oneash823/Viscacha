# Surveyor — Layer 02 Conceptual Plan
## Candidate Specification and Evidence Assembly

**Pipeline context:** VISCACHA (Layer 0/1) → **Surveyor (Layer 2)** → AlphaFold (Layer 3) → Discovery Studio / CDOCKER (Layer 4/5)

**Version:** v01-draft  
**Status:** Conceptual — no code  
**Last updated:** 2026-06-13

---

## 1. Purpose and Design Philosophy

Surveyor is the automated candidate specification layer. It sits between the statistical output of VISCACHA and the structural biology steps that follow. Its job is to take each candidate isoform pair — defined entirely by a gene name, two ENST IDs, and a set of DTU statistics — and construct a complete, structured evidence dossier that answers three foundational questions before any compute is committed to structure prediction or docking:

1. **What is this protein, and where does the splice event fall relative to its functional architecture?**
2. **Is this protein already a drug target, and is any existing drug relevant to Alzheimer's disease?**
3. **What is the biological context — pathway membership, interaction partners, disease association — that justifies or challenges this candidate's inclusion in the docking study?**

The layer is fully automatic by design. When VISCACHA or the prioritization criteria change and new candidates are produced, Surveyor re-runs without manual intervention and produces updated dossiers. This automation is achieved by strict adherence to the Layer 1 output contract: Surveyor reads only from `layer1_candidates.yaml` and `layer1_expression.h5` and never requires human input mid-run.

Every algorithmic decision in Surveyor is logged with its rationale. The final output of each module is machine-readable JSON (consumed by downstream layers) and human-readable HTML (consumed by the researcher). No information is silently discarded — if a database query returns ambiguous or conflicting results, both values are preserved and flagged rather than arbitrarily resolved.

---

## 2. Input Contract

Surveyor depends on two files produced by Layer 1. These files must exist and pass a format validation check before any module runs.

### 2.1 `layer1_candidates.yaml`

The primary candidate manifest. Contains, for each candidate:

- Gene symbol and UniProt accession
- Tier classification (1 or 2)
- All significant transcripts with their ENST IDs, roles (disease_dominant / control_dominant / canonical_reference), PSI values per condition, ΔPSI, and adjusted p-values
- Cell type of significance and robust_to_braak flag
- Gene-level DE summary (log2FC, padj at gene level — separate from DTU)

### 2.2 `layer1_expression.h5`

The expression archive. Contains, for each candidate gene and each of its significant transcripts:

- Raw pseudo-bulk counts: shape (n_donors × n_cell_types)
- CP10K normalized counts: shape (n_donors × n_cell_types)
- Per-donor PSI values: shape (n_donors × n_cell_types), computed from raw counts, NaN where gene count is zero

Both files are versioned by the `viscacha_version` field in the YAML metadata. If the version field does not match the Surveyor-expected version, the run halts with an explicit error rather than silently processing mismatched data.

---

## 3. Output Contract

Surveyor produces three outputs per candidate and two global outputs per run.

### 3.1 Per-candidate: structured JSON dossier

`output/candidates/{gene}_dossier.json`

Machine-readable. Contains the full structured output of all eight modules. This is the file consumed by Layer 3 (AlphaFold submission preparation) and Layer 5 (docking setup). It must be complete and self-contained — downstream layers do not re-query databases.

Key top-level fields:

```
gene
enst_disease / enst_canonical
tier
splice_classification (Type A / B / C / D)
domain_overlap (list of affected domains with remapped residue coordinates)
docking_readiness (go / conditional / no-go with reason)
drug_target_status (known / novel)
ad_relevance (direct / indirect / none)
opentargets_score
brain_specificity_index
pathway_ad_flag
interaction_ad_genes (list)
```

### 3.2 Per-candidate: HTML report

`output/reports/{gene}_surveyor_report.html`

Human-readable. Self-contained single-file HTML with embedded CSS and inline SVG visualizations. Does not require a web server or external assets. Designed to be shared as a file attachment or archived alongside the manuscript.

Contains all eight module outputs organized into collapsible sections, three embedded visualizations (PSI strip plot, expression heatmap, isoform proportion bar), and the docking readiness checklist.

### 3.3 Per-candidate: AlphaFold submission package

`output/alphafold/{gene}/`

A directory containing the FASTA files and a submission metadata JSON ready for ColabFold input. Generated only for candidates with docking_readiness = go or conditional. Contains:

- `{gene}_{isoform}_disease.fasta`
- `{gene}_{isoform}_canonical.fasta` (if no PDB structure exists)
- `submission_params.json` (recommended ColabFold parameters, informed by sequence length and domain confidence)

### 3.4 Global: run summary

`output/surveyor_run_summary.md`

A single markdown file summarizing all candidates processed in the run: tier, splice classification, docking readiness verdict, and the one-line rationale for each verdict. Designed to be read in five minutes to understand the full output of the run.

### 3.5 Global: audit log

`output/surveyor_audit.log`

Timestamped log of every API call made, cache hits vs misses, any warnings raised, and the version of each database queried. Essential for reproducibility — the audit log allows any result to be traced back to a specific database state.

---

## 4. Module Architecture

Surveyor is organized as eight independent modules (M01–M08) plus an orchestrator. The modules are independent in the sense that each can succeed or fail without blocking the others, except for M03 which depends on M01 and M02. The orchestrator handles execution order and assembles the final output.

### Execution dependency graph

```
M01 (Sequence retrieval)  ──┐
                             ├──► M03 (Splice classification) ──► M08 (Report assembly)
M02 (Domain annotation)   ──┘                                          ▲
                                                                        │
M04 (Drug target)  ─────────────────────────────────────────────────►  │
M05 (Disease assoc.)  ──────────────────────────────────────────────►  │
M06 (Pathway)  ─────────────────────────────────────────────────────►  │
M07 (Expression)  ──────────────────────────────────────────────────►  │
M_VIS (Visualizations)  ────────────────────────────────────────────►  │
```

M01 and M02 run sequentially first. M03 runs after both complete. M04, M05, M06, M07, and M_VIS run in parallel using thread-based concurrency (I/O-bound API calls). M08 waits for all modules to complete before assembling the final output.

---

## 5. Module Specifications

---

### M01 — Sequence Retrieval

**Database:** Ensembl REST API  
**Fallback:** BioMart for bulk queries if REST rate limit is hit

#### What it retrieves

For every ENST ID in the candidate list (both isoforms of each candidate pair, plus the MANE Select canonical if not already included):

- Full protein sequence in FASTA format
- Exon structure: for each exon, genomic start, genomic end, CDS start phase, CDS end phase, exon rank within transcript
- Transcript biotype (protein_coding / retained_intron / nonsense_mediated_decay / etc.)
- MANE Select transcript ID for the gene (used when no canonical comparator is specified in the YAML)
- Protein length in amino acids
- CDS length in nucleotides

#### Key algorithmic decisions

**Exon comparison method — coordinate-based symmetric difference, not sequence alignment**

To identify which exons differ between the disease and canonical isoform, M01 performs a symmetric difference of exon sets by genomic coordinate rather than sequence alignment. The reasoning:

Sequence alignment can misplace exon boundaries when skipped exons share partial sequence similarity with adjacent exons — a known failure mode in alternatively spliced genes where exon sequences evolve under similar selective pressures. Coordinate-based comparison is unambiguous: an exon either exists at a genomic position in both transcripts or it does not. This is the correct approach for long-read data where transcript structure is directly resolved rather than inferred from short reads.

Matching is done by genomic overlap (not exact coordinate match) to handle cases where transcript annotation versions differ slightly between releases. An exon in transcript A is considered to match an exon in transcript B if their genomic intervals overlap by at least 80% of the shorter exon's length.

**CDS phase at exon boundaries**

For each exon unique to one isoform, M01 records the CDS phase at the splice boundary. Phase 0 means the exon boundary falls between codons — the splice event is in-frame and results in insertion or deletion of complete amino acids. Phase 1 or 2 means the boundary falls within a codon — the splice event is potentially frame-shifting downstream.

Frame-shifting events are not necessarily non-functional (they may resolve via compensatory exon inclusion) but they must be flagged explicitly, as a frame-shifted isoform may not produce a translatable protein. If the disease isoform has a frame-shifting exon boundary and is nonetheless highly expressed (high PSI_AD), this suggests either NMD escape or that the annotation is incomplete.

**Biotype check**

If the disease isoform has biotype `retained_intron` or `nonsense_mediated_decay`, M01 flags the isoform as potentially non-protein-coding. Surveyor does not automatically discard such isoforms — long-read sequencing frequently detects isoforms that short-read-based annotation has miscategorized. The flag is passed to M08 for inclusion in the readiness checklist.

#### Output

Per-transcript:
- FASTA protein sequence file
- Exon structure table (TSV)
- Exon symmetric difference table relative to comparator transcript
- Biotype and frame-shift flags

---

### M02 — Domain Annotation

**Database:** UniProt REST API (primary), InterPro API (supplementary for novel domain types)

#### What it retrieves

For the canonical UniProt accession of each candidate gene:

- All annotated sequence features: DOMAIN, BINDING, ACT_SITE, MOTIF, MOD_RES (post-translational modification sites), REGION, DISULFID, SIGNAL, TRANSIT
- For each feature: type, description, start residue, end residue, evidence type (experimental / by similarity / predicted)
- Active site residues with catalytic roles (where annotated)
- Known disease variants overlapping functional features (from UniProt variant annotations)
- Reviewed status of the UniProt entry (Swiss-Prot reviewed vs TrEMBL unreviewed)

#### Key algorithmic decisions

**Residue coordinate remapping to isoform sequences**

UniProt domain annotations are defined on the canonical UniProt sequence, which corresponds to the MANE Select or longest annotated isoform. Applying these residue numbers directly to a non-canonical isoform is incorrect when the alternatively spliced region falls before the domain in the linear sequence — any insertion or deletion shifts all downstream residue numbers.

M02 remaps domain boundaries from canonical to isoform coordinates using pairwise global alignment (BLOSUM62 substitution matrix, gap opening penalty -10, gap extension penalty -0.5). The alignment produces a position mapping between the two sequences, which is used to translate each UniProt feature's start/end residues into the isoform coordinate system.

When a domain feature maps to a region that is deleted in the disease isoform (i.e., the alternatively spliced exon encodes part of the domain), M02 records the domain as partially or fully absent in the disease isoform. This is a direct structural consequence that must be reported prominently in the M08 output.

**Evidence type weighting**

UniProt features annotated as "experimental evidence" are more reliable than those annotated "by similarity" or "predicted." M02 records the evidence type for each feature and weights it accordingly in the M03 splice classification:

- Experimental binding site overlapping the splice event → confirmed Type A
- "By similarity" binding site overlapping the splice event → probable Type A, flagged for manual review
- Predicted binding site → Type A candidate, lower confidence

**PTM site proximity**

Post-translational modification sites (phosphorylation, ubiquitination, acetylation) within 5 residues of the splice event boundary are flagged separately. A splice event that removes a known phosphorylation site changes the protein's regulatory properties independently of its binding pocket geometry. For kinase candidates (CAMK2B), this is especially relevant.

#### Output

Per-gene:
- Full feature table with original UniProt coordinates
- Remapped feature table in disease isoform coordinates
- Domain absence/truncation flags
- PTM proximity flags

---

### M03 — Splice Event Classification

**Dependencies:** M01 output, M02 output  
**External databases:** None — internal computation only

#### Purpose

M03 is the convergence point of sequence and domain information. It determines where the alternatively spliced exon(s) fall relative to the protein's functional architecture and assigns a classification that directly informs docking priority.

#### Classification scheme

**Type A — Splice event directly overlaps binding pocket**

The alternatively spliced exon encodes residues that are annotated as BINDING or ACT_SITE in UniProt (or InterPro where UniProt annotation is absent). The disease isoform has a structurally different drug binding pocket. This is the strongest possible signal for the hypothesis — the docking comparison is directly testing whether the pocket geometry change affects drug affinity.

Priority for docking: highest. Proceed regardless of other factors.

**Type B — Splice event adjacent to binding pocket**

The alternatively spliced exon does not directly encode binding residues, but falls within 15 residues (in the linear sequence) of a BINDING or ACT_SITE feature, OR the exon encodes part of a DOMAIN that contains the binding site (even if the binding residues themselves are retained).

The 15-residue threshold corresponds approximately to one alpha-helix turn. Residue substitutions within this range can propagate conformational changes to the binding site through secondary structure rearrangement. This is a well-established principle in allosteric pharmacology and is supported by mutational sensitivity analyses of binding pocket accessibility.

Priority for docking: high. Proceed with note that binding affinity change may be indirect.

**Type C — Splice event within a domain, distal from binding pocket**

The alternatively spliced exon encodes part of a DOMAIN feature, but the nearest BINDING or ACT_SITE is more than 15 residues away. The splice event changes the domain's structure but the binding pocket residues are fully conserved in both isoforms.

This classification does not preclude a docking result — allosteric effects can propagate across entire domains — but the expected ΔΔG between isoforms is smaller and less predictable. The scientific claim must be framed as allosteric modulation rather than direct pocket alteration.

Priority for docking: medium. Proceed if Tier 1, conditional if Tier 2.

**Type D — Splice event in unannotated or disordered region**

The alternatively spliced exon falls in a region with no UniProt functional annotation. This includes linker regions, intrinsically disordered regions, and N/C-terminal extensions. The binding pocket is fully conserved between isoforms.

A Type D classification does not mean the isoform switch is biologically irrelevant — the splice event may affect protein-protein interaction surfaces, localization signals, or degradation signals. But it means the canonical docking comparison (same ligand, same pocket, different isoform) is unlikely to produce a meaningful ΔΔG difference.

Priority for docking: low. Flag for PPI surface analysis instead. May still proceed if M04 identifies a drug that binds the alternatively spliced region specifically.

#### Distance calculation

The distance between the splice event and the nearest binding feature is calculated as the minimum residue distance in the disease isoform coordinate system (after M02 remapping), using the closest residue of the alternatively spliced exon to the closest annotated binding residue.

When the splice event removes a domain entirely (domain absent in disease isoform), the distance is recorded as 0 and the classification is automatically Type A, regardless of whether the binding site was within the removed region.

#### Output

Per-candidate:
- Type A / B / C / D classification
- Specific exon(s) responsible, with genomic and protein coordinates
- Distance to nearest binding feature (residues)
- Confidence level (based on UniProt evidence type)
- Narrative justification (one sentence, included verbatim in M08 report)

---

### M04 — Drug Target and Known Compounds

**Databases:** DrugBank (primary), ChEMBL (bioactivity data)

#### What it retrieves from DrugBank

- All drugs targeting the gene, organized by group: approved, investigational, experimental, withdrawn
- For each drug: name, DrugBank ID, mechanism of action, pharmacological action (inhibitor / activator / modulator), indication(s), approval status by region
- Known binding site on the protein target (where annotated in DrugBank)
- Drug-drug interaction liabilities for each compound (relevant for repurposing feasibility)

#### What it retrieves from ChEMBL

- All bioactivity records for the target: Kd, IC50, Ki, EC50 values with assay conditions
- Assay type: binding assay vs functional assay (binding assays are more directly relevant to docking validation)
- Compound structures (SMILES) for compounds with Kd or IC50 < 1 µM — these are the candidates for the docking study
- Clinical phase of each compound (Phase I / II / III / approved)

#### Key algorithmic decisions

**AD relevance classification — two-tier approach**

AD relevance is determined by two separate criteria applied independently:

*Tier 1 — Direct indication match:* The drug's indication in DrugBank contains "Alzheimer," "dementia," or a related MeSH term. These are drugs already in clinical use or trials for AD. They are the primary candidates for the canonical vs. disease isoform comparison.

*Tier 2 — Indirect mechanistic relevance:* The drug targets a pathway (cross-referenced with M06 Reactome output) directly implicated in AD pathology. The AD-relevant pathway list used for this cross-reference:

- Amyloid precursor protein processing and amyloid-β production
- Tau protein phosphorylation and aggregation
- Mitochondrial quality control and mitophagy
- Endosomal-lysosomal trafficking
- Neuroinflammation (microglial activation, cytokine signaling)
- Axonal transport (kinesin/dynein-mediated)
- Synaptic plasticity and LTP

A drug classified as Tier 2 is a repurposing candidate. It was not designed for AD but its mechanism of action engages an AD-relevant process through the target gene. Tier 2 candidates are included in the repurposing screen (Layer 5 / Layer 6) and labeled as such in the report.

**Rationale for two-tier approach:** Many mechanistically relevant compounds are not labeled for AD in DrugBank because they were developed for other indications (cancer, pain, psychiatric conditions). A kinase inhibitor developed for glioblastoma that targets CAMK2B is not labeled as an AD drug, but if CAMK2B shows isoform switching in AD neurons, that inhibitor is a legitimate repurposing candidate. The two-tier approach captures these without conflating them with directly approved AD drugs.

**ChEMBL potency threshold**

The module retrieves all bioactivity records but only flags compounds with measured affinity (Kd or Ki) below 1 µM as docking candidates. The 1 µM threshold is standard for defining a compound as a specific binder rather than a promiscuous or non-specific interactor. Compounds with IC50 values only (not Kd/Ki) are included but de-prioritized, as IC50 reflects functional inhibition which may not correspond to direct binding at the pocket being docked.

**Binding site annotation cross-reference**

Where DrugBank provides a known binding site for a drug-target interaction, M04 cross-references this with the M03 splice classification. If the drug binds at a site that is structurally altered by the splice event (Type A or B classification), this is flagged as a high-priority docking pair: the drug was designed for the canonical pocket, and the disease isoform has an altered version of that same pocket.

#### Output

Per-gene:
- Drug table: name, group, indication, AD relevance tier, mechanism, phase
- Bioactivity table: compound, assay type, Kd/IC50/Ki, units, ChEMBL assay ID
- Docking candidate list: SMILES + compound ID for all sub-µM binders
- AD relevance summary: n_direct, n_indirect, n_novel (no existing drugs)

---

### M05 — Disease Association

**Database:** OpenTargets Platform API (GraphQL endpoint)

#### What it retrieves

For each candidate gene against Alzheimer's disease (EFO:0000249):

- Overall association score (0–1, weighted across evidence types)
- Evidence breakdown by data type:
  - Genetic association (GWAS, rare variant studies)
  - Somatic mutation (COSMIC, cancer genomics — low relevance for AD but included for completeness)
  - Differential expression (GTEx, Expression Atlas)
  - Animal model evidence (MGI, RGD)
  - Known drug evidence (clinical trials in AD)
  - Literature evidence (text-mining score)
- For genetic evidence specifically: lead SNPs, effect allele frequency, GWAS p-value, population
- For context: association scores against Parkinson's disease (EFO:0002508), frontotemporal dementia (EFO:0001627), and ALS (EFO:0000253) — to assess AD specificity vs general neurodegeneration

#### Key algorithmic decisions

**Score interpretation — no hard threshold applied**

OpenTargets scores are not applied as binary pass/fail thresholds. The reasoning is specific to this project: several candidates (MGRN1, RHOT1) are biologically plausible AD targets whose roles in lysosomal trafficking and mitochondrial dynamics are mechanistically linked to AD pathology, but they are understudied relative to canonical AD targets like APP, PSEN1, or MAPT. An understudied gene will have a low OpenTargets score not because evidence contradicts its AD relevance, but because the research community has not yet systematically investigated it.

Applying a minimum score threshold would systematically exclude potentially novel findings — exactly the targets that a long-read isoform study is best positioned to discover.

Instead, the score is reported with its evidence type breakdown, and a flag is applied:

- Score ≥ 0.5: established evidence → label "supported target"
- Score 0.1–0.49: partial evidence → label "emerging target"
- Score < 0.1 with no genetic evidence: label "novel target — low prior evidence" (scientifically interesting, not disqualifying)
- Score < 0.1 with contradicting evidence (e.g. protective variants): label "conflicting evidence — review manually"

**AD specificity ratio**

For each candidate, M05 computes:

```
AD_specificity = opentargets_score_AD / mean(opentargets_score_other_NDs)
```

Where other_NDs includes Parkinson's, FTD, and ALS. A ratio > 2 indicates the gene is more strongly associated with AD specifically than with neurodegeneration in general. A ratio near 1 indicates a general neurodegeneration target. Both are scientifically interesting but frame the candidate's narrative differently in the report and in any eventual manuscript.

**Genetic evidence extraction**

Where GWAS evidence exists, M05 extracts the specific SNP(s) and their genomic coordinates, and checks whether any lead SNPs fall within the gene's exon boundaries or splice sites. A GWAS hit within a splice site is direct genetic evidence that splicing at this locus is disease-relevant — this would substantially strengthen the case for the candidate and is flagged explicitly in the M08 report.

#### Output

Per-gene:
- OpenTargets association score and evidence type breakdown
- AD target label (supported / emerging / novel / conflicting)
- AD specificity ratio vs other neurodegenerative diseases
- GWAS SNP table (where available) with splice site proximity flag
- Clinical trial evidence in AD (cross-reference with M04)

---

### M06 — Pathway and Interaction Analysis

**Databases:** STRING API (protein-protein interactions), Reactome Content Service (pathway membership)

#### What it retrieves from STRING

- Top interaction partners ranked by STRING combined score, filtered to combined score ≥ 0.7 (high confidence)
- For each interaction: evidence channels present (experimental, co-expression, database, co-occurrence, text-mining, homology)
- Functional enrichment of the interaction network: GO Biological Process terms (FDR < 0.05), KEGG pathways, Reactome pathways, Disease terms

#### What it retrieves from Reactome

- All pathways containing the gene, organized by top-level pathway hierarchy
- For each pathway: pathway ID, name, species, evidence type
- AD-relevant pathway membership specifically (see pathway list in M04 section)
- Whether the gene is in the same Reactome pathway as any known AD drug target (APP, BACE1, PSEN1, MAPT, APOE, BIN1, CLU, CR1, PICALM)

#### Key algorithmic decisions

**STRING confidence threshold — 0.7 (high confidence), not 0.4 (medium)**

The STRING default threshold of 0.4 includes many text-mining-only interactions that are not experimentally validated. For a drug target specification, the interaction partners you report must be experimentally supported, because they define the functional context you use to justify why a drug targeting this protein would affect AD pathology.

At 0.7, STRING interactions are predominantly backed by experimental evidence (co-immunoprecipitation, yeast two-hybrid, proximity ligation assay) or multiple independent co-expression datasets. False positive interactions at this threshold are substantially lower than at 0.4.

The tradeoff is that rare or understudied interactions may be excluded. To mitigate this, M06 also runs a lower-confidence query (threshold 0.4) and reports the additional interactions separately, labeled "moderate confidence — not experimentally validated" so they are available but not conflated with high-confidence interactions.

**AD-gene interaction flag**

M06 cross-references the top 20 high-confidence interaction partners against a curated list of established AD-related proteins. The curated list includes:

- Core genetic risk genes: APP, PSEN1, PSEN2, APOE, TREM2, BIN1, CLU, CR1, PICALM, ABCA7
- Core pathological proteins: MAPT (tau), SNCA (alpha-synuclein)
- Major AD drug targets: BACE1, ADAM10, GSK3B, CDK5
- Established AD pathway proteins: ADAM17, NCSTN, APH1A, PSENEN (gamma-secretase complex)

If any of the candidate's top interaction partners appear in this list, M06 flags the interaction with the AD gene name and the STRING evidence type. A candidate that directly interacts with MAPT (tau) or APP — even if it has a low OpenTargets score itself — is immediately placed in a different scientific category: it is mechanistically positioned to affect AD pathology through its interaction network regardless of its own genetic association.

**Reactome pathway hierarchy traversal**

Reactome pathways are hierarchical. A gene may appear in a very specific sub-pathway (e.g., "Kinesin-mediated transport of APP-containing vesicles") that belongs to a broader pathway (e.g., "Axonal transport") that belongs to a top-level pathway (e.g., "Neuronal System"). M06 records the full hierarchy for each pathway membership so the report can present both the specific mechanistic context and the broader functional category.

#### Output

Per-gene:
- High-confidence interaction partner table (top 20, score ≥ 0.7) with evidence channels
- Moderate-confidence interactions (score 0.4–0.69) labeled separately
- AD-gene interaction flags with interaction type
- Reactome pathway membership list with full hierarchy
- AD-relevant pathway flags
- GO Biological Process enrichment (top 10 terms by FDR) for the interaction network

---

### M07 — Expression Context

**Primary data source:** `layer1_expression.h5` (from Layer 1 — no external API for primary data)  
**Supplementary database:** GTEx Portal API (bulk tissue expression for off-target risk assessment)

#### What it computes from Layer 1 HDF5

**Per-donor PSI per cell type**

For the disease isoform and the canonical comparator, M07 extracts the full (n_donors × n_cell_types) PSI matrix from the HDF5. From this it computes:

- Per-condition mean PSI ± standard deviation for each cell type
- Per-donor PSI values for the significant cell type, stratified by condition (AD / Control / Active Control)
- Cell-type specificity of the isoform switch: is the ΔPSI pattern consistent across all cell types, or confined to the statistically significant cell type?
- Donor-level variance within each condition: high within-condition variance suggests the isoform switch may be driven by a subset of donors (e.g., those with more advanced pathology) rather than being a uniform disease effect

**Gene-level expression context**

From the gene-level expression data in the HDF5:

- Per-condition mean expression (CP10K) across all cell types
- Log2FC at the gene level (already computed in Layer 1 and stored in YAML, retrieved here for report integration)
- Classification of the candidate as one of three expression patterns:
  - *Pure isoform switch:* No significant gene-level DE (log2FC < 0.5, padj > 0.05) but significant DTU — the gene's total output is constant, but the isoform proportion shifts. This is the cleanest signal for the hypothesis.
  - *Combined DE + DTU:* Both gene-level expression and isoform proportion change in AD. The isoform switch is real but occurs in the context of broader transcriptional dysregulation at this locus.
  - *Complex:* Gene expression changes in one direction while the disease isoform proportion changes in the other — the disease isoform is becoming proportionally dominant even as total gene expression falls, or vice versa. Requires careful narrative framing.

**Cell-type expression breadth**

M07 identifies which cell types express the candidate gene above a minimum threshold (mean CP10K ≥ 1 across donors in any condition). This defines the cell-type expression breadth of the gene — whether it is narrowly expressed in the statistically significant cell type or broadly expressed across all brain cell types. A gene expressed only in excitatory neurons is a more specific therapeutic target than one expressed ubiquitously across all cell types.

#### What it retrieves from GTEx

- Median TPM expression across all GTEx v10 tissues (for off-target risk assessment)
- Brain-region-specific expression: frontal cortex BA9, hippocampus, anterior cingulate cortex BA24, cortex, putamen, caudate, cerebellum
- Cardiac tissue expression: left ventricle, atrial appendage (for cardiac off-target risk)
- Hepatic expression: liver (for metabolic off-target risk)

**Brain specificity index (BSI)**

```
BSI = mean(brain region TPM) / mean(all tissue TPM)
```

- BSI > 5: brain-specific → favorable target profile, low peripheral off-target risk
- BSI 2–5: brain-enriched → moderate off-target risk, manageable with appropriate formulation
- BSI 1–2: ubiquitous → high off-target risk, requires careful therapeutic window analysis
- BSI < 1: peripherally enriched → unfavorable target profile for CNS drug development

**Limitation explicitly noted in report:** GTEx is bulk tissue RNA-seq, not single-cell, and does not have isoform resolution. The BSI is computed on gene-level expression and does not distinguish between isoforms. A gene may be brain-enriched at the gene level but the disease-specific isoform may be expressed in other tissues — this cannot be resolved from GTEx alone. The Layer 1 single-cell data is the authoritative source for isoform-specific expression in brain cell types; GTEx provides tissue-level off-target risk context only.

#### Output

Per-gene:
- PSI matrix (disease isoform): donors × cell types
- PSI matrix (canonical isoform): donors × cell types
- Per-condition PSI summary table (mean ± SD per cell type)
- Expression pattern classification (pure switch / combined DE+DTU / complex)
- Cell-type expression breadth map
- GTEx brain specificity index
- Cardiac and hepatic expression flags (off-target risk)
- Three visualization datasets (consumed by M_VIS — see below)

---

### M_VIS — Visualization Generation

**Dependencies:** M07 output (expression data), M01 output (sequence annotation)  
**External libraries:** matplotlib, seaborn (for SVG generation to embed in HTML report)

M_VIS runs in parallel with M04–M07 and generates the three visualizations embedded in the HTML report. Visualizations are generated as SVG (scalable, embeddable, accessible) rather than PNG.

**Visualization 1 — Per-donor PSI strip plot**

For the significant cell type only. Each donor is a point. X-axis: condition (Control / Active Control / AD). Y-axis: PSI value for the disease isoform (0–1). Points are colored by condition. The mean ± SD per condition is overlaid as a horizontal line with error bars. Individual donor IDs are available on hover (implemented as SVG title elements).

This is the primary visualization for communicating the isoform switch. It shows both the group-level effect (mean shift) and the individual-level variance (point scatter). A reader should be able to assess from this plot alone whether the switch is driven by all AD donors or a subset.

**Visualization 2 — Gene expression heatmap**

X-axis: all cell types. Y-axis: all donors, grouped and sorted by condition. Color: CP10K normalized gene expression (log scale, capped at 99th percentile). Cell types with significant DTU hits are indicated by a border on the column.

This shows at a glance where the gene is expressed, whether expression varies by condition, and which cell types are the primary expressers.

**Visualization 3 — Isoform proportion stacked bar**

Three bars (Control / Active Control / AD), stacked to 100%. Each segment represents one isoform's mean PSI for the significant cell type. The disease isoform is colored distinctly from the canonical/control isoform. Minor isoforms (PSI < 5% in all conditions) are collapsed into a single "other" segment.

This replicates the VISCACHA output format for the specific candidate, making the Surveyor report self-contained without requiring reference back to VISCACHA output files.

---

### M08 — Report Assembly

**Dependencies:** All modules (M01–M07, M_VIS)  
**Libraries:** Jinja2 (HTML templating), PyYAML, h5py, json

M08 is the final step. It waits for all modules to complete, then assembles the per-candidate JSON dossier and HTML report.

#### JSON dossier assembly

The JSON dossier is a structured merge of all module outputs under a standardized schema. Every field has a defined key, type, and whether it is required or optional. Downstream layers (AlphaFold submission, docking setup) read specific fields from this JSON — if a field is missing because a module failed, the downstream layer raises an explicit error rather than silently using a default value.

#### HTML report structure

The HTML report is a single self-contained file organized into the following sections:

**Header:** Gene name, tier, cell type of significance, ΔPSI, padj, date of Surveyor run, database versions queried.

**Executive summary:** Three-sentence plain-language summary of what this candidate is, why it was selected, and what the docking hypothesis is for this specific isoform pair. Auto-generated from structured fields — not a language model output.

**Docking readiness checklist:** The most important section for workflow automation. A checklist of conditions required before proceeding to AlphaFold and CDOCKER, with automated pass/fail/warning for each:

```
✓ Protein sequence retrieved — {n} aa (disease), {n} aa (canonical)
✓ / ✗ Biotype: protein-coding confirmed / WARNING: potential NMD isoform
✓ / ✗ Exon diff computed — {n} exons unique to disease isoform
✓ Splice classification: Type {A/B/C/D} — {one-line rationale}
✓ / ✗ / ⚠ Canonical structure: PDB {ID} available / AlphaFold required / low-confidence region
✓ / ✗ Known drugs found: {n} direct AD / {n} indirect / none
✓ / ⚠ OpenTargets: {score} ({label})
✓ / ⚠ Brain specificity: BSI = {value} ({interpretation})
⚠ robust_to_braak = FALSE — Braak stage confounding, note in manuscript
→ VERDICT: PROCEED / CONDITIONAL / NO-GO
   Reason: {one line}
```

Verdict logic:
- **PROCEED:** Type A or B splice classification, protein-coding confirmed, at least one known drug or reasonable druggability evidence
- **CONDITIONAL:** Type C classification, OR uncertain biotype, OR no known drugs but strong disease association — proceed with explicit caveats noted
- **NO-GO:** Type D classification AND no known drugs AND BSI < 1 — flag for manual review, do not automatically proceed to AlphaFold

**Section 1 — Isoform architecture:** Exon structure diagram (schematic, not to genomic scale), protein domain map with the alternatively spliced region highlighted, Type A/B/C/D classification with distance to nearest binding feature.

**Section 2 — Drug target profile:** Drug table (DrugBank), bioactivity table (ChEMBL), AD relevance classification for each compound.

**Section 3 — Disease association:** OpenTargets score panel, evidence type breakdown, GWAS SNPs (if any), AD specificity ratio.

**Section 4 — Biological context:** STRING interaction network summary (top 10 partners, AD-gene flags), Reactome pathway table, GO enrichment summary.

**Section 5 — Expression profile:** Three embedded SVG visualizations (PSI strip plot, expression heatmap, isoform proportion bar), expression pattern classification, brain specificity index, off-target risk flags.

**Section 6 — AlphaFold submission parameters:** Recommended ColabFold parameters for the disease isoform, informed by protein length and domain confidence. PDB structure recommendation for canonical isoform (download link if available).

**Footer:** Audit trail — list of all API calls made, cache status, database versions, any warnings raised during the run.

---

## 6. Caching Strategy

All external API responses are cached locally to avoid redundant calls when the pipeline re-runs with the same or overlapping candidate lists.

**Cache location:** `cache/` directory in the Surveyor working directory.

**Cache key:** MD5 hash of the full request URL and query parameters. This ensures different queries to the same endpoint produce different cache keys.

**Cache TTL by database:**

| Database | TTL | Rationale |
|---|---|---|
| Ensembl REST | 90 days | Annotation changes slowly; new releases quarterly |
| UniProt | 90 days | Swiss-Prot entries are manually curated; infrequent changes |
| DrugBank | 14 days | Approval statuses and clinical trial phases can change monthly |
| ChEMBL | 30 days | New bioactivity data deposited regularly |
| OpenTargets | 30 days | Platform updated with new GWAS data regularly |
| STRING | 90 days | Interaction database rebuilt with major releases only |
| Reactome | 90 days | Pathway annotations change slowly |
| GTEx | 180 days | GTEx v10 is a versioned static dataset; changes only with new releases |

**Cache invalidation:** A `--force-refresh` flag in the orchestrator bypasses all cache lookups and re-queries all databases. This should be used when a new database release is announced that may affect the analysis.

**Cache validation:** On each cache read, the module checks whether the database version embedded in the cached response matches the current expected version. If a version mismatch is detected, the cached entry is discarded and a fresh query is made. The expected version for each database is defined in `config.yaml`.

---

## 7. Error Handling and Partial Failure Policy

Because Surveyor queries seven external databases, partial failure (one database unreachable, rate-limited, or returning malformed data) is expected in routine operation. The failure policy is:

**Non-blocking failures (pipeline continues):**
- M04 DrugBank or ChEMBL API failure → report drug section as "unavailable — API error" but continue; docking readiness checklist uses "drug data unavailable" rather than "no drugs found"
- M05 OpenTargets API failure → report association as "unavailable" with neutral label
- M06 STRING or Reactome failure → report pathway section as "unavailable"
- M07 GTEx API failure → report BSI as "unavailable"; off-target risk assessment skipped

**Blocking failures (pipeline halts for this candidate):**
- M01 Ensembl failure → protein sequence unavailable; cannot proceed to M03 or any downstream layer. Candidate is logged as "failed — sequence retrieval error" in the run summary.
- M02 UniProt failure → domain annotation unavailable; M03 cannot classify splice event. Candidate is logged as "failed — domain annotation error." M04–M07 still run and produce output for the partial dossier, but the docking readiness verdict is set to "CONDITIONAL — manual splice classification required."

**All errors are logged** to `surveyor_audit.log` with the full error message, HTTP status code (for API errors), and timestamp. The run summary reports the number of successful, partial, and failed candidates.

---

## 8. Reproducibility Requirements

Surveyor is part of a pipeline intended to support publication. The following requirements ensure that the analysis can be reproduced from the Layer 1 output files alone:

1. **Database version recording:** Every API response includes version metadata (Ensembl release, UniProt release date, ChEMBL version, etc.). These are recorded in the audit log and embedded in the JSON dossier under a `database_versions` field.

2. **Parameter recording:** All thresholds, confidence cutoffs, and classification rules used in the run are recorded in the JSON dossier under a `parameters` field. If a threshold is changed between runs, the parameter record reflects the value used for that specific run.

3. **Deterministic output:** Given the same input YAML, the same cached API responses, and the same parameter configuration, Surveyor produces identical output. There are no random or stochastic steps.

4. **Input file hashing:** The MD5 hashes of `layer1_candidates.yaml` and `layer1_expression.h5` are recorded in the audit log at the start of each run, so it is always possible to verify which Layer 1 output produced a given Surveyor run.

---

## 9. Configuration File Structure

`config.yaml` — the single configuration file for all Surveyor parameters. All thresholds, API endpoints, and behavioral flags live here. No hardcoded values in module code.

```yaml
# surveyor_config.yaml

version: "01"
viscacha_version_expected: "v01"

# API endpoints
apis:
  ensembl: "https://rest.ensembl.org"
  uniprot: "https://rest.uniprot.org/uniprotkb"
  drugbank: "https://api.drugbank.com/v1"         # requires API key
  chembl: "https://www.ebi.ac.uk/chembl/api/data"
  opentargets: "https://api.platform.opentargets.org/api/v4/graphql"
  string: "https://string-db.org/api"
  reactome: "https://reactome.org/ContentService"
  gtex: "https://gtexportal.org/api/v2"

# API authentication
auth:
  drugbank_api_key: "${DRUGBANK_API_KEY}"  # loaded from environment variable

# Thresholds
thresholds:
  string_confidence_high: 0.7
  string_confidence_moderate: 0.4
  chembl_potency_uM: 1.0                  # Kd/Ki threshold for docking candidates
  psi_ad_minimum: 0.20                    # minimum PSI_AD for disease isoform
  delta_psi_minimum: 0.10                 # minimum |ΔPSI| for inclusion
  splice_domain_proximity_residues: 15    # Type B/C boundary
  gtex_expression_minimum_tpm: 1.0       # minimum for "expressed" in GTEx
  brain_specificity_enriched: 2.0        # BSI threshold for brain-enriched
  brain_specificity_specific: 5.0        # BSI threshold for brain-specific
  cell_type_expression_minimum_cp10k: 1.0

# Classification labels
opentargets_labels:
  supported: 0.5
  emerging: 0.1

# AD pathway list (used in M04 and M06)
ad_pathways:
  - "amyloid precursor protein processing"
  - "tau protein phosphorylation"
  - "mitophagy"
  - "endosomal-lysosomal trafficking"
  - "neuroinflammation"
  - "axonal transport"
  - "synaptic plasticity"
  - "long-term potentiation"

# AD gene list (used in M06 interaction flagging)
ad_genes:
  - APP
  - PSEN1
  - PSEN2
  - APOE
  - TREM2
  - BIN1
  - CLU
  - CR1
  - PICALM
  - ABCA7
  - MAPT
  - SNCA
  - BACE1
  - ADAM10
  - GSK3B
  - CDK5

# Docking readiness verdict rules
verdict_rules:
  proceed:
    - splice_type: [A, B]
    - biotype: protein_coding
    - drugs_or_druggability: true
  conditional:
    - splice_type: [C]
    - biotype_uncertain: true
    - no_known_drugs: true
  no_go:
    - splice_type: [D]
    - no_known_drugs: true
    - brain_specificity_index_below: 1.0

# Caching
cache:
  directory: "cache"
  ttl_days:
    ensembl: 90
    uniprot: 90
    drugbank: 14
    chembl: 30
    opentargets: 30
    string: 90
    reactome: 90
    gtex: 180
  force_refresh: false

# Concurrency
execution:
  max_parallel_modules: 4     # M04-M07 run in parallel; limited by API rate limits
  request_timeout_seconds: 30
  max_retries: 3
  retry_backoff_seconds: 5

# Output
output:
  candidates_dir: "output/candidates"
  reports_dir: "output/reports"
  alphafold_dir: "output/alphafold"
  run_summary: "output/surveyor_run_summary.md"
  audit_log: "output/surveyor_audit.log"
```

---

## 10. Open Questions Requiring Resolution Before Implementation

The following questions must be answered before coding begins. Each unresolved question has a potential impact on implementation that is noted.

**Q1 — DrugBank access method**

Does the lab have a DrugBank academic API key, or will the analysis use the monthly XML bulk download? The API is more convenient for automated queries. The XML dump is more complete and not rate-limited. Impact: M04 implementation differs substantially between the two approaches.

**Q2 — Node212 external network access**

The server has confirmed HTTPS blocking for Anthropic domains. Is this blocking specific to Anthropic (i.e., all other external HTTPS is permitted), or is there a broader policy? Test:

```bash
curl -s "https://rest.ensembl.org/lookup/id/ENST00000445352?content-type=application/json" | head -c 200
```

If this fails, all seven API modules are blocked and the entire retrieval strategy shifts to pre-downloaded database dumps queried locally. Impact: complete architectural change to offline mode.

**Q3 — Gene-level DE results from VISCACHA**

Does the existing VISCACHA pipeline produce gene-level differential expression statistics (log2FC, padj at the gene level, separate from DTU)? The `gene_expression` block in the Layer 1 YAML includes these values. If they do not exist in current VISCACHA output, Layer 1 must be extended with a DESeq2 gene-level DE step before the export function runs. Impact: Layer 1 modification required.

**Q4 — Conda environment specification**

What Python version and conda environment name should Surveyor use on node212? This determines the `requirements.txt` and environment setup instructions. Impact: environment setup and documentation.

---

## 11. Relationship to Adjacent Layers

| Layer | Direction | Interface |
|---|---|---|
| VISCACHA (Layer 0/1) | Upstream → Surveyor | `layer1_candidates.yaml` + `layer1_expression.h5` |
| AlphaFold submission (Layer 3) | Surveyor → Downstream | `output/candidates/{gene}_dossier.json` + `output/alphafold/{gene}/` |
| Discovery Studio preparation (Layer 4) | Surveyor → Downstream | `output/candidates/{gene}_dossier.json` (binding site coordinates, domain map) |
| Repurposing screen (Layer 6) | Surveyor → Downstream | Docking candidate SMILES list from M04 + ChEMBL potency data |

Surveyor does not communicate with layers downstream of Layer 3 directly. Layers 4–6 read the dossier JSON produced here. Any change to the dossier schema must be versioned and communicated to downstream layer implementations.

---

*End of Surveyor Layer 02 Conceptual Plan v01-draft*

# Surveyor — Layer 02 Conceptual Plan
## Candidate Specification and Evidence Assembly

**Pipeline context:** VISCACHA (Layer 0/1) → **Surveyor (Layer 2)** → AlphaFold (Layer 3) → Discovery Studio / CDOCKER (Layer 4/5)

**Version:** v2.2-draft  
**Status:** In build — foundation + M01 implemented and validated against live data  
**Last updated:** 2026-06-13  
**Supersedes:** v2.1-draft (same file)

---

## Changelog

### v2.0 (from v01)
| Section | Change |
|---|---|
| §2 Input contract | Replaced YAML/HDF5 wrapper with direct reads from existing Layer 0/1 files |
| §2.5 | Defined the one required Layer 1 addition (step15 gene-level DE) |
| §3 | Tier removed from inputs; added to M03 output |
| §5 M01 | Added comparator-selection rule for single-hit genes (MANE Select fallback) |
| §5 M03 | Tier assigned here, not read from input |
| §5 M07 | Multi-cell-type handling; per-donor PSI from Layer 0 pseudobulk CSVs |
| §5 M_VIS | One visualization set per significant cell type |
| §5 M08 | Multi-cell-type layout; Braak robustness rule defined |
| §10 | Q2 (network) removed — confirmed open; Q3 resolved → step15 |

### v2.2 (empirical validation against the 12 hits, during build)
| Change | Section | Rationale |
|---|---|---|
| Added **Type N** (no protein-level change) | §5 M03, M08 verdict | M01 run on real data showed 5/12 hits encode a protein identical to the comparator (or are non-coding). They have no docking hypothesis; routed to a regulatory analysis path (full dossier, no AlphaFold). |
| Dropped exon-length **frame-shift proxy** | §5 M01 | False-positives on UTR exons (EIF4G2-208). Protein-sequence identity + length-delta is exact and authoritative. |
| **Broadened functional-feature set** | §5 M03 | Real UniProt entries annotate functional cores variously: KLC1 = Repeat/Region, MGRN1 = Zinc finger/Motif, not DOMAIN/BINDING. Keying only on DOMAIN/BINDING/ACT_SITE would misclassify them as Type D. |
| Comparator fallback when **sig tx == MANE Select** | §5 M01 | RHOT1-205 is itself the MANE Select; self-comparison is meaningless. Falls through to highest-expressed control transcript (RHOT1-201). |
| Output path convention → `outputs/layer2/` | §9 | Matches existing `outputs/layer0`, `outputs/layer1` repo layout. |

### v2.1 (gap patches)
| Gap | Section changed | Fix |
|---|---|---|
| No PDB lookup | §5 M02, §3.3 | PDB cross-references retrieved from UniProt in M02; used in checklist and AlphaFold package |
| Active Control PSI for comparator missing | §2.3, §5 M07 | Sensitivity CSVs added to input contract; M07 reads them for Active Control PSI of all transcripts |
| MANE Select comparator may not be in pseudobulk | §5 M07 | Explicit fallback defined: comparator-absent case handled gracefully |
| Same-sign ΔPSI edge case | §5 M01 | Rule added: both transcripts positive/negative ΔPSI → treat each as single-hit, MANE Select as comparator |
| 3+ significant transcripts per gene | §5 M01 | Rule added: MANE Select is structural reference for all; all compared against it |
| Gene-level vs. cell-type-level execution undefined | §4 (new §4.1) | Orchestrator fully specified: gene-level modules run once per gene, cell-type modules per (gene, cell_type) |
| Orchestrator unspecified | §4 (new §4.1) | Orchestrator loop, deduplication, and thread pool described |
| Input validation absent | §2.6 (new) | Pre-run validation step specified |
| Gene-level index page unspecified | §3.2 | Formally specified as output with path, content, and trigger condition |
| SMILES → 3D conversion unaddressed | §5 M04 | Layer 2 owns this; RDKit ETKDG added to M04 post-processing |
| Verdict logic priority conflicts | §5 M08 | Explicit priority order defined |
| r = 0.914 hardcoded | §5 M08 | Computed at runtime from metadata; stored in audit log |
| Sensitivity CSVs absent from §2 | §2.3 | Added as formal input |

---

## 1. Purpose and Design Philosophy

Surveyor is the automated candidate specification layer. It sits between the statistical output of VISCACHA and the structural biology steps that follow. Its job is to take each candidate isoform — defined by a gene name, one or two ENST IDs, and a set of DTU statistics — and construct a complete, structured evidence dossier that answers three foundational questions before any compute is committed to structure prediction or docking:

1. **What is this protein, and where does the splice event fall relative to its functional architecture?**
2. **Is this protein already a drug target, and is any existing drug relevant to Alzheimer's disease?**
3. **What is the biological context — pathway membership, interaction partners, disease association — that justifies or challenges this candidate's inclusion in the docking study?**

The layer is fully automatic. Surveyor reads directly from files that Layer 0 and Layer 1 already produce. When Layer 0 or Layer 1 is re-run with different data or thresholds, Surveyor re-runs on the new files without modification. No intermediate manifest needs to be regenerated or manually edited between pipeline runs.

Every algorithmic decision is logged with its rationale. No information is silently discarded — if a database query returns ambiguous or conflicting results, both values are preserved and flagged rather than arbitrarily resolved.

---

## 2. Input Contract

Surveyor reads from five sources. The first four exist after any complete Layer 0/1 run. The fifth requires a new step in Layer 1 (see §2.5). A validation step (§2.6) runs before any module starts.

### 2.1 `outputs/layer1/dtu_significant_all_celltypes.csv`

The primary candidate list. One row per significant transcript. Columns used by Surveyor:

- `transcript_id` — GENCODE transcript name (e.g. `CAMK2B-204`)
- `gene_id` — gene symbol (e.g. `CAMK2B`)
- `padj_gene`, `padj_tx` — stageR-adjusted p-values
- `psi_AD`, `psi_ctrl`, `delta_psi` — condition-level mean PSI
- `psi_active_ctrl` — Active Control descriptive PSI
- `robust_to_braak` — logical
- `cell_type` — cell type in which significance was detected

If a gene appears in more than one cell type, it produces multiple rows. Surveyor treats each unique `(gene_id, cell_type)` combination as one analysis unit for cell-type-specific modules, and each unique `gene_id` as one unit for gene-level modules (see §4.1).

### 2.2 `outputs/layer0/filtered_adata/adata_tx_step03.h5ad`

Used exclusively to map transcript names to Ensembl IDs. Surveyor reads only `adata.var`, never the expression matrix. Columns used:

- `ENST_ID` — Ensembl transcript ID
- `ENSG_ID` — Ensembl gene ID

### 2.3 `outputs/layer0/pseudobulk/counts_{cell_type}.csv` + `metadata_{cell_type}.csv` + `counts_{cell_type}_sensitivity.csv` + `metadata_{cell_type}_sensitivity.csv`

Both the primary pseudobulk files (AD + Control donors) and the sensitivity files (all conditions including Active Control donors) are inputs. M07 uses the primary files for per-donor PSI under the two tested conditions, and the sensitivity files for Active Control PSI of all transcripts (both significant and comparator). Surveyor reads only cell types present in the significant transcript list.

Not all cell types have sensitivity files (Lymphocyte does not). If a sensitivity file is absent for a cell type, Active Control PSI is set to `unavailable` and Visualization 3 shows only two bars (Control / AD).

### 2.4 `outputs/layer1/gene_level_de_results.csv` *(requires Layer 1 step15)*

Gene-level differential expression results. Columns: `gene_id`, `cell_type`, `log2FC`, `lfcSE`, `pval`, `padj`, `baseMean`. Used by M07 for expression pattern classification. If absent, this classification is marked `unavailable` — non-blocking.

### 2.5 Required Layer 1 addition: `step15_gene_de.R`

**This is the only required change to Layer 1.** A new step running DESeq2 gene-level DE on the same pseudobulk data, using the same primary model formula (`~ condition + age + sex + median_pct_mt`), processing each cell type in the existing loop. Runs after step14, before plots and the QC log. Output: `outputs/layer1/gene_level_de_results.csv`. Does not modify any existing step outputs.

### 2.6 Input Validation (pre-run)

Before any module starts, the orchestrator runs a validation pass:

| Check | Failure mode |
|---|---|
| `dtu_significant_all_celltypes.csv` exists and has required columns | Abort run with explicit error |
| h5ad exists and `adata.var` contains `ENST_ID` column | Abort run with explicit error |
| At least one pseudobulk CSV exists for each cell type in the significant list | Abort run for affected cell types; others proceed |
| Sensitivity CSVs: missing is non-fatal | Log warning; Active Control PSI marked unavailable |
| `gene_level_de_results.csv`: missing is non-fatal | Log warning; expression pattern marked unavailable |
| Significant CSV contains at least one row | Abort run — nothing to process |

Validation errors are written to the audit log before any API call is made. A run with blocking validation failures does not produce partial output for the failed candidates — it fails cleanly.

---

## 3. Output Contract

Surveyor produces outputs at three scopes: per (gene, cell_type), per gene (multi-cell-type only), and per run.

### 3.1 Per (gene, cell_type): structured JSON dossier

`output/candidates/{gene}_{cell_type}_dossier.json`

Machine-readable. Contains the full structured output of all modules. Consumed by Layer 3 (AlphaFold submission), Layer 4 (Discovery Studio preparation), Layer 5 (CDOCKER docking setup), and Layer 6 (repurposing screen). Complete and self-contained — downstream layers do not re-query databases.

Key top-level fields:

```
gene
cell_type
significant_transcripts         list of {enst_id, transcript_name, delta_psi, role}
comparator_transcript           {enst_id, transcript_name, source: "paired_hit" | "mane_select" | "highest_ctrl"}
tier                            1 | 2 | null (Type D)
splice_classification           A | B | C | D
domain_overlap                  list of affected domains with remapped residue coordinates
pdb_structures                  list of {pdb_id, method, resolution_A, chains} for comparator protein
docking_readiness               go | conditional | no-go
docking_readiness_reason        one-line string
drug_target_status              known | novel
ad_relevance                    direct | indirect | none
opentargets_score
opentargets_label               supported | emerging | novel | conflicting
ad_specificity_ratio
brain_specificity_index
expression_pattern              pure_switch | combined_de_dtu | complex | unavailable
robust_to_braak                 boolean
braak_condition_correlation     float (computed at runtime from metadata)
pathway_ad_flag                 boolean
interaction_ad_genes            list of strings
database_versions               {ensembl, uniprot, drugbank, chembl, opentargets, string, reactome, gtex}
parameters                      snapshot of config.yaml thresholds used in this run
input_file_checksums            {dtu_csv_md5, h5ad_md5, gene_de_csv_md5}
```

### 3.2 Per (gene, cell_type): HTML report

`output/reports/{gene}_{cell_type}_surveyor_report.html`

Self-contained single-file HTML with embedded CSS and inline SVG. No web server or external assets required. Content described in M08.

### 3.3 Per gene (multi-cell-type only): gene-level index page

`output/reports/{gene}_index.html`

Generated only when a gene is significant in more than one cell type. A lightweight single-page HTML file containing:

- Gene name, ENST IDs of all significant transcripts
- A table listing each significant cell type with its ΔPSI, padj, tier, splice classification, and docking readiness verdict — one row per (gene, cell_type) report
- The shared gene expression heatmap (Visualization 2 from M_VIS), embedded once here rather than duplicated in each per-cell-type report
- Links to each `{gene}_{cell_type}_surveyor_report.html`

This page is the entry point for reviewing a gene that was significant across multiple cell types.

### 3.4 Per (gene, cell_type): AlphaFold submission package

`output/alphafold/{gene}_{cell_type}/`

Generated only for `docking_readiness = go` or `conditional`. Contains:

- `{gene}_{enst_id}_ad_enriched.fasta` — protein sequence of the AD-enriched isoform
- `{gene}_{enst_id}_comparator.fasta` — comparator protein sequence (omitted if at least one PDB structure with resolution ≤ 3.5 Å exists in `pdb_structures`)
- `submission_params.json` — recommended ColabFold parameters (sequence length, MSA mode, template use)
- `docking_compounds/` — subdirectory containing one SDF file per docking candidate compound (3D conformers generated by M04; see §5 M04)

### 3.5 Global: run summary

`output/surveyor_run_summary.md`

One row per (gene, cell_type): gene, cell type, tier, splice classification, docking readiness verdict, one-line rationale. Readable in under five minutes.

### 3.6 Global: audit log

`output/surveyor_audit.log`

Timestamped log of every API call, cache hits vs. misses, validation warnings, runtime-computed values (e.g. Braak correlation), and database versions. Every dossier field is traceable to a specific log entry.

---

## 4. Module Architecture

Eight modules (M01–M08) plus an orchestrator. The orchestrator manages candidate iteration, gene deduplication, and parallel execution.

### Dependency graph

```
M01 (Sequence retrieval)  ──┐
                             ├──► M03 (Splice classification + tier) ──► M08 (Report assembly)
M02 (Domain + PDB lookup) ──┘                                                  ▲
                                                                                │
M04 (Drug target + 3D conformers)  ──────────────────────────────────────────► │
M05 (Disease association)  ──────────────────────────────────────────────────► │
M06 (Pathway + PPI)  ────────────────────────────────────────────────────────► │
M07 (Expression context)  ───────────────────────────────────────────────────► │
M_VIS (Visualizations)  ─────────────────────────────────────────────────────► │
```

M01 and M02 run sequentially first. M03 runs after both complete. M04, M05, M06, M07, and M_VIS run concurrently. M08 waits for all to finish.

### 4.1 Orchestrator specification

The orchestrator is the entry point for Surveyor. It performs the following steps in order:

**Step 1 — Input validation (§2.6).** Abort if blocking failures are found.

**Step 2 — Candidate table construction.** Read `dtu_significant_all_celltypes.csv`. Compute the Braak-condition correlation from the metadata files (Spearman correlation between numeric Braak stage and binary condition indicator across all donors with both values present). Store the correlation value in the audit log; embed it in each dossier's `braak_condition_correlation` field.

**Step 3 — Gene deduplication.** Build two work lists:
- `gene_work_list`: one entry per unique `gene_id`. Gene-level modules (M01, M02, M03, M04, M05, M06) are run once per entry, regardless of how many cell types the gene is significant in.
- `cell_type_work_list`: one entry per unique `(gene_id, cell_type)`. Cell-type-level modules (M07, M_VIS, M08) are run once per entry.

This prevents redundant API calls for genes significant in multiple cell types. The gene-level module outputs are stored in a shared in-memory cache keyed by `gene_id`, accessible to all subsequent cell-type-level runs for the same gene.

**Step 4 — Gene-level module loop.** For each gene in `gene_work_list`, run M01 → M02 → M03 sequentially, then submit M04, M05, M06 to the thread pool (max 4 concurrent threads across all genes). Wait for all gene-level work to complete.

**Step 5 — Cell-type-level module loop.** For each (gene, cell_type) in `cell_type_work_list`, run M07, then submit M_VIS to the thread pool. When both complete, run M08 to assemble the JSON dossier and HTML report for that (gene, cell_type).

**Step 6 — Multi-cell-type index pages.** After all (gene, cell_type) pairs are complete, identify genes with more than one significant cell type and generate the gene-level index HTML (§3.3) for each.

**Step 7 — Global outputs.** Write `surveyor_run_summary.md` and finalize `surveyor_audit.log`.

---

## 5. Module Specifications

---

### M01 — Sequence Retrieval

**Database:** Ensembl REST API  
**Fallback:** BioMart bulk query if REST rate limit is hit repeatedly

#### Comparator selection

M01 first determines which transcripts need sequences, using all significant transcripts for the gene across all cell types (not just the current cell type). Comparator selection is gene-level.

**Case: exactly two significant transcripts with opposite ΔPSI signs**

The transcript with positive ΔPSI is the AD-enriched isoform; the one with negative ΔPSI is the control-enriched isoform. Each serves as the other's structural comparator. This is the cleanest case.

**Case: one significant transcript**

The significant transcript plus a comparator are retrieved. Comparator selected by priority:

1. **MANE Select** — retrieved from Ensembl. Correct baseline: represents the transcript used in most existing structural data and drug target annotations.
2. **Highest-expressed in control donors** — if no MANE Select exists, Surveyor selects the transcript with the highest mean count across control donors in the most-significant cell type, using Layer 0 pseudobulk data, excluding the significant transcript itself.
3. **Error** — if neither is available, M01 logs a blocking failure for this gene.

**Case: two significant transcripts with the same ΔPSI sign**

Both transcripts are gaining in AD (or both losing) simultaneously, implying a more complex splicing rearrangement rather than a simple isoform switch. In this case, treat each significant transcript as if it were a single-hit: retrieve MANE Select as the structural reference for each, independently. Log a warning in the audit log describing the same-sign situation. Both transcripts proceed through M02–M08 with MANE Select as their comparator.

**Case: three or more significant transcripts**

MANE Select is the single structural reference for all. Each significant transcript is compared against MANE Select. Each transcript is labeled by its ΔPSI sign (AD-enriched or control-enriched). The dossier's `significant_transcripts` list contains all of them; the `comparator_transcript` field is set to the MANE Select entry with `source: "mane_select"`. Log a note in the audit log.

#### What it retrieves

For every ENST ID in the candidate set (significant transcripts + comparator):

- Full protein sequence (FASTA)
- Exon structure: genomic start, end, CDS start phase, CDS end phase, exon rank per transcript
- Transcript biotype
- Protein length (aa); CDS length (nt)
- UniProt accession — via Ensembl cross-reference endpoint (`/xrefs/id/{ENST_ID}`), filtered to UniProt/Swiss-Prot. Passed to M02.

#### Key algorithmic decisions

**Exon comparison: coordinate-based symmetric difference**

Exon differences are identified by genomic coordinate overlap (≥ 80% of shorter exon's length), not sequence alignment. Unambiguous for directly resolved long-read transcript structures; robust to minor annotation-version differences.

**CDS phase recording**

M01 retrieves the full protein sequence for both the significant transcript and its comparator and compares them directly. This is the authoritative protein-level signal (`protein_identical`, `protein_length_delta`), and it feeds M03's Type N gate. The earlier v2.1 exon-length frame-shift proxy (mod-3 on unique exon lengths) was **dropped in v2.2**: it false-positives on UTR exons (e.g. EIF4G2-208, where a unique non-coding exon flagged a "frame-shift" while the encoded protein is identical to the comparator). Exon coordinates are still reported for the M08 architecture schematic, but the frame/coding consequence is read from the protein comparison, not exon arithmetic.

**Biotype check**

`retained_intron` or `nonsense_mediated_decay` biotype → flagged as potentially non-protein-coding, passed to M08. Surveyor does not discard such isoforms.

#### Output

Per-transcript:
- FASTA protein sequence
- Exon structure table (TSV)
- Exon symmetric difference table relative to comparator
- Biotype and frame-shift flags
- UniProt accession (passed to M02)

---

### M02 — Domain Annotation and PDB Lookup

**Databases:** UniProt REST API (primary), InterPro API (supplementary)

#### What it retrieves from UniProt

For the canonical UniProt accession of each candidate gene:

- All annotated sequence features: DOMAIN, BINDING, ACT_SITE, MOTIF, MOD_RES, REGION, DISULFID, SIGNAL, TRANSIT
- For each feature: type, description, start residue, end residue, evidence type (experimental / by similarity / predicted)
- Active site residues with catalytic roles where annotated
- Known disease variants overlapping functional features
- Reviewed status (Swiss-Prot reviewed vs. TrEMBL unreviewed)
- **PDB cross-references**: all PDB entries listed in the UniProt entry's "3D structure" section, with method (X-ray crystallography, cryo-EM, NMR), resolution in Å, and chain identifiers

The PDB cross-references are retrieved from the same UniProt API call as the domain annotations — no additional API call is needed. The list is filtered to structures of the comparator isoform (or MANE Select) only. Structures of the AD-enriched isoform are noted separately if they exist.

#### PDB availability output

For each PDB entry:
- PDB ID, experimental method, resolution (Å), deposited chains
- Flag: `resolution_ok` = true if resolution ≤ 3.5 Å (threshold for reliable docking template)
- Flag: `covers_splice_region` = true if the structure's sequence coverage includes residues encoded by the alternatively spliced exon (checked against M01 exon coordinates)

This output populates `pdb_structures` in the dossier JSON and drives:
- M08 checklist item: `✓/✗/⚠ Canonical structure: PDB {ID} available / AlphaFold required / low-confidence region`
- §3.4 AlphaFold package: comparator FASTA omitted if at least one `resolution_ok` structure exists

#### Key algorithmic decisions

**Residue coordinate remapping**

UniProt domain annotations are defined on the canonical UniProt sequence. M02 remaps domain boundaries to AD-enriched isoform coordinates via pairwise global alignment (BLOSUM62, gap open −10, gap extend −0.5). When a domain feature maps to a region deleted in the AD-enriched isoform, it is recorded as partially or fully absent.

**Evidence type weighting**

Propagated to M03:
- Experimental binding site overlapping splice event → confirmed Type A
- "By similarity" binding site → probable Type A, flagged for manual review
- Predicted → Type A candidate, lower confidence

**PTM site proximity**

PTM sites (phosphorylation, ubiquitination, acetylation) within 5 residues of the splice event boundary are flagged separately. Appears in the M08 report regardless of splice classification.

#### Output

Per-gene:
- Full feature table with original UniProt coordinates
- Remapped feature table in AD-enriched isoform coordinates
- Domain absence/truncation flags
- PTM proximity flags
- PDB structure list with resolution and splice-region coverage flags

---

### M03 — Splice Event Classification and Tier Assignment

**Dependencies:** M01, M02  
**External databases:** None

#### Pre-classification gate: protein-level change (Type N)

**Added in v2.2 after empirical validation against the 12 hits.** Before any splice classification, M03 checks M01's protein-sequence comparison between the significant transcript and its comparator. If the two encode an **identical protein**, or the significant transcript is **non-coding** (no protein product, e.g. `retained_intron`), there is no structural docking hypothesis — the DTU is a regulatory event (alternative UTR / non-coding exon usage / NMD), not a protein-altering one.

**Type N — No protein-level change**

The significant transcript produces a protein identical to its comparator, or produces no protein at all. The isoform switch is real and may be biologically important (UTR-mediated translational control, transcript stability, NMD-coupled regulation, localization), but a canonical same-pocket docking comparison is meaningless because the protein is unchanged.

Priority for docking: none. No tier assigned. Routed to the **regulatory analysis path**: the candidate still receives a full dossier (M04–M07 run normally — drug target status, disease association, pathways, expression are all still relevant), but the AlphaFold package is **not** generated and the verdict is NO-GO for docking with reason "no protein-level change — regulatory candidate."

This gate is checked first. A candidate that fails it (Type N) skips Type A/B/C/D classification entirely. In the current dataset 5 of 12 hits are Type N: CAMK2B-204/-216 (identical 542 aa), CSRNP3-202, EIF4G2-208 (identical to MANE Select), GRIA2-231 (`retained_intron`, no protein).

#### Functional feature set

**Broadened in v2.2.** UniProt does not annotate every protein's functional core as `DOMAIN`/`BINDING`/`ACT_SITE`. Empirically among the candidates: RHOT1 uses `Domain` + `Binding site` (classic), KLC1 uses `Repeat` + `Region` (TPR cargo domain, no `Domain` feature), MGRN1 uses `Zinc finger` + `Motif` (RING E3 ligase, no `Domain` or `Binding site`). Keying only on DOMAIN/BINDING/ACT_SITE would misclassify KLC1 and MGRN1 as Type D.

M03 therefore uses two feature tiers:

- **Pocket features** (drive Type A/B): `Binding site`, `Active site`, `Site`
- **Structural-domain features** (drive Type C, and Type B by containment): `Domain`, `Repeat`, `Zinc finger`, `DNA binding`, `Coiled coil`, `Motif`, `Region` (the last two with lower confidence, recorded as such)

#### Classification scheme

**Type A — Splice event directly overlaps binding pocket**

Alternatively spliced exon encodes residues annotated as a pocket feature (Binding site / Active site / Site), or the splice event removes an entire structural-domain feature (distance = 0, auto-Type A). AD-enriched isoform has a structurally different drug binding pocket.

Priority for docking: highest. Proceed regardless of other factors.

**Type B — Splice event adjacent to binding pocket**

Alternatively spliced exon does not directly encode pocket residues but falls within 15 residues of a pocket feature, OR the exon encodes part of a structural-domain feature that contains a pocket feature even if the pocket residues are retained.

Priority: high. Proceed with note that binding affinity change may be indirect.

**Type C — Splice event within a domain, distal from binding pocket**

Alternatively spliced exon encodes part of a structural-domain feature but the nearest pocket feature is > 15 residues away (or the protein has no annotated pocket feature at all). Pocket residues, where they exist, are conserved in both isoforms.

Priority: medium. Verdict is always at least CONDITIONAL (Type C = Tier 2 by definition).

**Type D — Splice event in unannotated or disordered region**

No UniProt functional annotation (neither pocket nor structural-domain feature) in the alternatively spliced region. May affect PPI surfaces, localization signals, or degradation signals, but canonical docking comparison unlikely to produce meaningful ΔΔG.

Priority: low. Flag for PPI surface analysis. No tier assigned.

#### Tier assignment

Derived from splice classification — not read from any external file.

- **Tier 1:** Type A or B
- **Tier 2:** Type C
- **Type D:** no tier; candidate routed to PPI surface analysis path
- **Type N:** no tier; candidate routed to regulatory analysis path (no docking, no AlphaFold)

#### Distance calculation

Minimum residue distance in AD-enriched isoform coordinate system (after M02 remapping), from the closest residue of the alternatively spliced exon to the closest annotated binding residue.

#### Output

Per-gene:
- Type A / B / C / D classification
- Tier 1 / 2 / null
- Specific exon(s) with genomic and protein coordinates
- Distance to nearest binding feature (residues)
- Confidence level (from UniProt evidence type)
- Narrative justification (one sentence, included verbatim in M08)

---

### M04 — Drug Target, Known Compounds, and 3D Conformer Generation

**Databases:** DrugBank (primary), ChEMBL  
**Library:** RDKit (3D conformer generation — local, no API call)

*Note: DrugBank access method (academic API key vs. monthly XML dump) is unresolved — see §10.*

#### What it retrieves from DrugBank

- All drugs targeting the gene, by group: approved, investigational, experimental, withdrawn
- For each drug: name, DrugBank ID, mechanism of action, pharmacological action, indication(s), approval status
- Known binding site on the protein target where annotated
- Drug-drug interaction liabilities

#### What it retrieves from ChEMBL

- All bioactivity records: Kd, IC50, Ki, EC50 with assay conditions
- Assay type: binding vs. functional
- Compound SMILES for Kd or Ki < 1 µM
- Clinical phase of each compound

#### 3D conformer generation (post-retrieval, local)

For each docking candidate compound (Kd or Ki < 1 µM), M04 generates a 3D conformer using RDKit:

1. Parse SMILES → RDKit molecule object
2. Add hydrogens (`Chem.AddHs`)
3. Generate up to 100 conformers using ETKDG (`AllChem.EmbedMultipleConfs`, `ETKDGv3`)
4. Minimize each conformer with MMFF94 force field (`AllChem.MMFFOptimizeMoleculeConfs`)
5. Select the lowest-energy conformer
6. Write to SDF format: `output/alphafold/{gene}_{cell_type}/docking_compounds/{compound_id}.sdf`

If RDKit fails to embed a conformer (typically due to complex ring systems), the SMILES is preserved in the dossier with a flag `conformer_generation_failed: true`. Layer 5 is responsible for handling these manually.

This step runs locally with no API call and does not block any parallel modules.

#### Key algorithmic decisions

**AD relevance — two-tier approach**

*Direct (Tier 1):* Drug indication in DrugBank contains "Alzheimer," "dementia," or a related MeSH term.

*Indirect (Tier 2 — repurposing candidates):* Drug targets a pathway cross-referenced with M06 Reactome output that is directly implicated in AD pathology. AD-relevant pathway list: amyloid precursor protein processing, tau phosphorylation and aggregation, mitophagy, endosomal-lysosomal trafficking, neuroinflammation, axonal transport, synaptic plasticity and LTP.

**ChEMBL potency threshold**

Kd or Ki < 1 µM → docking candidate (specific binder). IC50-only compounds included but de-prioritized (functional inhibition ≠ direct pocket binding).

**Binding site cross-reference with M03**

Where DrugBank provides a binding site annotation, M04 cross-references it against the M03 splice classification. Drug binding at a site altered by the splice event (Type A or B) is flagged as a high-priority docking pair.

#### Output

Per-gene:
- Drug table: name, group, indication, AD relevance tier, mechanism, phase
- Bioactivity table: compound, assay type, Kd/IC50/Ki, ChEMBL assay ID
- Docking candidate list: SMILES + compound ID for all sub-µM binders
- AD relevance summary: n_direct, n_indirect, n_novel
- SDF files: one per successfully embedded docking candidate (written to §3.4 directory)

---

### M05 — Disease Association

**Database:** OpenTargets Platform API (GraphQL)

#### What it retrieves

For each candidate gene against Alzheimer's disease (EFO:0000249):

- Overall association score (0–1)
- Evidence breakdown: genetic association (GWAS, rare variant), somatic mutation, differential expression, animal model, known drug, literature
- For genetic evidence: lead SNPs, effect allele frequency, GWAS p-value, population
- For context: scores against Parkinson's disease (EFO:0002508), frontotemporal dementia (EFO:0001627), ALS (EFO:0000253)

#### Key algorithmic decisions

**Score labels — no hard threshold as pass/fail**

- ≥ 0.5 → `supported target`
- 0.1–0.49 → `emerging target`
- < 0.1, no genetic evidence → `novel target — low prior evidence`
- < 0.1, contradicting evidence → `conflicting evidence — review manually`

Understudied genes (MGRN1, RHOT1) will have low scores not because evidence contradicts their AD relevance, but because they are underinvestigated. A cutoff would exclude exactly the novel findings this pipeline is designed to surface.

**AD specificity ratio**

```
AD_specificity = opentargets_score_AD / mean(score_PD, score_FTD, score_ALS)
```

Ratio > 2: AD-specific. Near 1: general neurodegeneration.

**GWAS SNP splice site check**

M05 checks whether any lead GWAS SNPs fall within exon boundaries or splice sites (using M01 exon coordinates). A GWAS hit within a splice site is flagged explicitly — direct genetic evidence that splicing at this locus is disease-relevant.

#### Output

Per-gene:
- OpenTargets score and evidence breakdown
- AD target label
- AD specificity ratio
- GWAS SNP table with splice site proximity flag
- Clinical trial evidence (cross-referenced with M04)

---

### M06 — Pathway and Interaction Analysis

**Databases:** STRING API, Reactome Content Service

#### What it retrieves from STRING

- Top interaction partners, combined score ≥ 0.7 (high confidence)
- Evidence channels: experimental, co-expression, database, co-occurrence, text-mining, homology
- Functional enrichment of interaction network: GO Biological Process (FDR < 0.05), KEGG, Reactome

#### What it retrieves from Reactome

- All pathways containing the gene, with full hierarchy
- AD-relevant pathway membership (pathway list from M04 section)
- Whether the gene is in the same Reactome pathway as any known AD drug target: APP, BACE1, PSEN1, MAPT, APOE, BIN1, CLU, CR1, PICALM

#### Key algorithmic decisions

**STRING threshold — 0.7 primary, 0.4 secondary**

At 0.7, interactions are predominantly experimentally backed. At 0.4, text-mining dominates. Both are queried and reported separately; never merged. Enrichment analysis run only on the high-confidence (0.7) network. If the high-confidence network has fewer than 5 partners, a note is added that enrichment results may be unreliable due to small network size.

**AD-gene interaction flag**

Top 20 high-confidence partners cross-referenced against curated AD-related protein list:
- Core genetic risk: APP, PSEN1, PSEN2, APOE, TREM2, BIN1, CLU, CR1, PICALM, ABCA7
- Core pathological: MAPT, SNCA
- Major AD drug targets: BACE1, ADAM10, GSK3B, CDK5
- Gamma-secretase complex: ADAM17, NCSTN, APH1A, PSENEN

Direct high-confidence interaction with any of these is flagged with the AD gene name and evidence type.

**Reactome hierarchy traversal**

Full hierarchy recorded: specific sub-pathway and top-level category both reported.

#### Output

Per-gene:
- High-confidence interaction table (top 20, score ≥ 0.7) with evidence channels
- Moderate-confidence interactions (0.4–0.69) labeled separately
- AD-gene interaction flags
- Reactome pathway membership with full hierarchy
- AD-relevant pathway flags
- GO Biological Process enrichment (top 10 by FDR) — only reported if ≥ 5 high-confidence partners

---

### M07 — Expression Context

**Primary data:** Layer 0 pseudobulk CSVs (primary + sensitivity)  
**Supplementary:** GTEx Portal API

#### Per-donor PSI computation

M07 reads both `counts_{cell_type}.csv` and `counts_{cell_type}_sensitivity.csv`. From the primary CSV, per-donor PSI is computed for AD and Control conditions. From the sensitivity CSV, Active Control PSI is computed.

For the significant transcript(s): PSI values for all three conditions come directly from the significant CSV (`psi_AD`, `psi_ctrl`, `psi_active_ctrl`). M07 uses the pseudobulk CSVs to get **per-donor** values, not just condition means.

For the comparator transcript: per-donor PSI for AD and Control is computed from the primary CSV; Active Control PSI from the sensitivity CSV.

**Comparator absent from pseudobulk CSV:** If the comparator transcript is not in the pseudobulk CSV (e.g., MANE Select below Layer 0 prevalence threshold), M07 handles this explicitly:
- The per-donor strip plot (Vis 1) shows the significant transcript's PSI only, labeled accordingly
- The stacked bar (Vis 3) shows the significant transcript's PSI plus an "other isoforms" segment computed as 1 − PSI of the significant transcript, with a note: "Comparator transcript below prevalence threshold — not individually shown"
- The dossier JSON records `comparator_in_pseudobulk: false`

#### Multi-cell-type handling

Separate PSI matrices computed per significant cell type. Stored under distinct keys in the dossier JSON. Each cell type produces its own visualization set. No cross-cell-type aggregation.

#### Gene-level expression context

Requires `gene_level_de_results.csv` (step15). If available:

- **Pure isoform switch:** |log2FC| < 0.5 AND padj > 0.05 (gene-level) but significant DTU. Total gene output constant; only isoform proportion shifts. Cleanest signal for the drug target hypothesis.
- **Combined DE + DTU:** Both gene-level expression and isoform proportion change in AD.
- **Complex:** Gene expression and isoform proportion shift in opposite directions.

If `gene_level_de_results.csv` is absent: `expression_pattern: unavailable`.

#### Cell-type expression breadth

Identifies which cell types express the gene above mean CP10K ≥ 1 across donors in any condition, using Layer 0 pseudobulk data.

#### GTEx retrieval

Median TPM queried for:
- All tissues (for BSI denominator)
- Brain regions: frontal cortex BA9, hippocampus, anterior cingulate BA24, cortex, putamen, caudate, cerebellum
- Cardiac: left ventricle, atrial appendage
- Hepatic: liver

**Brain Specificity Index:**

```
BSI = mean(brain region TPM) / mean(all tissue TPM)
```

- BSI > 5: brain-specific
- BSI 2–5: brain-enriched
- BSI 1–2: ubiquitous
- BSI < 1: peripherally enriched

Limitation noted in report: GTEx is bulk tissue, no isoform resolution. BSI is gene-level only.

#### Output

Per-gene per-significant-cell-type:
- Per-donor PSI matrix
- Per-condition PSI summary (mean ± SD)
- `comparator_in_pseudobulk` flag
- Expression pattern classification or `unavailable`
- Cell-type expression breadth map
- GTEx brain specificity index
- Cardiac and hepatic expression flags
- Visualization datasets consumed by M_VIS

---

### M_VIS — Visualization Generation

**Dependencies:** M07 output, M01 output  
**Libraries:** matplotlib, seaborn (SVG output)

One complete set of three visualizations per significant cell type. When a gene is significant in two cell types, six SVGs are produced. The shared gene expression heatmap (Vis 2) is produced once per gene, not per cell type, and embedded in the gene-level index page (§3.3) rather than duplicated in each per-cell-type report.

#### Visualization 1 — Per-donor PSI strip plot (per significant cell type)

X-axis: condition (Control / Active Control / AD). Y-axis: PSI of the AD-enriched isoform (0–1). Each donor is one point, colored by condition. Mean ± SD overlaid. Donor IDs on hover via SVG title elements. If Active Control PSI is unavailable (no sensitivity CSV), the strip plot shows two groups only.

#### Visualization 2 — Gene expression heatmap (per gene, shared)

X-axis: all cell types. Y-axis: all donors, grouped and sorted by condition. Color: CP10K normalized gene expression (log scale, capped at 99th percentile). Columns with significant DTU have a distinct border. Embedded in the gene-level index page for multi-cell-type genes; embedded in the per-cell-type HTML report for single-cell-type genes.

#### Visualization 3 — Isoform proportion stacked bar (per significant cell type)

Three bars: Control / Active Control / AD (two bars if Active Control unavailable). Stacked to 100%. AD-enriched isoform distinctly colored. Minor isoforms (PSI < 5% in all conditions) collapsed to "other." If comparator is not in the pseudobulk CSV, "other isoforms" segment labeled accordingly with a note.

---

### M08 — Report Assembly

**Dependencies:** All modules  
**Libraries:** Jinja2, json, pandas

M08 waits for all modules to complete, then assembles the JSON dossier and HTML report for one (gene, cell_type).

#### JSON dossier assembly

Structured merge of all module outputs. Every field has a defined key, type, and required/optional status. Fields from failed non-blocking modules are set to `null` with a `_status` field. Downstream layers raise explicit errors on null required fields.

#### HTML report structure

**Header:** Gene, cell type(s) of significance, ΔPSI, padj, tier, run date, database versions.

**Executive summary:** Three auto-generated sentences — what this candidate is, why it was selected, what the docking hypothesis is.

**Docking readiness checklist:**

```
✓     Protein sequence retrieved — {n} aa (AD-enriched), {n} aa (comparator)
✓/✗   Biotype: protein-coding confirmed / WARNING: potential NMD/retained-intron isoform
✓/✗   Protein-level change: differs from comparator ({±n} aa) / TYPE-N: protein identical — no docking hypothesis
✓     Exon diff computed — {n} exon(s) unique to AD-enriched isoform
✓     Splice classification: Type {A/B/C/D} — {one-line rationale from M03}
✓/✗/⚠ Canonical structure: PDB {ID} ({method}, {resolution} Å) / AlphaFold required / only low-resolution structures (> 3.5 Å)
✓/✗   Known drugs: {n} direct AD / {n} indirect / none
✓/⚠   OpenTargets: {score} ({label})
✓/⚠   Brain specificity: BSI = {value} ({interpretation})
⚠     robust_to_braak = FALSE — see note [Braak-condition r = {runtime value}]
→ VERDICT: PROCEED / CONDITIONAL / NO-GO
  Reason: {one line}
```

**Braak robustness note (rendered in checklist):**

> `robust_to_braak = FALSE` for this candidate. The Braak sensitivity model (which adds braak_stage to the primary formula) produces an underdetermined design matrix for this dataset because Braak stage and diagnosis condition are highly correlated (Spearman r = {braak_condition_correlation}, computed from this dataset's metadata). The hit disappearing under the Braak model reflects collinearity instability, not evidence that Braak stage is confounding the DTU result. This should be stated explicitly in the manuscript methods. It does not affect the docking readiness verdict.

The `{braak_condition_correlation}` value is the runtime-computed Spearman r stored in the dossier JSON, not a hardcoded constant.

**Verdict logic — explicit priority order:**

Rules are evaluated top to bottom. The first matching rule determines the verdict. A lower rule never upgrades a verdict assigned by a higher rule.

| Priority | Condition | Verdict |
|---|---|---|
| 1 | M01 blocking failure (no sequence) | NO-GO — sequence unavailable |
| 2 | **Type N (protein identical to comparator, or non-coding / no protein)** | **NO-GO (docking) — no protein-level change; routed to regulatory analysis path. Full dossier produced; no AlphaFold package.** |
| 3 | M02 blocking failure (no domain data) | CONDITIONAL — manual splice classification required |
| 4 | Type A or B AND biotype confirmed protein-coding AND ≥ 1 known drug (direct or indirect) | PROCEED |
| 5 | Type A or B AND biotype confirmed protein-coding AND no known drugs AND BSI ≥ 1 | CONDITIONAL — no known drugs, target is novel |
| 6 | Type A or B AND biotype confirmed protein-coding AND BSI < 1 | CONDITIONAL — off-target risk (peripherally enriched) |
| 7 | Type C | CONDITIONAL — domain-level change, indirect docking signal |
| 8 | Type D AND known drug binding in alternatively spliced region | CONDITIONAL — drug-specific override |
| 9 | Type D AND no known drugs AND BSI < 1 | NO-GO — no docking hypothesis; flag for PPI surface analysis |
| 10 | Type D AND no known drugs AND BSI ≥ 1 | CONDITIONAL — PPI surface analysis recommended |

Type N is checked immediately after the M01 sequence-availability gate and before everything else: a protein that does not change cannot have a docking hypothesis regardless of drug or disease evidence. The 5 current Type N hits (CAMK2B-204/-216, CSRNP3-202, EIF4G2-208, GRIA2-231) take this row.

All current 12 candidates have `robust_to_braak = FALSE`. This does not appear in the verdict logic — it is informational only.

**Section 1 — Isoform architecture:** Exon structure schematic, domain map with splice region highlighted, Type classification and distance to nearest binding feature, PDB availability.

**Section 2 — Drug target profile:** DrugBank table, ChEMBL bioactivity table, AD relevance classification, docking compound count.

**Section 3 — Disease association:** OpenTargets score panel, evidence breakdown, GWAS SNPs, AD specificity ratio.

**Section 4 — Biological context:** STRING interaction summary (top 10 partners, AD-gene flags), Reactome pathway table, GO enrichment (if ≥ 5 high-confidence partners).

**Section 5 — Expression profile (one panel per significant cell type):**

One collapsible panel per cell type. Each panel contains its three SVG visualizations (Vis 1 and Vis 3 are cell-type-specific), PSI summary table, and expression pattern classification. For single-cell-type genes, Vis 2 (gene expression heatmap) appears at the top of this section. For multi-cell-type genes, Vis 2 is omitted here and shown instead in the gene-level index page.

**Section 6 — AlphaFold submission:** ColabFold recommended parameters, PDB recommendation for comparator, docking compound SDF file list.

**Footer:** All API calls, cache status, database versions, input file checksums, Braak-condition correlation value, run timestamp.

---

## 6. Caching Strategy

All external API responses cached locally. Cache key = MD5 hash of full request URL + query parameters.

| Database | TTL | Rationale |
|---|---|---|
| Ensembl REST | 90 days | Quarterly releases |
| UniProt | 90 days | Swiss-Prot manually curated; infrequent changes |
| DrugBank | 14 days | Approval statuses change monthly |
| ChEMBL | 30 days | New bioactivity data deposited regularly |
| OpenTargets | 30 days | Updated with new GWAS data |
| STRING | 90 days | Rebuilt with major releases only |
| Reactome | 90 days | Pathway annotations change slowly |
| GTEx | 180 days | v10 is a versioned static dataset |

`--force-refresh` flag bypasses all cache. Version mismatch between cached response and expected version in `config.yaml` discards the cache entry and triggers a fresh query.

---

## 7. Error Handling and Partial Failure Policy

**Non-blocking (pipeline continues, field marked unavailable):**
- M04 DrugBank or ChEMBL failure → drug section `unavailable — API error`; checklist uses `drug data unavailable`
- M04 RDKit conformer generation failure → SMILES preserved with `conformer_generation_failed: true`
- M05 OpenTargets failure → `unavailable`, neutral label
- M06 STRING or Reactome failure → pathway section `unavailable`
- M06 < 5 high-confidence partners → GO enrichment `unavailable — network too small`
- M07 GTEx failure → BSI `unavailable`
- M07 gene-level DE file absent → expression pattern `unavailable`
- M07 sensitivity CSV absent → Active Control PSI `unavailable`; visualizations show two conditions only
- M07 comparator not in pseudobulk CSV → `comparator_in_pseudobulk: false`; visualizations degrade gracefully

**Blocking (candidate halts, logged in run summary):**
- M01 Ensembl failure → `failed — sequence retrieval error`
- M01 comparator selection failure (no MANE Select, no eligible control transcript) → `failed — no comparator available`
- M02 UniProt failure → M03 cannot classify; M04–M07 still run; verdict forced to CONDITIONAL

All errors logged to `surveyor_audit.log` with full error message, HTTP status code, and timestamp.

---

## 8. Reproducibility Requirements

1. **Database version recording:** Every API response's version metadata recorded in audit log and in `dossier.database_versions`.
2. **Parameter recording:** All thresholds in `dossier.parameters` as a snapshot of `config.yaml` at run time.
3. **Input file checksums:** MD5 of all input files recorded in `dossier.input_file_checksums` at run start.
4. **Runtime-computed values:** Braak-condition Spearman r recorded in audit log and in `dossier.braak_condition_correlation`.
5. **Deterministic output:** Same inputs + same cache + same `config.yaml` → identical output. No random or stochastic steps (RDKit conformer generation uses a fixed random seed specified in `config.yaml`).

---

## 9. Configuration File Structure

```yaml
version: "2.1"

# Input paths (relative to pipeline root)
inputs:
  dtu_significant_csv:        "outputs/layer1/dtu_significant_all_celltypes.csv"
  h5ad_path:                  "outputs/layer0/filtered_adata/adata_tx_step03.h5ad"
  pseudobulk_dir:             "outputs/layer0/pseudobulk"
  gene_de_csv:                "outputs/layer1/gene_level_de_results.csv"  # optional

# API endpoints
apis:
  ensembl:     "https://rest.ensembl.org"
  uniprot:     "https://rest.uniprot.org/uniprotkb"
  drugbank:    "https://api.drugbank.com/v1"
  chembl:      "https://www.ebi.ac.uk/chembl/api/data"
  opentargets: "https://api.platform.opentargets.org/api/v4/graphql"
  string:      "https://string-db.org/api"
  reactome:    "https://reactome.org/ContentService"
  gtex:        "https://gtexportal.org/api/v2"

# API authentication
auth:
  drugbank_api_key: "${DRUGBANK_API_KEY}"

# Thresholds
thresholds:
  string_confidence_high:            0.7
  string_confidence_moderate:        0.4
  string_min_partners_for_enrichment: 5
  chembl_potency_uM:                 1.0
  pdb_resolution_max_angstrom:       3.5
  splice_domain_proximity_residues:  15
  exon_overlap_min_fraction:         0.80
  gtex_expression_minimum_tpm:       1.0
  brain_specificity_enriched:        2.0
  brain_specificity_specific:        5.0
  cell_type_expression_min_cp10k:    1.0
  gene_de_log2fc_threshold:          0.5
  gene_de_padj_threshold:            0.05

# RDKit conformer generation
rdkit:
  random_seed:        42
  n_conformers:       100
  force_field:        "MMFF94"

# Tier assignment (derived in M03)
tier_assignment:
  tier1: [A, B]
  tier2: [C]

# AD pathway list (M04, M06)
ad_pathways:
  - "amyloid precursor protein processing"
  - "tau protein phosphorylation"
  - "mitophagy"
  - "endosomal-lysosomal trafficking"
  - "neuroinflammation"
  - "axonal transport"
  - "synaptic plasticity"
  - "long-term potentiation"

# AD gene list (M06 interaction flagging)
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
  - ADAM17
  - NCSTN
  - APH1A
  - PSENEN

# Braak robustness
braak_robustness:
  affects_verdict:  false
  checklist_symbol: "⚠"

# Caching
cache:
  directory: "cache"
  ttl_days:
    ensembl:     90
    uniprot:     90
    drugbank:    14
    chembl:      30
    opentargets: 30
    string:      90
    reactome:    90
    gtex:        180
  force_refresh: false

# Concurrency
execution:
  max_parallel_modules:    4
  request_timeout_seconds: 30
  max_retries:             3
  retry_backoff_seconds:   5

# Output
output:
  candidates_dir:  "output/candidates"
  reports_dir:     "output/reports"
  alphafold_dir:   "output/alphafold"
  run_summary:     "output/surveyor_run_summary.md"
  audit_log:       "output/surveyor_audit.log"
```

---

## 10. Open Questions Requiring Resolution Before Implementation

Two questions remain unresolved.

**Q1 — DrugBank access method**

Academic API key or monthly XML bulk download? API is convenient for per-compound automated queries. XML dump is more complete and not rate-limited. Impact: M04 implementation differs substantially. With the XML dump, M04 becomes a local parsing module rather than a REST client, and the DrugBank cache TTL becomes irrelevant.

**Q2 — Conda environment specification**

Python version and conda environment name on node212. Required packages: `requests`, `anndata`, `pandas`, `numpy`, `biopython` (pairwise alignment in M02), `jinja2`, `matplotlib`, `seaborn`, `pyyaml`, `rdkit`. Package versions should be pinned for reproducibility. Impact: environment setup and `requirements.txt` / `environment.yml`.

---

## 11. Relationship to Adjacent Layers

| Layer | Direction | Interface |
|---|---|---|
| VISCACHA Layer 0 | Upstream → Surveyor | `adata_tx_step03.h5ad`, pseudobulk CSVs (primary + sensitivity) |
| VISCACHA Layer 1 | Upstream → Surveyor | `dtu_significant_all_celltypes.csv`, `gene_level_de_results.csv` |
| AlphaFold (Layer 3) | Surveyor → Downstream | `output/alphafold/{gene}_{cell_type}/` (FASTA + submission params) |
| Discovery Studio (Layer 4) | Surveyor → Downstream | `dossier.json` (binding site coordinates, domain map, PDB IDs) |
| CDOCKER (Layer 5) | Surveyor → Downstream | `dossier.json` + `docking_compounds/*.sdf` (3D conformers ready for docking) |
| Repurposing screen (Layer 6) | Surveyor → Downstream | ChEMBL potency table + docking candidate SMILES from `dossier.json` |

Any change to the dossier JSON schema must be versioned and communicated to downstream layer implementations. Downstream layers read specific fields and raise explicit errors on null required fields.

---

*End of Surveyor Layer 02 Conceptual Plan v2.1-draft*

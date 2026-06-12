# Viscacha-v01
**AD drug target isoform divergence — full analysis pipeline**

---

## Overview

Viscacha is a serial pipeline that tests whether Alzheimer's disease drug target genes express a non-canonical dominant isoform in diseased brain tissue, whether that isoform is structurally divergent at the drug binding site, and whether existing AD drugs bind it with reduced affinity.

The pipeline runs in two stages in series:

- **Pipeline 1** (Layers 0–2): Statistical DTU identification using pseudo-bulk DRIMSeq + stageR
- **Pipeline 2** (Layers 3–4): Functional annotation and structural pharmacology via IsoformSwitchAnalyzeR + AlphaFold2 + docking

A parallel sensitivity check using Allos runs on top hits only and is reported in supplementary.

**Central hypothesis:** AD drug target genes show significant differential transcript usage (DTU) between AD and control post-mortem brain tissue, and the dominant disease-state isoform encodes sequence differences that reduce drug binding affinity (ΔΔG > 1.5 kcal/mol) relative to the canonical isoform used in drug development.

---

## Data inputs

| Object | Dimensions | Contents |
|---|---|---|
| `adata_sr.h5ad` | 138,385 cells × 34,903 genes | Short-read reference. Cell type labels, UMAP, doublet scores, `pct_counts_mt`, `total_counts`, Leiden clusters |
| `adata_transcript_loose_filtering_for_bulk_analysis.h5ad` | 95,487 cells × 63,994 transcripts | Transcript-level pseudo-bulk source. 75 transcripts/cell average. `sample`, `donor`, `condition`, `cell_type`, `leiden` |
| `삼성서울_알츠하이머_샘플_list.xlsx` | 18 SMC donors | Age, sex, A/B/C ADNC scores, Braak stage, Thal phase, CERAD, clinical Dx, primary Dx |
| `김준표교수님_분양_건_임상정보.pptx` | 9 PO donors | Age, sex, APOE genotype, ADNC grade, Braak stage, Thal phase, CERAD |

**Cohort:** 25 donors total — 13 AD, 8 Control, 4 Active Control.

**Disease definition:** Neuropathological ADNC grade (High = A3B3C3), not clinical diagnosis. AD group is clinically heterogeneous (typical LOAD, EOAD, LvPPA, bvFTD with AD neuropathology) — ADNC grade is the operationally consistent definition.

**Active Controls:** Cognitively impaired patients with non-AD primary diagnosis (PSP, CBD) and low/intermediate ADNC. Treated as a separate sensitivity cohort, not collapsed with Controls or included in primary modeling.

---

## Layer 0 — Pre-aggregation QC and metadata construction

### Purpose
Build a complete, harmonized donor metadata table and attach it to `adata_transcript_loose` before any modeling. Filter barcodes and transcripts to a quality floor.

### Step 0.1 — Metadata join and ID harmonization

**Inputs:** xlsx (SMC donors) + pptx (PO donors)

**Problem:** Sample ID formats are inconsistent across sources.
- xlsx uses `SMC-036`, pptx uses `P_005`, AnnData uses `SMC036` / `PO05`

**Harmonization map:**
- `SMC-XXX` → `SMCXXX` (strip hyphen)
- `P_00X` → `PO0X` (reformat prefix and zero-pad)

**Output:** Unified metadata table with one row per donor:

| Field | Source | Notes |
|---|---|---|
| `donor_id` | both files | harmonized to AnnData format |
| `condition` | both files | AD / Control / Active Control |
| `age` | both files | age at death (years) |
| `sex` | both files | M / F |
| `braak_stage` | both files | ordinal 0–VI |
| `thal_phase` | both files | amyloid phase |
| `cerad_score` | both files | neuritic plaque density |
| `apoe` | pptx (PO donors) | E3/E3, E3/E4, E4/E4 |
| `adnc_grade` | both files | Low / Intermediate / High |

**Flag any donor with missing values** before proceeding. Do not drop silently — record in a QC log.

### Step 0.2 — Barcode merge and obs enrichment

Match barcodes between `adata_sr` and `adata_transcript_loose`. Transfer the following fields from `adata_sr.obs` into `adata_transcript_loose.obs`:

- `cell_type` — high-resolution cell type label (Microglia, Astrocyte, Oligodendrocyte, Excitatory neuron, Inhibitory neuron, OPC, Vascular cell, Lymphocyte)
- `doublet_score` — Scrublet doublet probability
- `pct_counts_mt` — mitochondrial read fraction (RNA quality proxy; transfers here because `adata_transcript_loose.obs` does not contain it)

Then join the unified metadata table from Step 0.1 onto `adata_transcript_loose.obs` via `donor` column.

**Code pattern:**
```python
import scanpy as sc
import pandas as pd

# Transfer fields from adata_sr via barcode index
transfer_cols = ['cell_type', 'doublet_score', 'pct_counts_mt']
adata_tx.obs = adata_tx.obs.join(
    adata_sr.obs[transfer_cols],
    how='left'
)

# Join external metadata
adata_tx.obs = adata_tx.obs.join(
    metadata_df.set_index('donor_id'),
    on='donor',
    how='left'
)
```

Only barcodes with a 1:1 match in `adata_sr` are retained. Unmatched barcodes are dropped and logged.

### Step 0.3 — Barcode confidence filter

Exclude barcodes that fail either criterion:

| Filter | Threshold | Rationale |
|---|---|---|
| `doublet_score` | < 0.3 | Remove likely doublets before aggregation |
| `cell_type` transfer | must be non-null | Drop barcodes with no label match |

Log the percentage of barcodes dropped per cell type per donor. If any cell type loses >20% of barcodes for a given donor, flag that donor-cell type combination.

### Step 0.4 — Transcript prevalence filter

Before aggregation, remove transcripts that are too sparse to model reliably.

**Threshold:** A transcript must have count > 0 in at least 40% of donors within a given cell type's analysis.

This is applied per cell type, not globally — a transcript absent in most Vascular cells may still be well-detected in Excitatory neurons.

```python
def prevalence_filter(count_matrix, donors, min_prevalence=0.40):
    detected = (count_matrix > 0).sum(axis=0)
    n_donors = count_matrix.shape[0]
    return count_matrix.loc[:, detected / n_donors >= min_prevalence]
```

### Step 0.5 — Covariate completeness audit

Verify the following covariates are present and non-null for all donors before handing to DRIMSeq:

| Covariate | Expected | Status |
|---|---|---|
| `age` | all 25 donors | from xlsx + pptx |
| `sex` | all 25 donors | from xlsx + pptx |
| `braak_stage` | all 25 donors | from xlsx + pptx |
| `pct_counts_mt` | all 25 donors | computed from adata_sr transfer |
| RIN | — | **Confirmed absent** from all sources. `pct_counts_mt` used as RNA quality proxy. Document in methods. |

Any donor with missing covariates after join: flag in QC log, decide exclude vs. impute before proceeding.

### Step 0.6 — Pseudo-bulk aggregation

For each cell type separately, aggregate transcript counts per donor:

```python
def pseudobulk(adata, cell_type, min_cells=10):
    sub = adata[adata.obs['cell_type'] == cell_type]
    donors = sub.obs['donor'].unique()
    
    pb_counts = {}
    pb_meta = {}
    
    for donor in donors:
        cells = sub[sub.obs['donor'] == donor]
        if len(cells) < min_cells:
            continue  # skip — log as excluded
        pb_counts[donor] = cells.X.sum(axis=0)
        pb_meta[donor] = {
            'condition': cells.obs['condition'].iloc[0],
            'age': cells.obs['age'].iloc[0],
            'sex': cells.obs['sex'].iloc[0],
            'braak_stage': cells.obs['braak_stage'].iloc[0],
            'median_pct_mt': cells.obs['pct_counts_mt'].median(),
            'n_cells': len(cells)
        }
    
    return pd.DataFrame(pb_counts).T, pd.DataFrame(pb_meta).T
```

**Minimum cell threshold:** 10 cells per donor per cell type. Donor-cell type combinations below this are excluded from that cell type's analysis and logged. They are not excluded from other cell types.

**Output per cell type:**
- Count matrix: rows = donors, columns = transcripts (post-prevalence filter)
- Sample metadata table: rows = donors, columns = covariates

---

## Layer 1 — DRIMSeq Dirichlet-multinomial modeling

### Purpose
Fit a Dirichlet-multinomial generalized linear model to detect differential transcript usage (DTU) between AD and Control donors, controlling for age, sex, disease severity, and RNA quality.

### Step 1.1 — Active Control treatment

Active Controls (n=4) are **excluded from the primary DRIMSeq run**.

Primary comparison: AD (n=13) vs. Control (n=8) only.

Active Controls are retained for a separate sensitivity analysis run after primary results are obtained. Expected behavior: if an isoform switch is driven by AD neuropathology, Active Controls (low ADNC) should show PSI values intermediate between AD and Control. This is a post-hoc confirmation, not a primary test.

### Step 1.2 — dmFilter

Apply DRIMSeq's built-in filter before model fitting. Use conservative thresholds appropriate for sparse long-read data:

```r
library(DRIMSeq)

d <- dmDSdata(counts = count_matrix, samples = sample_metadata)

d <- dmFilter(
  d,
  min_samps_gene_expr   = 0.6 * ncol(count_matrix),  # 60% of donors
  min_samps_feature_expr = 0.4 * ncol(count_matrix), # 40% of donors
  min_gene_expr          = 10,
  min_feature_expr       = 5
)
```

These are more stringent than DRIMSeq defaults. The 75 transcripts/cell average in the transcript object makes default settings too permissive.

### Step 1.3 — Model formula

```r
design <- model.matrix(
  ~ condition + age + sex + median_pct_mt + braak_stage,
  data = samples(d)
)

d <- dmPrecision(d, design = design)
d <- dmFit(d, design = design)
d <- dmTest(d, coef = "conditionAD")
```

**Formula:** `~condition + age + sex + median_pct_mt + braak_stage`

| Covariate | Type | Role |
|---|---|---|
| `condition` | factor (AD / Control) | primary effect of interest |
| `age` | continuous | confound — AD donors typically older |
| `sex` | factor (M / F) | confound — sex-differential neuroinflammation |
| `median_pct_mt` | continuous | RNA quality proxy (RIN unavailable; see Layer 0) |
| `braak_stage` | ordinal (0–VI) | disease severity gradient; positions donors more precisely than binary condition label |

**Note on braak_stage vs. RIN:** RIN scores are absent from all available metadata sources (xlsx, pptx, AnnData obs, AnnData uns). `median_pct_mt` — median mitochondrial read fraction per donor per cell type, computed during pseudo-bulk aggregation from `adata_sr.obs.pct_counts_mt` — is used as the RNA quality surrogate. High mitochondrial fraction is a direct molecular indicator of cell damage and RNA degradation.

### Step 1.4 — stageR two-stage FDR control

Pass DRIMSeq results into stageR for two-stage false discovery rate control across 63,994 transcripts.

```r
library(stageR)

# Stage 1: gene-level p-values
pScreen <- results(d, level = "gene")$pvalue
names(pScreen) <- results(d, level = "gene")$gene_id

# Stage 2: transcript-level p-values
pConfirmation <- matrix(
  results(d, level = "feature")$pvalue,
  ncol = 1
)
rownames(pConfirmation) <- results(d, level = "feature")$feature_id

tx2gene <- results(d, level = "feature")[, c("feature_id", "gene_id")]

stageRObj <- stageRTx(
  pScreen       = pScreen,
  pConfirmation = pConfirmation,
  pScreenAdjusted = FALSE,
  tx2gene       = tx2gene
)
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)
```

**Stage 1:** Gene-level screening — does this gene have any DTU at all? Controls FWER at gene level.

**Stage 2:** Transcript-level confirmation — which specific isoforms shift? Only genes passing Stage 1 are tested at Stage 2. Dramatically reduces the multiple testing burden across 63,994 transcripts.

**Output:** Per-cell-type table of significant DTU events with gene ID, transcript ID, ΔPSI (AD mean PSI − Control mean PSI), and stageR-adjusted p-values.

---

## Layer 2 — Prioritization and cross-cell-type concordance

### Purpose
Filter stageR hits to a tractable high-confidence candidate list focused on AD drug target genes, and characterize each candidate's behavior across cell types.

### Step 2.1 — Drug target mandatory gate

Only genes present in the predefined AD drug target gene list advance to IsoformSwitchAnalyzeR. All other DTU hits are retained in supplementary tables but do not proceed to structural analysis.

The AD drug target gene list must be defined before analysis begins. Candidate genes include, but are not limited to: BACE1, PSEN1, PSEN2, APP, MAPT, APOE, ACHE, NMDA receptor subunits (GRIN2A, GRIN2B), GSK3B, CDK5.

This gate is mandatory. A gene failing it does not advance regardless of DTU significance or effect size.

### Step 2.2 — Scoring and ranking

For genes passing the mandatory gate, compute a priority score:

| Criterion | Weight |
|---|---|
| |ΔPSI| magnitude of dominant isoform | high |
| Number of cell types showing concordant DTU | high |
| Switching isoform is protein-coding (not NMD/retained intron) | mandatory sub-gate |
| Gene has known PDB structural data | high |
| ΔPSI direction consistent (not discordant across cell types) | moderate |

Top-scoring genes from each cell type analysis advance to Layer 3. Practical target: 5–10 candidate gene × isoform pairs total.

### Step 2.3 — Cross-cell-type concordance classification

For each drug target gene with ≥1 significant DTU hit, compute ΔPSI direction and magnitude across all 7 cell types and classify:

| Class | Definition | Interpretation |
|---|---|---|
| Pan-cellular switcher | Consistent direction in ≥5 cell types | Robust disease-wide effect |
| Cell-type restricted | Significant in ≤2 cell types | Cell-type-specific splicing regulation |
| Discordant | Opposite direction in different cell types | Most biologically complex; interpret with caution |

Discordant switchers are flagged but not excluded — they represent the most interesting cases for cell-type-specific drug target biology.

---

## Layer 3 — IsoformSwitchAnalyzeR functional annotation

### Purpose
Convert DTU hits into biological mechanism narratives: determine which candidate isoforms are protein-coding, what functional domains they gain or lose, and whether their sequence differences fall at the drug binding interface.

### Step 3.1 — TransDecoder ORF verification

Before IsoformSwitchAnalyzeR, verify ORF completeness for all candidate isoform FASTA sequences using TransDecoder with BLASTP homology evidence.

```bash
TransDecoder.LongOrfs -t candidate_isoforms.fasta
blastp -query candidate_isoforms.fasta.transdecoder_dir/longest_orfs.pep \
       -db uniprot_sprot.fasta -outfmt 6 -max_target_seqs 1 \
       -out blastp.outfmt6 -num_threads 8
TransDecoder.Predict -t candidate_isoforms.fasta \
                     --retain_blastp_hits blastp.outfmt6 --single_best_only
```

Feed TransDecoder's GFF3 output into ISA via `addORFfromGTF()` to override CPC2 predictions for non-canonical transcripts where reference annotation is unreliable.

Any isoform without a complete ORF (intact start codon, stop codon, no premature stop) is **excluded from structural modeling**. It is retained in functional tables but flagged.

### Step 3.2 — NMD sensitivity prediction

ISA predicts nonsense-mediated decay sensitivity for each candidate isoform. Isoforms predicted NMD-sensitive are **deprioritized** — they are unlikely to produce stable protein and building AlphaFold2 models for them is computationally wasteful and logically incoherent.

NMD-sensitive candidates are moved to a supplementary "potential regulatory isoforms" table. They may still be biologically interesting (NMD-based regulatory switching is a real mechanism) but they do not advance to structural analysis.

### Step 3.3 — Functional consequence annotation

Run the full consequence set in `analyzeSwitchConsequences()`:

```r
library(IsoformSwitchAnalyzeR)

switchAnalyzeRlist <- analyzeSwitchConsequences(
  switchAnalyzeRlist,
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'signal_peptide_identified',
    'tss',    # alternative transcription start site
    'tes'     # alternative transcription end site
  )
)
```

**Consequence types and their relevance:**

| Consequence | Relevance to drug binding |
|---|---|
| `domains_identified` | Pfam domain gain/loss — coarse but fast |
| `intron_retention` | Retained introns often create frameshifts near binding regions |
| `tss` / `tes` | Alternative N/C termini can alter membrane anchoring or signal peptides for receptor targets |
| `NMD_status` | Gate for structural arm advancement |
| `signal_peptide_identified` | Subcellular localization changes affect drug accessibility |

### Step 3.4 — Manual binding interface residue check

ISA's Pfam annotation is gene-level and coarse. For each candidate that passes NMD and ORF filters, perform a manual binding interface check:

1. Retrieve the canonical isoform crystal structure from PDB
2. Identify binding interface residues: any residue within 4Å of the co-crystallized ligand
3. Align canonical vs. disease-state isoform protein sequences using Biopython PairwiseAligner (BLOSUM62, structure-guided anchoring)
4. Flag any binding interface residue that is: altered (substitution), deleted, or repositioned by an insertion in the disease-state isoform

```python
from Bio import pairwise2
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")

alignments = pairwise2.align.globalds(
    canonical_seq,
    disease_seq,
    blosum62,
    -10, -0.5  # gap open, extend
)
```

Only candidates with ≥1 binding interface residue affected advance to AlphaFold2 modeling. This filter converts the candidate list from "has functional consequences" to "has consequences specifically at the pharmacologically relevant site."

---

## Layer 4 — Structural arm

### Purpose
Model the 3D structure of each disease-state isoform, compare its binding pocket geometry to the canonical isoform, and estimate the effect on drug binding affinity through ensemble docking and MM-GBSA rescoring.

### Step 4.1 — AlphaFold2 prediction

Run AlphaFold2 monomer prediction for each candidate disease-state isoform FASTA sequence. Use full MSA (not reduced).

```bash
python run_alphafold.py \
  --fasta_paths=disease_isoform.fasta \
  --model_preset=monomer \
  --db_preset=full_dbs \
  --output_dir=af2_output/
```

**pLDDT gate (mandatory):** After prediction, extract per-residue pLDDT scores for the alternatively spliced region specifically. If the spliced region has median pLDDT < 70, the model is flagged as low-confidence at the region of interest.

- pLDDT ≥ 70 in spliced region → advance to docking
- pLDDT < 70 in spliced region → flag; report pLDDT score; do not report docking ΔG as a primary result

Report pLDDT scores for all binding site residues explicitly. Reviewers with structural biology backgrounds will check this.

### Step 4.2 — Pocket comparison (Fpocket)

Run Fpocket on both canonical (PDB structure) and disease-state (AF2 model) structures:

```bash
fpocket -f canonical.pdb
fpocket -f disease_isoform.pdb
```

**Report per isoform:**
- Pocket volume (Å³)
- Druggability score
- Key residue composition at pocket wall

Differences in pocket volume >20% or druggability score >0.2 are flagged as structurally significant.

### Step 4.3 — Structural alignment (TM-align)

Compute global and local structural similarity between canonical and disease-state isoforms:

```bash
TMalign canonical.pdb disease_isoform.pdb -o tmalign_output
```

**Report:**
- Global TM-score (>0.5 = same fold; <0.5 = different fold)
- RMSD restricted to binding site residues only (more sensitive than global RMSD for drug binding assessment)

### Step 4.4 — Ensemble docking

Dock each existing AD drug against both the canonical (PDB) and disease-state (AF2) structures using AutoDock Vina or Gnina.

```bash
# Prepare receptor
prepare_receptor4.py -r canonical.pdb -o canonical.pdbqt
prepare_receptor4.py -r disease_isoform.pdb -o disease_isoform.pdbqt

# Prepare ligand
prepare_ligand4.py -l drug.mol2 -o drug.pdbqt

# Dock — same box coordinates for both structures
vina --receptor canonical.pdbqt --ligand drug.pdbqt \
     --center_x X --center_y Y --center_z Z \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 32 --num_modes 10 --out canonical_poses.pdbqt

vina --receptor disease_isoform.pdbqt --ligand drug.pdbqt \
     --center_x X --center_y Y --center_z Z \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 32 --num_modes 10 --out disease_poses.pdbqt
```

**Critical:** Use identical box coordinates (center and size) for both structures, defined from the canonical PDB co-crystallized ligand position. This ensures a like-for-like comparison.

**Run 10 independent docking runs per drug per structure.** Report the best pose from 10 runs — AutoDock Vina uses stochastic search and a single run may miss the global minimum.

### Step 4.5 — MM-GBSA rescoring

Rescore the top docking pose from each run using MM-GBSA implicit solvent to obtain a more reliable ΔG estimate than raw docking scores alone.

Use Gnina's built-in CNN rescoring or OpenMM with an implicit solvent force field (GBSA).

**Report per drug × isoform pair:**

| Metric | Description |
|---|---|
| ΔG canonical | Binding free energy estimate, canonical isoform |
| ΔG disease | Binding free energy estimate, disease-state isoform |
| ΔΔG | ΔG disease − ΔG canonical (positive = weaker binding in disease) |
| Pose RMSD | RMSD of disease-state top pose vs. canonical binding mode |
| Pharmacophore contacts | H-bond donors/acceptors retained, hydrophobic contacts retained |

**H3 threshold:** ΔΔG > 1.5 kcal/mol is used as the primary significance criterion, with the caveat that docking RMSE is typically 1.5–2.0 kcal/mol. Supplement with pose RMSD > 2Å and pharmacophore contact loss as orthogonal evidence. Do not rely on ΔΔG alone.

---

## Layer 5 — Validation

### Purpose
Provide at least one orthogonal line of evidence that the computationally identified disease-state isoforms exist as translated protein in AD brain tissue.

### Option A — AMP-AD proteomics rescue (no wet lab required)

Query the AMP-AD Knowledge Portal (Synapse) for existing brain proteomics datasets. Search for isoform-specific junction peptides — tryptic peptides that span a splice junction unique to the disease-state isoform.

**Steps:**
1. Build a custom isoform-aware FASTA database from the transcript catalog
2. Search raw MS data using MSFragger or Philosopher with the custom FASTA
3. Report any detected junction peptides with spectral evidence

Even one junction peptide detected in AD brain MS data converts the claim from "transcriptomically predicted" to "proteomically supported."

### Option B — RT-PCR isoform confirmation (one day, if brain cDNA is available)

Design primers flanking the alternatively spliced exon in the top candidate gene, with one primer in the exon-specific sequence.

Expected result: canonical isoform produces band of size X, disease-state isoform produces band of size Y. Report band ratio AD vs. Control.

**Methods statement for either option:**

> "The dominant disease-state isoform identified computationally was validated at the [protein / RNA] level using [AMP-AD proteomics / RT-PCR] on independent AD post-mortem brain samples."

---

## Parallel sensitivity check — Allos

Allos (McAndrew et al., 2026, bioRxiv doi:10.64898/2026.03.24.713944) is a Python-native isoform-level single-cell transcriptomics toolkit integrated with the scverse ecosystem.

**Role in Viscacha:** Supplementary sensitivity check only. Run Allos on the top 5 DTU hits identified by DRIMSeq + stageR to assess concordance between pipelines. Concordance between independent frameworks strengthens confidence in primary findings.

**Do not use as primary pipeline.** Allos is a 2026 preprint and has not been peer-reviewed. Its results should not be reported as primary evidence in a graduation thesis.

---

## Complete pipeline summary

```
LAYER 0 — Pre-aggregation QC
  Step 0.1  ID harmonization + metadata join (xlsx + pptx → unified table)
  Step 0.2  Barcode merge: cell_type, doublet_score, pct_counts_mt from adata_sr
            + donor metadata from unified table → adata_transcript_loose.obs
  Step 0.3  Barcode confidence filter (doublet_score < 0.3, 1:1 label match)
  Step 0.4  Transcript prevalence filter (≥40% donors per cell type)
  Step 0.5  Covariate completeness audit (flag missing; confirm RIN absent)
  Step 0.6  Pseudo-bulk aggregation (sum per donor × cell type, min 10 cells)
            → compute median_pct_mt per donor × cell type

LAYER 1 — DRIMSeq statistical modeling
  Step 1.1  Active Control exclusion from primary run (→ sensitivity cohort)
  Step 1.2  dmFilter (conservative thresholds for sparse LR data)
  Step 1.3  DRIMSeq GLM: ~condition + age + sex + median_pct_mt + braak_stage
  Step 1.4  stageR two-stage FDR (gene-level screen → transcript confirmation)
            → significant DTU table per cell type

LAYER 2 — Prioritization
  Step 2.1  Drug target mandatory gate (predefined gene list)
  Step 2.2  Priority scoring (|ΔPSI|, concordance, coding status, PDB availability)
  Step 2.3  Cross-cell-type concordance classification
            → 5–10 candidate isoform pairs

LAYER 3 — IsoformSwitchAnalyzeR
  Step 3.1  TransDecoder ORF verification (BLASTP-supported)
  Step 3.2  NMD sensitivity gate (NMD-sensitive → deprioritized)
  Step 3.3  analyzeSwitchConsequences() (full consequence set)
  Step 3.4  Manual 4Å binding interface residue check (Biopython + PDB)
            → candidates with binding interface alterations only

LAYER 4 — Structural arm
  Step 4.1  AlphaFold2 monomer prediction + pLDDT gate (≥70 at spliced region)
  Step 4.2  Fpocket pocket comparison (volume, druggability, residue composition)
  Step 4.3  TM-align structural comparison (global TM-score + binding site RMSD)
  Step 4.4  Ensemble docking (10 runs × drug × isoform, identical box coordinates)
  Step 4.5  MM-GBSA rescoring → ΔΔG, pose RMSD, pharmacophore contact retention

LAYER 5 — Validation
  Option A  AMP-AD proteomics: junction peptide search (MSFragger + custom FASTA)
  Option B  RT-PCR: isoform-discriminating primers on AD brain cDNA

PARALLEL   Allos sensitivity check on top 5 hits (supplementary only)
```

---

## Key methodological decisions and rationale

| Decision | Rationale |
|---|---|
| Pseudo-bulk rather than single-cell DTU | 75 transcripts/cell average makes single-cell zero-inflation intractable; pseudo-bulk bypasses this while preserving biological replication |
| DRIMSeq over satuRn | Dirichlet-multinomial model directly accounts for multinomial nature of isoform proportions; appropriate for this data structure |
| stageR two-stage FDR | 63,994 transcripts makes transcript-level testing without gene-level pre-screening statistically indefensible |
| Braak stage as disease severity covariate | More precise than binary condition label; positions Active Controls correctly on the severity axis |
| pct_counts_mt as RIN proxy | RIN confirmed absent from all metadata sources; mitochondrial fraction is a direct molecular RNA quality signal derivable from existing obs columns |
| Active Controls as sensitivity cohort | n=4 is underpowered for a three-level DRIMSeq factor; neuropathological staging (low ADNC) confirms they are intermediate, not equivalent to either primary group |
| TransDecoder over CPC2 alone | More robust ORF prediction for non-canonical transcripts using homology evidence; prevents structural modeling of non-translated isoforms |
| pLDDT gate at spliced region | AF2 produces a structure regardless of confidence; low pLDDT in the alternatively spliced region makes docking results unreliable and reviewer-challengeable |
| Ensemble docking (10 runs) | AutoDock Vina is stochastic; single-run results are not reproducible; 10 runs report the global minimum reliably |
| MM-GBSA rescoring | Raw docking scores have RMSE ~1.5–2.0 kcal/mol; rescoring with implicit solvent improves estimate quality for ΔΔG comparison |
| Condition defined by ADNC grade | Clinical diagnoses are heterogeneous (EOAD, LvPPA, bvFTD with AD neuropathology); neuropathological ADNC grade is the consistent operationally defined variable |

---

## Software and tools

| Tool | Version | Purpose |
|---|---|---|
| Python / scanpy | ≥1.9 | AnnData manipulation, barcode merge, metadata join |
| pandas | ≥1.5 | Metadata harmonization and tabular operations |
| Biopython | ≥1.80 | Sequence alignment, binding interface residue check |
| R / DRIMSeq | Bioconductor | Dirichlet-multinomial DTU modeling |
| R / stageR | Bioconductor | Two-stage FDR control |
| R / IsoformSwitchAnalyzeR | Bioconductor ≥2.0 | Functional consequence annotation |
| TransDecoder | ≥5.5 | ORF prediction with BLASTP homology support |
| AlphaFold2 | v2.3 | Disease-state isoform structure prediction |
| Fpocket | ≥4.0 | Binding pocket geometry comparison |
| TM-align | — | Global and local structural similarity |
| AutoDock Vina / Gnina | ≥1.2 / latest | Ensemble molecular docking |
| OpenMM | ≥7.7 | MM-GBSA rescoring (implicit solvent) |
| MSFragger / Philosopher | latest | Proteomics junction peptide search (Layer 5A) |
| Allos | 2026 preprint | Sensitivity check (supplementary only) |

---

*Pipeline name: Viscacha-v01*
*Cohort: Samsung Medical Center post-mortem brain, 25 donors (13 AD, 8 Control, 4 Active Control)*
*Data: Single-cell long-read RNA-seq, AnnData h5ad format*

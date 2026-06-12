# Viscacha Pipeline — Layer 1: Differential Transcript Usage
### Detailed Explanation of Every Step, Equation, and Decision

---

## Table of Contents

1. [Biological Question](#1-biological-question)
2. [Input Data](#2-input-data)
3. [Step 11 — Count-Based Transcript Filter](#3-step-11--count-based-transcript-filter)
4. [Step 12 — Fitting the DTU Model (satuRn)](#4-step-12--fitting-the-dtu-model-saturn)
5. [Step 13 — Two-Stage FDR Control (stageR) and ΔPSI](#5-step-13--two-stage-fdr-control-stager-and-δpsi)
6. [Step 14 — Active Control Descriptive PSI](#6-step-14--active-control-descriptive-psi)
7. [The Braak Collinearity Problem and Option C](#7-the-braak-collinearity-problem-and-option-c)
8. [Cell Type Loop and Output](#8-cell-type-loop-and-output)
9. [Configuration Decisions Reference](#9-configuration-decisions-reference)
10. [Output Column Glossary](#10-output-column-glossary)

---

## 1. Biological Question

The pipeline tests whether **Alzheimer's Disease (AD) drug target genes** exhibit **differential transcript usage (DTU)** in post-mortem brain tissue. DTU means the *proportion* of transcript isoforms from a gene changes between conditions — not just whether the gene is more or less expressed overall (which is differential gene expression, DGE). A gene could be differentially expressed without DTU (all isoforms go up equally), or it could show DTU without any change in total expression (one isoform replaces another).

This distinction matters because many drugs target specific protein domains. If a gene switches from producing isoform A to isoform B in AD brain tissue, the drug's binding site may be present or absent depending on which isoform dominates. Layer 1 asks: which transcripts in AD drug target genes change their relative abundance in AD-affected brain?

---

## 2. Input Data

Inputs come from Layer 0, which aggregated single-cell RNA-seq read counts into **pseudobulk** matrices — one count per transcript per donor, within each cell type. A pseudobulk aggregation is used because single-cell models that treat each cell as an independent replicate inflate statistical power and produce false positives; pseudobulk collapses cells from the same donor into a single count that respects the true unit of biological replication (the donor).

**Files per cell type (from Layer 0):**
- `counts_<cell_type>.csv` — rows = donors, columns = transcripts; values = total read counts across all cells of that type in that donor
- `metadata_<cell_type>.csv` — rows = donors; columns include `condition`, `age`, `sex`, `braak_stage`, `median_pct_mt`
- `counts_<cell_type>_sensitivity.csv` / `metadata_<cell_type>_sensitivity.csv` — same format but includes **Active Control** donors (used in Step 14 only)

The matrix is transposed after loading (`counts_mat = t(as.matrix(counts_raw))`) so that it becomes **transcripts × donors** for all downstream steps.

**Cell types analyzed:**
Excitatory neuron, Inhibitory neuron, Oligodendrocyte, OPC, Astrocyte, Microglia, Vascular cell.

**Lymphocyte is excluded** because only 6 transcripts pass the prevalence filter (Step 11). Running a DTU model on 6 transcripts would be meaningless and would produce unreliable FDR estimates.

---

## 3. Step 11 — Count-Based Transcript Filter

Before fitting any model, transcripts that are too sparse or unreliably quantified are removed. This is done in three sequential stages. The order matters: Stage 1 applies first, and only transcripts surviving Stage 1 enter Stage 2, and so on.

### Why filter at all?

Long-read RNA-seq pseudobulk data is sparse. Many transcripts will have reads in only a handful of donors or in very low numbers. Including them in the model would:
- Produce near-zero count estimates with extreme variance, destabilizing the GLM
- Inflate the number of tests, increasing the multiple-testing burden
- Generate unreliable p-values for transcripts whose expression is mostly noise

### Minimum samples threshold

Before any stage, compute the minimum number of samples (donors) that must meet each threshold:

```
min_samps = max(2, floor(n_samples × MIN_SAMPS_FRAC))
```

With `MIN_SAMPS_FRAC = 0.40` and, e.g., 25 donors: `min_samps = max(2, floor(25 × 0.40)) = max(2, 10) = 10`.

**Why 40%?** A transcript that is expressed in fewer than 40% of donors within a group is too sparse to estimate a reliable mean or to compare between conditions. The `max(2, ...)` floor ensures that even for very small datasets (< 5 donors), at least 2 must express the transcript.

The threshold is computed once and applies identically to all three filter stages.

### Stage 1: Transcript count threshold

```
Keep transcript t  if:  |{j : count_tj ≥ MIN_TX_COUNT}| ≥ min_samps
```

Where `count_tj` is the pseudobulk count for transcript `t` in donor `j`, and `MIN_TX_COUNT = 5`.

In words: transcript `t` passes if it has **at least 5 counts** in **at least 40% of donors**. A count below 5 is too close to sequencing noise to confidently call "expressed."

### Stage 2: Gene total threshold

After Stage 1, compute gene-level totals per donor by summing all passing transcripts of each gene:

```
gene_total_gj = Σ_{t ∈ gene g} count_tj
```

Then keep transcript `t` (in gene `g`) only if:

```
|{j : gene_total_gj ≥ MIN_GENE_COUNT}| ≥ min_samps
```

With `MIN_GENE_COUNT = 10`. A transcript can individually pass Stage 1 (≥5 counts) but belong to a gene that is poorly covered overall. If the gene total is below 10, computing PSI = count_t / gene_total would be dividing by a very small number, making the proportion estimate unreliable (e.g., 5/7 = 71% vs 5/10 = 50% — a 3-count difference in the denominator changes the proportion by 21 percentage points). Stage 2 ensures the denominator of the PSI ratio is meaningful.

### Stage 3: Multi-transcript gene filter

After Stages 1 and 2, keep only transcripts belonging to genes with **at least 2 surviving transcripts**:

```
Keep transcript t  if:  |{t' ∈ gene g : t' survived stages 1 & 2}| ≥ 2
```

**This is the most conceptually important filter.** DTU is defined as a change in the *relative proportion* of isoforms. If only one isoform of a gene survives the count filters, there is nothing to compare proportions between — the single surviving transcript will always have PSI ≈ 1.0 in all samples. Testing this would produce p-values of 1.0 and waste statistical power. A gene must have at least 2 measurable isoforms to be a candidate for DTU.

### Filter breakdown tracking

The pipeline records how many transcripts fail at each stage (`filter_counts_breakdown()`). This is used to plot the step 11 retention figure, showing that most exclusions come from Stage 1 (low individual counts) and Stage 3 (single-isoform genes), not Stage 2.

---

## 4. Step 12 — Fitting the DTU Model (satuRn)

### What satuRn does

satuRn (Differential Transcript Usage in single-cell RNA-seq) fits a **quasi-binomial generalized linear model (GLM)** for each transcript separately. The quasi-binomial is appropriate because:

1. PSI (proportion of splicing) is naturally bounded in [0, 1], so a binomial family with logit link is mathematically correct.
2. RNA-seq counts are over-dispersed relative to a strict binomial — the variance exceeds the binomial expectation. The quasi-binomial family adds a dispersion parameter φ to model this excess variance:

```
Var(Y_ij) = φ × μ_ij × (1 - μ_ij) / n_ij
```

where `n_ij` is the total gene count for transcript `i` in sample `j`, and `μ_ij` is the expected PSI.

### The model

For each transcript `i`, the model is:

```
logit(E[PSI_ij]) = β₀ + β₁ × conditionAD_j + β₂ × age_j + β₃ × sex_j + β₄ × median_pct_mt_j
```

Where `logit(p) = log(p / (1−p))`. The model is estimated in the count space: satuRn passes `(count_ij, gene_total_ij − count_ij)` as the (success, failure) pair to the GLM, which is mathematically equivalent to modeling PSI.

The data entering the model is:
- **successes** = reads mapping to transcript `i` in donor `j`
- **trials** = total reads mapping to any transcript of the same gene in donor `j`

This naturally encodes the constraint that all transcript proportions of a gene must sum to 1.

### Covariates and why each is included

| Covariate | Rationale |
|---|---|
| `conditionAD` | The effect of interest. Encodes AD = 1, Control = 0 with Control as the reference level. |
| `age` | Post-mortem brain tissue from elderly donors. Age independently affects splicing. Omitting it would confound the AD effect if AD donors are systematically older. |
| `sex` | Sex-linked differences in both splicing and AD risk are well-documented. Omitting it would introduce unexplained variance. |
| `median_pct_mt` | The median fraction of mitochondrial reads per cell, aggregated to donor level. High mitochondrial transcript fraction is a marker of cell stress or low-quality nuclei. It captures technical variation in library quality across donors. Including it removes this batch-like noise from the DTU estimate. |

`condition` is encoded as a factor with levels `c("Control", "AD")`, making **Control the reference level**. This means the `conditionAD` coefficient represents the change *in AD relative to Control*.

### The design matrix and contrast

The design matrix `design` is built with `model.matrix(formula, data = meta)`. Its columns are the model coefficients: `(Intercept)`, `conditionAD`, `age`, `sexMale` (or similar), `median_pct_mt`.

To test for DTU, the pipeline constructs a **contrast matrix** `L`:

```
L = n_coefficients × 1 matrix of zeros
L["conditionAD", "AD_vs_Control"] = 1
```

This selects only the `conditionAD` coefficient for testing. The test asks: *is β_conditionAD significantly different from zero?* A non-zero β_conditionAD means the logit-transformed PSI differs between AD and Control after adjusting for age, sex, and mitochondrial fraction.

### satuRn's empirical FDR

When a cell type has ≥ 500 transcripts passing the filter, satuRn computes **empirical p-values** by permuting sample labels and re-testing. The resulting permuted t-statistics define the null distribution, and the observed t-statistic is compared against it. This corrects for the fact that the quasi-binomial t-statistic may not follow a theoretical t-distribution exactly when the model is misspecified or the data is highly structured.

For smaller cell types (< 500 transcripts), satuRn falls back on theoretical p-values from the quasi-binomial t-test.

The pipeline uses whichever type was computed (`empirical_pval` if it exists, otherwise `pval`) in the volcano plot, but **stageR always uses the raw `pval`** (not empirical, not FDR-adjusted) because stageR performs its own FDR adjustment.

### SummarizedExperiment requirements

satuRn requires the data to be in a `SummarizedExperiment` object with specific structure:
- `assays$counts` — the transcript × donor count matrix
- `colData` — donor metadata (with `condition` and `sex` as factors)
- `rowData` — must contain columns named exactly `isoform_id` and `gene_id` (satuRn parses these internally to group transcripts by gene)

The `gene_id` is extracted from transcript names by removing the last hyphen-suffixed number:

```
"CAMK2B-204"  →  gene_id = "CAMK2B",  isoform_id = "CAMK2B-204"
```

This convention assumes transcript IDs follow the `GENE-NNN` format from GENCODE/Ensembl.

---

## 5. Step 13 — Two-Stage FDR Control (stageR) and ΔPSI

### Why standard FDR correction is problematic for DTU

A standard Benjamini-Hochberg (BH) correction applied to all transcript p-values treats each transcript as independent. But transcripts within the same gene are correlated — they share the same donors, the same total count denominator, and are often co-regulated. Testing 5 isoforms of gene X and 5 isoforms of gene Y at α = 0.05 and then applying BH globally treats all 10 tests as exchangeable, which over-penalizes genes with many isoforms (more tests = harder to pass BH).

stageR implements a **hierarchical two-stage testing procedure** that respects the gene → transcript structure.

### Stage 1: Gene-level screening (pScreen)

For each gene `g`, compute a gene-level p-value:

```
pScreen_g = min_{t ∈ gene g} p_t
```

The gene-level p-value is the **minimum transcript p-value** across all transcripts in that gene. This is a form of the Simes test: if any transcript of a gene shows a significant change, the gene is flagged as potentially having DTU.

Apply **Benjamini-Hochberg correction** to `pScreen` across all genes. A gene passes screening if `pScreen_g_adj < ALPHA = 0.05`. Only genes passing screening advance to Stage 2.

**Why minimum p-value?** DTU requires at least one transcript to change proportion. Taking the minimum captures this "any significant change" criterion at the gene level.

### Stage 2: Transcript-level confirmation (pConfirmation)

For each transcript `t` within a gene that passed Stage 1, the transcript-level p-value is tested using the **DTU method** (`stageWiseAdjustment(method = "dtu")`).

The DTU method adjusts the per-transcript significance threshold `α_tx` such that the overall false discovery rate (OFDR) is controlled at ALPHA across the entire two-stage procedure. Specifically, the method from Van den Berge et al. (2017) ensures:

```
OFDR = P(at least one false positive at either stage) ≤ ALPHA
```

The adjusted transcript-level p-value `padj_tx` is reported. A transcript is declared **significant** if:

```
padj_gene < ALPHA  AND  padj_tx < ALPHA
```

Both conditions must hold: the gene must first pass screening, *then* the individual transcript must pass confirmation. This two-gate requirement maintains strong control of the OFDR.

**Advantage over global BH:** Genes with many isoforms are not penalized more than genes with few isoforms, because the gene-level screening step pools all isoform evidence together, and the per-transcript confirmation step only needs to compete within the confirmed gene's budget, not globally.

### ΔPSI: the effect size

The stageR p-values tell you *whether* a change is real. ΔPSI tells you *how large* the change is and in which direction. It is computed separately from the model, directly from raw counts:

**Step 1: Compute PSI per transcript per donor**

```
PSI_tj = count_tj / gene_total_gj     if gene_total_gj > 0
PSI_tj = NA                            if gene_total_gj = 0
```

**Step 2: Average within each condition**

```
PSI_AD_t    = mean( PSI_tj  for all AD donors j,  ignoring NA )
PSI_ctrl_t  = mean( PSI_tj  for all Control donors j,  ignoring NA )
```

**Step 3: Compute the difference**

```
ΔPSI_t = PSI_AD_t − PSI_ctrl_t
```

A positive ΔPSI means transcript `t` makes up a larger fraction of the gene's output in AD compared to Control. A negative ΔPSI means it makes up a smaller fraction.

**Why not use the model's β coefficient as the effect size?** The model coefficient is on the logit scale (log-odds), which is not directly interpretable. ΔPSI is a linear difference in proportion, bounded in [−1, +1], and immediately meaningful: "this isoform goes from 20% to 66% of the gene's reads in AD" is a ΔPSI of +0.46.

**Why are PSIs computed on filtered counts only?** The PSI computation in Step 13 uses the same filtered count matrix that entered the model (from the SummarizedExperiment after filtering). This ensures that the denominator (gene total) is computed over the same set of transcripts that were modeled. Including transcripts that failed the filter in the denominator would deflate PSI values inconsistently.

---

## 6. Step 14 — Active Control Descriptive PSI

### What is an Active Control?

Active Controls are **non-demented individuals** whose brains are nevertheless pathologically active or show early-stage changes. They serve as a comparison group that is not simply "normal aging." Including Active Controls in the visualization allows us to ask: does the isoform switch we observe in AD also appear in cognitively intact individuals with some pathology, or is it specific to clinical AD?

**Crucially, no statistical test is run on Active Controls.** Step 14 is purely descriptive. The rationale is:
1. The Active Control sample size is typically too small for reliable DTU testing.
2. We already have our significance calls from the AD vs. Control comparison. The Active Control PSI is an *interpretive aid*, not a separate test.
3. Running a second statistical test on Active Controls would require another round of multiple testing correction, which would reduce power and complicate interpretation.

### Computation

For each transcript `t` in the sensitivity dataset:

```
PSI_active_tj = count_tj / gene_total_gj    for all Active Control donors j
psi_active_ctrl_t = mean( PSI_active_tj,  ignoring NA )
```

The sensitivity count matrix contains **Active Control donors only** (the primary count matrix contains AD + Control donors).

### The "% of switch" metric (Step 14 plot)

In the dumbbell plot visualization, a "% of switch" metric is computed:

```
% of switch = (psi_active_ctrl − psi_ctrl) / (psi_AD − psi_ctrl) × 100
```

This asks: *of the total isoform proportion change between Control and AD, how much of it is already present in Active Controls?*

- **0%** → Active Controls look exactly like Controls; the switch is specific to AD.
- **100%** → Active Controls look exactly like AD cases; the switch may be driven by aging or early pathology, not AD per se.
- **50%** → Active Controls are halfway between; the switch is partially explained by shared pathology.
- **Negative** → Active Controls go in the opposite direction; the switch is truly AD-specific and may even be counter-regulated in non-demented aging.
- **>100%** → Active Controls overshoot AD; the isoform is at an extreme in non-demented pathology.

This metric is reported in the step 14 visualization (dumbbell plot) to flag whether hits from the main analysis are AD-specific or potentially confounded by aging or shared pathology.

---

## 7. The Braak Collinearity Problem and Option C

### What is Braak stage?

Braak stage (0–6) is a histological score of neurofibrillary tangle (tau pathology) burden in the brain. It is the most widely used post-mortem staging system for AD severity. Higher Braak = more tau pathology.

### The collinearity problem

In this dataset:

```
r(braak_stage, condition) ≈ 0.914
```

AD donors almost universally have high Braak scores (5–6), and Control donors almost universally have low Braak scores (0–2). This near-perfect correlation between `braak_stage` and `condition` creates **multicollinearity** if both are included in the same regression model.

**What multicollinearity does to the GLM:**
- The design matrix becomes nearly rank-deficient — its columns are nearly linearly dependent
- The coefficient estimates for `conditionAD` and `braak_stage` become unstable; small changes in the data produce large swings in each coefficient
- Standard errors of both coefficients inflate dramatically (variance inflation factor, VIF >> 10)
- The test for `conditionAD` loses power; p-values become unreliable

**Three options were considered:**

| Option | Approach | Problem |
|---|---|---|
| A | Drop braak_stage entirely | Cannot separate AD diagnosis from tau burden; braak is a potential confounder |
| B | Include braak_stage in all models | VIF inflates SEs, reduces power, unstable estimates |
| C | Two separate models | Primary model without braak; sensitivity model with braak |

**Option C was chosen** as the principled compromise:

1. **Primary model** (`FORMULA_PRIMARY = ~ condition + age + sex + median_pct_mt`): run without braak_stage. Full statistical power. Results represent the AD vs. Control difference after adjusting for age, sex, and cell quality.

2. **Braak sensitivity model** (`FORMULA_BRAAK = ~ condition + age + sex + median_pct_mt + braak_stage`): run on the same data. Note that because braak_stage has many NA values in some donors, samples with missing braak data are dropped for this model only.

3. **Robustness flag**: a transcript is marked `robust_to_braak = TRUE` if it is significant in *both* models. If it is significant in the primary model but not in the braak model, it means the signal may be driven by tau burden (braak) rather than AD diagnosis per se. This does not mean it is a false positive — it means its interpretation is more nuanced.

**Why not use a mixed model or include braak as a random effect?** The sample sizes per cell type (typically 15–30 donors) are too small for stable random effects estimation. Option C avoids this instability while providing honest transparency about the braak confound.

---

## 8. Cell Type Loop and Output

The pipeline runs independently for each cell type. This is correct because:
- Cells from different types have vastly different transcriptomic profiles
- Pooling across cell types would mix incompatible count distributions
- A transcript may be DTU in excitatory neurons but not in oligodendrocytes — these are separate biological questions

For each cell type, the pipeline:
1. Loads counts + metadata
2. Runs Step 11 filter
3. Runs primary model (Step 12 + 13)
4. Runs braak sensitivity model (Step 12 + 13 on braak-complete donors)
5. Marks `robust_to_braak`
6. Loads sensitivity data and adds `psi_active_ctrl` (Step 14)
7. Writes `dtu_results_<cell_type>.csv`

After all cell types, a combined `dtu_significant_all_celltypes.csv` is written containing only transcripts that pass both `padj_gene < 0.05` AND `padj_tx < 0.05`.

### QC log

For each cell type and model, the pipeline records:
- Number of transcripts before and after filtering
- Number of genes tested
- Number of significant genes and transcripts

This is written to `outputs/layer1/qc/layer1_qc.txt` and allows rapid diagnosis of whether a cell type had insufficient data or unexpected results.

---

## 9. Configuration Decisions Reference

| Parameter | Value | Rationale |
|---|---|---|
| `MIN_TX_COUNT` | 5 | Minimum counts to call a transcript "expressed" per donor. Below 5 is indistinguishable from sequencing noise in long-read data. |
| `MIN_GENE_COUNT` | 10 | Minimum gene-level total needed to compute a reliable PSI denominator. Ensures PSI ratios are not driven by near-zero denominators. |
| `MIN_SAMPS_FRAC` | 0.40 | Fraction of donors that must express the transcript. 40% ensures we test transcripts present across the group, not rare individual expression events. |
| `ALPHA` | 0.05 | Standard FDR threshold. Applied at both stageR stages. Controls OFDR at 5%. |
| `N_CORES` | 8 | Parallel cores for `fitDTU`. Reduces wall time; satuRn is transcript-parallelizable via BiocParallel. |
| `COND_LEVELS` | `c("Control", "AD")` | Ensures Control is the reference level in all factor encodings. The `conditionAD` coefficient then represents AD relative to Control (positive = higher in AD). |

---

## 10. Output Column Glossary

Each `dtu_results_<cell_type>.csv` file contains:

| Column | Source | Meaning |
|---|---|---|
| `transcript_id` | rownames of count matrix | Full transcript identifier (e.g. `CAMK2B-204`) |
| `gene_id` | Parsed from transcript_id | Parent gene (e.g. `CAMK2B`) |
| `padj_gene` | stageR Stage 1 | BH-adjusted gene-level p-value (minimum transcript pval per gene, BH-corrected across all genes) |
| `padj_tx` | stageR Stage 2 | stageR-adjusted transcript-level p-value (within-gene confirmation, OFDR-controlled) |
| `psi_AD` | Step 13 | Mean PSI across AD donors: `mean(count_t / gene_total, na.rm=TRUE)` |
| `psi_ctrl` | Step 13 | Mean PSI across Control donors |
| `delta_psi` | Step 13 | `psi_AD − psi_ctrl`. Positive = more abundant in AD. |
| `pval` | satuRn | Raw quasi-binomial t-test p-value for `conditionAD` coefficient |
| `empirical_pval` | satuRn | Permutation-based p-value (only for cell types with ≥500 transcripts) |
| `regular_FDR` | satuRn | BH-adjusted pval across transcripts within cell type (NOT used for significance calling; stageR is used instead) |
| `empirical_FDR` | satuRn | BH-adjusted empirical_pval (same caveat) |
| `psi_active_ctrl` | Step 14 | Mean PSI across Active Control donors (descriptive, no test) |
| `robust_to_braak` | Option C | `TRUE` if significant in both primary and braak sensitivity models. `NA` if braak model failed to run. |
| `cell_type` | Loop variable | Which cell type this row belongs to |

### What counts as a significant DTU hit?

```
padj_gene < 0.05  AND  padj_tx < 0.05
```

Both conditions must hold. `padj_gene` passes the gene-level screening stage; `padj_tx` passes the transcript-level confirmation stage. A transcript with `padj_gene = 0.03` but `padj_tx = 0.20` is a gene with some DTU signal but this specific transcript is not individually confirmed. A transcript with `padj_gene = 0.20` cannot have `padj_tx` evaluated at all (stageR returns NA for transcripts in genes that did not pass screening).

### What `robust_to_braak = TRUE` means for interpretation

A hit with `robust_to_braak = TRUE` means the isoform switch between AD and Control persists even after including Braak stage in the model. Given the near-perfect correlation between Braak and AD diagnosis, passing both models provides stronger evidence that the DTU signal is associated with the AD clinical diagnosis independently of tau burden severity. Hits with `robust_to_braak = FALSE` remain valid AD-associated findings from the primary analysis but should be interpreted with the caveat that tau burden may be a more parsimonious explanation than the clinical diagnosis.

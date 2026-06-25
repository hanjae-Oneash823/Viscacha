# ============================================================
# Step 12: Build SummarizedExperiment and fit/test DTU
# ============================================================
# Uses satuRn 1.18.0
# Returns: tested SummarizedExperiment + data.frame of raw results
# ============================================================

suppressPackageStartupMessages({
  library(satuRn)
  library(SummarizedExperiment)
  library(BiocParallel)
})


build_se <- function(counts_mat, meta, formula) {
  # counts_mat: transcripts × samples
  # meta: data.frame, rownames = sample IDs matching colnames(counts_mat)
  meta <- meta[colnames(counts_mat), , drop = FALSE]

  # Encode factors with explicit reference levels
  meta$condition <- factor(meta$condition, levels = COND_LEVELS)
  meta$sex       <- factor(meta$sex)

  # rowData requires gene_id and isoform_id (mandatory for fitDTU)
  tx_ids   <- rownames(counts_mat)
  gene_ids <- sub("-[^-]+$", "", tx_ids)  # GENE-NNN -> GENE
  row_data <- DataFrame(
    isoform_id = tx_ids,
    gene_id    = gene_ids,
    row.names  = tx_ids
  )

  SummarizedExperiment(
    assays  = list(counts = as.matrix(counts_mat)),
    colData = DataFrame(meta),
    rowData = row_data
  )
}

fit_and_test <- function(counts_mat, meta, formula, contrast_name = "AD_vs_Control",
                         diagplots = FALSE, verbose = TRUE) {
  filter_result     <- filter_counts(counts_mat, verbose = verbose)
  counts_structural <- filter_result$counts

  # Drop samples with NA in any covariate used by formula (applied here so
  # counts_structural and the fit set stay column-aligned)
  formula_vars <- all.vars(formula)
  keep_samps   <- complete.cases(meta[, formula_vars, drop = FALSE])
  if (sum(keep_samps) < ncol(counts_structural)) {
    message(sprintf("  Dropping %d sample(s) with NA covariates",
                    ncol(counts_structural) - sum(keep_samps)))
    counts_structural <- counts_structural[, keep_samps, drop = FALSE]
    meta              <- meta[keep_samps, , drop = FALSE]
  }

  # Lighter expression filter restricts which transcripts are actually
  # fit/tested (see step11's filter_counts_for_fit) — counts_structural is
  # returned below for PSI denominators, which stay unbiased by this filter.
  counts_filt <- filter_counts_for_fit(counts_structural, verbose = verbose)

  if (nrow(counts_filt) < 2) {
    warning("  Too few transcripts after filtering — skipping fit")
    return(NULL)
  }

  se <- build_se(counts_filt, meta, formula)

  # Design matrix — used to build the contrast vector
  meta_coldata <- as.data.frame(colData(se))
  meta_coldata$condition <- factor(meta_coldata$condition, levels = COND_LEVELS)
  meta_coldata$sex       <- factor(meta_coldata$sex)
  design <- model.matrix(formula, data = meta_coldata)

  if (!"conditionAD" %in% colnames(design)) {
    stop("conditionAD not found in design matrix columns: ",
         paste(colnames(design), collapse = ", "))
  }

  # Contrast: test conditionAD coefficient
  L <- matrix(0, nrow = ncol(design), ncol = 1,
               dimnames = list(colnames(design), contrast_name))
  L["conditionAD", contrast_name] <- 1

  # Fit
  if (verbose) message("  Fitting DTU (satuRn)...")
  bpparam <- MulticoreParam(N_CORES)
  se <- fitDTU(object = se, formula = formula, parallel = TRUE, BPPARAM = bpparam,
               verbose = verbose)

  # Test
  if (verbose) message("  Testing DTU...")
  se <- testDTU(object = se, contrasts = L, diagplot1 = diagplots, diagplot2 = diagplots,
                sort = FALSE)

  # Extract results into a data.frame
  result_slot <- paste0("fitDTUResult_", contrast_name)
  res <- as.data.frame(rowData(se)[[result_slot]])
  res$transcript_id <- rownames(res)
  res$gene_id       <- sub("-[^-]+$", "", res$transcript_id)

  fit_stats <- list(
    n_tx    = nrow(counts_filt),
    n_genes = length(unique(sub("-[^-]+$", "", rownames(counts_filt))))
  )

  list(se = se, results = res, design = design, filter_stats = filter_result,
       fit_stats = fit_stats, counts_structural = counts_structural)
}

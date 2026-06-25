#!/usr/bin/env Rscript
# ============================================================
# Step 16: Covariate diagnostics for the primary DTU model
# ============================================================
# Investigates whether each covariate in FORMULA_PRIMARY
# (~ condition + age + sex + median_pct_mt) is earning its degree
# of freedom, as a way to look for legitimate DTU testing power
# beyond just relaxing the significance/filtering thresholds.
#
# Three tests, each with its own plot(s):
#
#   A. Per-transcript covariate significance — refits the full
#      primary model once per cell type and tests EVERY
#      coefficient (not just conditionAD) via satuRn's own
#      testDTU() contrast mechanism. Reports what fraction of
#      transcripts have a significant (p < 0.05) coefficient for
#      age / sex / median_pct_mt. A covariate that's rarely
#      significant is "expensive" (costs a DF, shrinks power for
#      conditionAD) without buying explanatory power.
#
#   B. Reduced-model hit comparison — reruns the full pipeline
#      (filter -> fitDTU -> testDTU -> stageR) with each covariate
#      dropped one at a time, holding the sample set fixed across
#      variants, and compares the number of significant DTU hits
#      to the full model.
#
#   C. Covariate/condition confounding check — tests whether each
#      covariate is balanced between AD and Control donors. An
#      imbalanced covariate is a real confounder: dropping it is
#      risky even if (A) says it's rarely significant and (B)
#      says it doesn't cost hits.
#
# Standalone diagnostic script — not wired into run_layer1.R, does
# not modify any layer1 results. Mirrors step15_gene_de.R's pattern.
#
# Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/Rscript layer1/step16_covariate_diagnostics.R
#      (from /home/welcome3/Viscacha_pipeline)
# Output: outputs/layer1/covariate_diagnostics/
#   covariate_pvalues_all_celltypes.csv
#   reduced_model_hit_counts.csv
#   covariate_balance.csv
#   covariate_diagnostics_summary.csv
#   plots/covariate_pval_histograms.png
#   plots/covariate_frac_significant.png
#   plots/hit_counts_by_model_variant.png
#   plots/covariate_balance.png
# ============================================================

suppressPackageStartupMessages({
  library(satuRn)
  library(stageR)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
})

options(mc.cores = 4)

# Source layer1 modules using the script's own directory
script_arg  <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg)
script_dir  <- if (length(script_path) > 0 && nchar(script_path) > 0) {
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  "layer1"
}

source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "step11_filter.R"))
source(file.path(script_dir, "step12_fit.R"))
source(file.path(script_dir, "step13_stageR.R"))

DIAG_OUT_DIR <- file.path(OUT_DIR, "covariate_diagnostics")
PLOTS_DIR    <- file.path(DIAG_OUT_DIR, "plots")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Formula variants for Test B: full primary model vs. each covariate dropped
FORMULA_VARIANTS <- list(
  full     = FORMULA_PRIMARY,
  drop_age = ~ condition + sex + median_pct_mt,
  drop_sex = ~ condition + age + median_pct_mt,
  drop_mt  = ~ condition + age + sex
)

THEME_DIAG <- theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "#555555"),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    legend.position  = "bottom"
  )

save_diag_plot <- function(p, name, width = 8, height = 6) {
  path <- file.path(PLOTS_DIR, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  message("  Plot: ", basename(path))
  invisible(path)
}

load_data <- function(cell_type) {
  f_counts <- file.path(IN_DIR, paste0("counts_", cell_type, ".csv"))
  f_meta   <- file.path(IN_DIR, paste0("metadata_", cell_type, ".csv"))
  if (!file.exists(f_counts)) return(NULL)
  counts_raw <- read.csv(f_counts, row.names = 1, check.names = FALSE)
  meta_raw   <- read.csv(f_meta,   row.names = 1, check.names = FALSE)
  counts_mat <- t(as.matrix(counts_raw))
  storage.mode(counts_mat) <- "integer"
  list(counts = counts_mat, meta = meta_raw)
}

# ------------------------------------------------------------
# Test A: per-transcript covariate significance
# ------------------------------------------------------------
# Tests every non-intercept coefficient of FORMULA_PRIMARY using
# satuRn's own testDTU() contrast mechanism, so the p-values are
# computed exactly the way the main pipeline computes them for
# conditionAD.
test_covariate_significance <- function(counts_mat, meta, cell_type) {
  filter_result <- filter_counts(counts_mat, verbose = FALSE)
  counts_filt <- filter_result$counts
  if (nrow(counts_filt) < 2) return(NULL)

  se <- build_se(counts_filt, meta, FORMULA_PRIMARY)
  meta_coldata <- as.data.frame(colData(se))
  meta_coldata$condition <- factor(meta_coldata$condition, levels = COND_LEVELS)
  meta_coldata$sex       <- factor(meta_coldata$sex)
  design <- model.matrix(FORMULA_PRIMARY, data = meta_coldata)

  terms_to_test <- setdiff(colnames(design), "(Intercept)")
  L <- matrix(0, nrow = ncol(design), ncol = length(terms_to_test),
              dimnames = list(colnames(design), terms_to_test))
  for (term in terms_to_test) L[term, term] <- 1

  se <- suppressWarnings(fitDTU(object = se, formula = FORMULA_PRIMARY, parallel = TRUE,
               BPPARAM = MulticoreParam(N_CORES), verbose = FALSE))
  se <- testDTU(object = se, contrasts = L, diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE)

  out <- lapply(terms_to_test, function(term) {
    res <- as.data.frame(rowData(se)[[paste0("fitDTUResult_", term)]])
    data.frame(
      cell_type     = cell_type,
      transcript_id = rownames(res),
      term          = term,
      estimate      = res$estimates,
      pval          = res$pval,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

summarize_covariate_significance <- function(pval_df) {
  pval_df %>%
    filter(term != "conditionAD") %>%
    group_by(cell_type, term) %>%
    summarise(
      n_tested = sum(!is.na(pval)),
      n_sig    = sum(pval < 0.05, na.rm = TRUE),
      frac_sig = n_sig / n_tested,
      .groups  = "drop"
    )
}

# ------------------------------------------------------------
# Test B: reduced-model hit comparison
# ------------------------------------------------------------
run_model_variant <- function(counts_mat, meta, formula, cell_type, label) {
  fit_result <- tryCatch(
    fit_and_test(counts_mat, meta, formula, contrast_name = "AD_vs_Control", verbose = FALSE),
    error = function(e) { message("    ERROR (", label, "): ", conditionMessage(e)); NULL }
  )
  if (is.null(fit_result)) {
    return(data.frame(cell_type = cell_type, model = label,
                       n_tx_tested = 0, n_sig_genes = 0, n_sig_tx = 0))
  }

  res_df  <- fit_result$results
  padj_df <- tryCatch(run_stageR(res_df, alpha = ALPHA), error = function(e) NULL)
  if (is.null(padj_df)) {
    return(data.frame(cell_type = cell_type, model = label,
                       n_tx_tested = nrow(res_df), n_sig_genes = 0, n_sig_tx = 0))
  }

  sig_mask <- !is.na(padj_df$gene) & padj_df$gene < ALPHA &
              !is.na(padj_df$transcript) & padj_df$transcript < ALPHA
  data.frame(
    cell_type   = cell_type,
    model       = label,
    n_tx_tested = nrow(res_df),
    n_sig_genes = length(unique(padj_df$geneID[sig_mask])),
    n_sig_tx    = sum(sig_mask, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# Test C: covariate/condition confounding check
# ------------------------------------------------------------
test_covariate_balance <- function(meta, cell_type) {
  meta$condition <- factor(meta$condition, levels = COND_LEVELS)

  age_test <- tryCatch(t.test(age ~ condition, data = meta), error = function(e) NULL)
  mt_test  <- tryCatch(t.test(median_pct_mt ~ condition, data = meta), error = function(e) NULL)
  sex_test <- tryCatch(fisher.test(table(meta$sex, meta$condition)), error = function(e) NULL)

  data.frame(
    cell_type = cell_type,
    covariate = c("age", "median_pct_mt", "sex"),
    test      = c("Welch t-test", "Welch t-test", "Fisher's exact"),
    pval      = c(
      if (!is.null(age_test)) age_test$p.value else NA,
      if (!is.null(mt_test))  mt_test$p.value  else NA,
      if (!is.null(sex_test)) sex_test$p.value else NA
    ),
    n_samples = nrow(meta),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
plot_covariate_pval_hist <- function(pval_df) {
  p <- pval_df %>%
    filter(term != "conditionAD") %>%
    ggplot(aes(x = pval)) +
    geom_histogram(bins = 30, fill = "#4292c6", boundary = 0) +
    facet_grid(term ~ cell_type) +
    labs(title = "Covariate coefficient p-value distributions",
         subtitle = "A flat/uniform distribution means the covariate is mostly noise across transcripts",
         x = "p-value", y = "count") +
    THEME_DIAG
  save_diag_plot(p, "covariate_pval_histograms", width = 12, height = 7)
}

plot_covariate_frac_sig <- function(summary_df) {
  p <- summary_df %>%
    ggplot(aes(x = term, y = frac_sig, fill = cell_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(title = "Fraction of transcripts where covariate coefficient is significant (p < 0.05)",
         subtitle = "Low fraction across all cell types = covariate is costing power without buying signal",
         x = NULL, y = "fraction significant") +
    THEME_DIAG
  save_diag_plot(p, "covariate_frac_significant")
}

plot_hit_counts <- function(hits_df) {
  p <- hits_df %>%
    ggplot(aes(x = model, y = n_sig_tx, fill = cell_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(title = "Significant DTU transcripts by model variant",
         subtitle = "full = primary model; drop_X = primary model with covariate X removed",
         x = NULL, y = "# significant transcripts (padj_gene & padj_tx < ALPHA)") +
    THEME_DIAG
  save_diag_plot(p, "hit_counts_by_model_variant")
}

plot_covariate_balance <- function(balance_df) {
  p <- balance_df %>%
    ggplot(aes(x = covariate, y = -log10(pval), color = cell_type)) +
    geom_jitter(size = 3, width = 0.15, height = 0) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = "Covariate balance between AD and Control donors",
         subtitle = "Points above the dashed line are imbalanced (p<0.05) -- real confounders, risky to drop",
         x = NULL, y = "-log10(p-value)") +
    THEME_DIAG
  save_diag_plot(p, "covariate_balance")
}

# ============================================================
# Main loop
# ============================================================
all_pvals   <- list()
all_hits    <- list()
all_balance <- list()

cat(strrep("=", 60), "\n")
cat("Viscacha Layer 1 — covariate diagnostics\n")
cat(strrep("=", 60), "\n")

for (ct in CELL_TYPES) {
  cat(sprintf("\n%s\n%s\n", strrep("-", 60), ct))

  dat <- load_data(ct)
  if (is.null(dat)) {
    message("  Pseudobulk not found — skipping")
    next
  }
  counts_mat <- dat$counts
  meta       <- dat$meta

  # Fix the sample set across all model variants (based on the FULL formula's
  # covariates) so Test B compares model variants on identical data — the only
  # difference between variants should be the design matrix, not sample count.
  formula_vars <- all.vars(FORMULA_PRIMARY)
  keep_samps   <- complete.cases(meta[, formula_vars, drop = FALSE])
  counts_mat   <- counts_mat[, keep_samps, drop = FALSE]
  meta         <- meta[keep_samps, , drop = FALSE]

  message("  [Test A] covariate coefficient significance...")
  pval_df <- tryCatch(test_covariate_significance(counts_mat, meta, ct),
                       error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
  if (!is.null(pval_df)) all_pvals[[ct]] <- pval_df

  message("  [Test B] reduced-model hit counts...")
  for (label in names(FORMULA_VARIANTS)) {
    all_hits[[paste(ct, label)]] <- run_model_variant(
      counts_mat, meta, FORMULA_VARIANTS[[label]], ct, label
    )
  }

  message("  [Test C] covariate/condition balance...")
  all_balance[[ct]] <- tryCatch(test_covariate_balance(meta, ct),
                                 error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL })
}

# ============================================================
# Combine, write CSVs, plot
# ============================================================
pval_combined    <- do.call(rbind, all_pvals)
hits_combined    <- do.call(rbind, all_hits)
balance_combined <- do.call(rbind, all_balance)

write.csv(pval_combined, file.path(DIAG_OUT_DIR, "covariate_pvalues_all_celltypes.csv"), row.names = FALSE)
write.csv(hits_combined, file.path(DIAG_OUT_DIR, "reduced_model_hit_counts.csv"), row.names = FALSE)
write.csv(balance_combined, file.path(DIAG_OUT_DIR, "covariate_balance.csv"), row.names = FALSE)

sig_summary <- summarize_covariate_significance(pval_combined)

plot_covariate_pval_hist(pval_combined)
plot_covariate_frac_sig(sig_summary)
plot_hit_counts(hits_combined)
plot_covariate_balance(balance_combined)

# ------------------------------------------------------------
# Decision summary: combine all three tests per covariate
# ------------------------------------------------------------
full_hits <- hits_combined %>% filter(model == "full") %>%
  summarise(total = sum(n_sig_tx)) %>% pull(total)

decision_rows <- lapply(c("age", "sex", "median_pct_mt"), function(cov) {
  drop_label <- paste0("drop_", if (cov == "median_pct_mt") "mt" else cov)
  term_name  <- if (cov == "sex") "sexM" else cov

  drop_hits <- hits_combined %>% filter(model == drop_label) %>%
    summarise(total = sum(n_sig_tx)) %>% pull(total)

  mean_frac_sig <- sig_summary %>% filter(term == term_name) %>%
    summarise(m = mean(frac_sig, na.rm = TRUE)) %>% pull(m)

  mean_balance_pval <- balance_combined %>% filter(covariate == cov) %>%
    summarise(m = mean(pval, na.rm = TRUE)) %>% pull(m)

  data.frame(
    covariate                 = cov,
    mean_frac_transcripts_sig = mean_frac_sig,
    mean_balance_pval         = mean_balance_pval,
    hits_full_model           = full_hits,
    hits_with_covariate_dropped = drop_hits,
    hit_delta                 = drop_hits - full_hits,
    recommendation            = ifelse(
      mean_balance_pval < 0.05,
      "KEEP — imbalanced between AD/Control (real confounder)",
      ifelse(drop_hits > full_hits,
             "CONSIDER DROPPING — balanced and dropping gains hits",
             "KEEP — balanced but dropping doesn't gain hits")
    )
  )
})

decision_df <- do.call(rbind, decision_rows)
write.csv(decision_df, file.path(DIAG_OUT_DIR, "covariate_diagnostics_summary.csv"), row.names = FALSE)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Covariate diagnostics summary\n")
cat(strrep("=", 60), "\n")
print(decision_df, row.names = FALSE)

cat("\nOutputs written to: ", DIAG_OUT_DIR, "\n", sep = "")
cat(strrep("=", 60), "\n")
cat("Viscacha Layer 1 covariate diagnostics — COMPLETE\n")
cat(strrep("=", 60), "\n")

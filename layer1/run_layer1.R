#!/usr/bin/env Rscript
# ============================================================
# Viscacha Layer 1 — orchestrator
# Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/Rscript layer1/run_layer1.R
# (from /home/welcome3/Viscacha_pipeline)
# ============================================================
# Option C: two models per cell type
#   Primary:  ~ condition + age + sex + median_pct_mt  (main results)
#   Braak:    + braak_stage                            (sensitivity / robustness)
# Significant in both  → robust_to_braak = TRUE
# ============================================================

suppressPackageStartupMessages({
  library(satuRn)
  library(stageR)
  library(SummarizedExperiment)
  library(BiocParallel)
})

options(mc.cores = 8)

# Source layer1 modules using the script's own directory
script_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
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
source(file.path(script_dir, "step14_sensitivity.R"))
source(file.path(script_dir, "plots.R"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- QC log accumulator (base R, jsonlite not available) ---
qc_entries <- list()

log_qc <- function(ct, model, n_tx_before, n_tx_after, n_genes, n_sig_genes, n_sig_tx, note = "") {
  qc_entries[[length(qc_entries) + 1]] <<- list(
    cell_type    = ct,
    model        = model,
    n_tx_before  = n_tx_before,
    n_tx_after   = n_tx_after,
    n_genes      = n_genes,
    n_sig_genes  = n_sig_genes,
    n_sig_tx     = n_sig_tx,
    note         = note
  )
}

write_qc_log <- function(path) {
  lines <- c("# Viscacha Layer 1 QC Log", "")
  for (e in qc_entries) {
    lines <- c(lines,
      sprintf("cell_type:   %s", e$cell_type),
      sprintf("model:       %s", e$model),
      sprintf("tx_before:   %d", e$n_tx_before),
      sprintf("tx_after:    %d", e$n_tx_after),
      sprintf("n_genes:     %d", e$n_genes),
      sprintf("sig_genes:   %d", e$n_sig_genes),
      sprintf("sig_tx:      %d", e$n_sig_tx),
      if (nchar(e$note) > 0) sprintf("note:        %s", e$note) else NULL,
      "---"
    )
  }
  writeLines(lines, path)
  message("QC log written to: ", path)
}

# --- Helpers ---
load_data <- function(cell_type, suffix = "") {
  stem <- cell_type
  f_counts <- file.path(IN_DIR, paste0("counts_", stem, suffix, ".csv"))
  f_meta   <- file.path(IN_DIR, paste0("metadata_", stem, suffix, ".csv"))
  if (!file.exists(f_counts)) return(NULL)

  counts_raw <- read.csv(f_counts, row.names = 1, check.names = FALSE)
  meta_raw   <- read.csv(f_meta,   row.names = 1, check.names = FALSE)

  # counts CSV is samples × transcripts; SE needs transcripts × samples
  counts_mat <- t(as.matrix(counts_raw))
  storage.mode(counts_mat) <- "integer"

  list(counts = counts_mat, meta = meta_raw)
}

run_one_model <- function(counts_mat, meta, formula, cell_type, model_label) {
  message(sprintf("\n  [%s | %s]", cell_type, model_label))

  n_tx_before <- nrow(counts_mat)

  fit_result <- tryCatch(
    fit_and_test(counts_mat, meta, formula, contrast_name = "AD_vs_Control",
                 diagplots = FALSE, verbose = TRUE),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }
  )

  if (is.null(fit_result)) {
    log_qc(cell_type, model_label, n_tx_before, 0, 0, 0, 0, "fit failed")
    return(NULL)
  }

  res_df  <- fit_result$results
  fstats  <- fit_result$filter_stats

  # stageR FDR
  stageR_df <- tryCatch(
    run_stageR(res_df, alpha = ALPHA),
    error = function(e) { message("  stageR ERROR: ", conditionMessage(e)); NULL }
  )

  if (is.null(stageR_df)) {
    log_qc(cell_type, model_label, n_tx_before, fstats$n_tx, fstats$n_genes, 0, 0,
           "stageR failed")
    return(NULL)
  }

  # ΔPSI (from filtered counts in the SE)
  counts_filt <- assay(fit_result$se, "counts")
  meta_filt   <- as.data.frame(colData(fit_result$se))
  dpsi_df     <- compute_delta_psi(counts_filt, meta_filt)

  # Merge all result components
  final_df <- merge_results(stageR_df, dpsi_df, res_df)

  # Significant = gene-level padj < ALPHA AND transcript-level padj < ALPHA
  sig_mask <- !is.na(final_df$padj_gene) & final_df$padj_gene < ALPHA &
              !is.na(final_df$padj_tx)   & final_df$padj_tx   < ALPHA
  n_sig_tx    <- sum(sig_mask, na.rm = TRUE)
  n_sig_genes <- length(unique(final_df$gene_id[sig_mask]))

  log_qc(cell_type, model_label, n_tx_before, fstats$n_tx, fstats$n_genes,
         n_sig_genes, n_sig_tx)
  message(sprintf("  Significant: %d genes, %d transcripts", n_sig_genes, n_sig_tx))

  list(results = final_df, filter_stats = fstats)
}

# ============================================================
# Main loop
# ============================================================
all_results      <- list()
filter_stats_all <- list()

cat(strrep("=", 60), "\n")
cat("Viscacha Layer 1 — starting\n")
cat(strrep("=", 60), "\n")

for (ct in CELL_TYPES) {
  cat(sprintf("\n%s\n%s\n", strrep("-", 60), ct))

  # Load primary data (AD + Control)
  dat <- load_data(ct)
  if (is.null(dat)) {
    message("  Primary data not found — skipping")
    next
  }
  counts_mat <- dat$counts
  meta       <- dat$meta

  # ---- Primary model ----
  primary_result <- run_one_model(counts_mat, meta, FORMULA_PRIMARY, ct, "primary")
  primary_df     <- if (!is.null(primary_result)) primary_result$results else NULL
  if (!is.null(primary_result))
    filter_stats_all[[ct]] <- filter_counts_breakdown(
      counts_mat, MIN_TX_COUNT, MIN_GENE_COUNT, MIN_SAMPS_FRAC
    )

  # ---- Braak sensitivity model ----
  # Drop donors missing braak_stage
  braak_ok   <- !is.na(meta$braak_stage)
  if (sum(braak_ok) < ncol(counts_mat)) {
    message(sprintf("  Braak model: dropping %d donor(s) with NA braak_stage",
                    ncol(counts_mat) - sum(braak_ok)))
  }
  counts_braak <- counts_mat[, braak_ok, drop = FALSE]
  meta_braak   <- meta[braak_ok, , drop = FALSE]

  braak_result <- run_one_model(counts_braak, meta_braak, FORMULA_BRAAK, ct, "braak_sensitivity")
  braak_df     <- if (!is.null(braak_result)) braak_result$results else NULL

  # ---- Mark robustness ----
  if (!is.null(primary_df) && !is.null(braak_df)) {
    sig_primary <- primary_df$transcript_id[
      !is.na(primary_df$padj_gene) & primary_df$padj_gene < ALPHA &
      !is.na(primary_df$padj_tx)   & primary_df$padj_tx   < ALPHA
    ]
    sig_braak <- braak_df$transcript_id[
      !is.na(braak_df$padj_gene) & braak_df$padj_gene < ALPHA &
      !is.na(braak_df$padj_tx)   & braak_df$padj_tx   < ALPHA
    ]
    primary_df$robust_to_braak <- primary_df$transcript_id %in% sig_braak &
                                   primary_df$transcript_id %in% sig_primary
  } else if (!is.null(primary_df)) {
    primary_df$robust_to_braak <- NA
  }

  # ---- Active Control descriptive PSI ----
  dat_sens <- load_data(ct, suffix = "_sensitivity")
  if (!is.null(dat_sens) && !is.null(primary_df)) {
    # Subset sensitivity counts to transcripts in primary results
    tx_keep       <- intersect(rownames(dat_sens$counts), primary_df$transcript_id)
    counts_active <- dat_sens$counts[tx_keep, , drop = FALSE]
    meta_active   <- dat_sens$meta

    active_psi <- compute_active_psi(counts_active, meta_active)
    primary_df <- merge_active_psi(primary_df, active_psi)
  }

  # ---- Save per-cell-type results ----
  if (!is.null(primary_df)) {
    primary_df$cell_type <- ct
    out_path <- file.path(OUT_DIR, paste0("dtu_results_", ct, ".csv"))
    write.csv(primary_df, out_path, row.names = FALSE)
    message(sprintf("  Saved: %s", basename(out_path)))
    all_results[[ct]] <- primary_df
  }
}

# ============================================================
# Combined output: significant transcripts only
# ============================================================
if (length(all_results) > 0) {
  combined <- do.call(rbind, all_results)
  sig_combined <- combined[
    !is.na(combined$padj_gene) & combined$padj_gene < ALPHA &
    !is.na(combined$padj_tx)   & combined$padj_tx   < ALPHA,
  ]

  sig_path <- file.path(OUT_DIR, "dtu_significant_all_celltypes.csv")
  write.csv(sig_combined, sig_path, row.names = FALSE)
  message(sprintf("\nCombined significant results: %d transcripts across %d genes",
                  nrow(sig_combined),
                  length(unique(sig_combined$gene_id))))
  message("Saved: ", basename(sig_path))
}

# ============================================================
# Plots
# ============================================================
sig_combined_for_plots <- if (exists("sig_combined")) sig_combined else NULL
generate_all_plots(all_results, filter_stats_all, sig_combined_for_plots)

# ============================================================
# QC log
# ============================================================
dir.create(file.path(OUT_DIR, "qc"), showWarnings = FALSE)
write_qc_log(file.path(OUT_DIR, "qc", "layer1_qc.txt"))

cat(strrep("=", 60), "\n")
cat("Viscacha Layer 1 — COMPLETE\n")
cat(strrep("=", 60), "\n")

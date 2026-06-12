#!/usr/bin/env Rscript
# ============================================================
# Viscacha Layer 1 — standalone plot regeneration
# Loads existing per-cell-type result CSVs and regenerates all plots.
# Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/Rscript layer1/make_plots.R
# (from /home/welcome3/Viscacha_pipeline)
# ============================================================

script_arg  <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg)
script_dir  <- if (length(script_path) > 0 && nchar(script_path) > 0) {
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  "layer1"
}

source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "step11_filter.R"))
source(file.path(script_dir, "plots.R"))

# ---- Load per-cell-type results from existing CSVs ----
all_results <- list()
for (ct in CELL_TYPES) {
  path <- file.path(OUT_DIR, paste0("dtu_results_", ct, ".csv"))
  if (!file.exists(path)) {
    message("  Missing: ", basename(path), " — skipping")
    next
  }
  all_results[[ct]] <- read.csv(path, stringsAsFactors = FALSE)
  message("  Loaded: ", basename(path),
          "  (", nrow(all_results[[ct]]), " transcripts)")
}

if (length(all_results) == 0) {
  stop("No result CSVs found in ", OUT_DIR, " — run run_layer1.R first")
}

# ---- Recompute filter stats from raw counts ----
# (needed for step11 filter retention plot)
message("\nRecomputing filter stats from raw counts...")
filter_stats_all <- list()
for (ct in names(all_results)) {
  f <- file.path(IN_DIR, paste0("counts_", ct, ".csv"))
  if (!file.exists(f)) next
  counts_raw <- read.csv(f, row.names = 1, check.names = FALSE)
  counts_mat <- t(as.matrix(counts_raw))
  storage.mode(counts_mat) <- "integer"
  filter_stats_all[[ct]] <- filter_counts_breakdown(
    counts_mat, MIN_TX_COUNT, MIN_GENE_COUNT, MIN_SAMPS_FRAC
  )
  s <- filter_stats_all[[ct]]
  message(sprintf("  %-22s  before=%6d  passed=%6d  (low-tx=%d  low-gene=%d  single=%d)",
                  ct, s$n_start, s$n_passed,
                  s$n_failed_tx_count, s$n_failed_gene_count, s$n_single_isoform))
}

# ---- Build significant results table ----
combined <- do.call(rbind, all_results)
sig_combined <- combined[
  !is.na(combined$padj_gene) & combined$padj_gene < ALPHA &
  !is.na(combined$padj_tx)   & combined$padj_tx   < ALPHA,
]
message(sprintf("\nSignificant: %d transcripts in %d genes across %d cell type(s)",
                nrow(sig_combined),
                length(unique(sig_combined$gene_id)),
                length(unique(sig_combined$cell_type))))

# ---- Generate all plots ----
generate_all_plots(all_results, filter_stats_all,
                   if (nrow(sig_combined) > 0) sig_combined else NULL)

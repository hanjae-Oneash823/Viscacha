#!/usr/bin/env Rscript
# ============================================================
# Step 15 (GX variant): Gene-level differential expression (DESeq2)
# using the LONG-READ GENE-LEVEL pseudobulk (outputs/00_PreAggregation_QC/pseudobulk_gx/,
# built from adata_gene_loose_filtering_for_bulk_analysis.h5ad directly —
# no transcript->gene collapsing), for comparison against
# step15_gene_de_sr.R (short-read gene-level).
#
# Same primary model as the DTU analysis and the other step15 variants
# (~ condition + age + sex + median_pct_mt), Control as reference,
# primary donors only (AD + Control).
#
# Run: /home/welcome3/anaconda3/envs/oneash_dtu/bin/Rscript 01_ViscachaDTU_Analysis/step15_gene_de_gx.R
#      (from /home/welcome3/Viscacha_pipeline)
# Output: outputs/01_ViscachaDTU_Analysis/gene_level_de_results_GX.csv
#   columns: gene_id, cell_type, log2FC, lfcSE, pval, padj, baseMean
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

script_arg  <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg)
script_dir  <- if (length(script_path) > 0 && nchar(script_path) > 0) {
  dirname(normalizePath(script_path, mustWork = FALSE))
} else { "01_ViscachaDTU_Analysis" }
source(file.path(script_dir, "config.R"))

IN_DIR_GX <- "/home/welcome3/Viscacha_pipeline/outputs/00_PreAggregation_QC/pseudobulk_gx"
MIN_GENE_TOTAL <- 10

run_gene_de_one <- function(cell_type) {
  f_counts <- file.path(IN_DIR_GX, paste0("counts_", cell_type, ".csv"))
  f_meta   <- file.path(IN_DIR_GX, paste0("metadata_", cell_type, ".csv"))
  if (!file.exists(f_counts) || !file.exists(f_meta)) {
    message("  [", cell_type, "] GX pseudobulk not found — skipping")
    return(NULL)
  }

  counts_raw <- read.csv(f_counts, row.names = 1, check.names = FALSE)
  meta       <- read.csv(f_meta,   row.names = 1, check.names = FALSE)

  # samples x genes -> genes x samples (already gene-level, no collapsing)
  gene_mat <- t(as.matrix(counts_raw))
  storage.mode(gene_mat) <- "integer"

  gene_mat <- gene_mat[rowSums(gene_mat) >= MIN_GENE_TOTAL, , drop = FALSE]

  meta <- meta[colnames(gene_mat), , drop = FALSE]
  model_vars <- all.vars(FORMULA_PRIMARY)
  keep <- complete.cases(meta[, model_vars, drop = FALSE])
  if (sum(keep) < ncol(gene_mat)) {
    message(sprintf("  [%s] dropping %d donor(s) with NA covariates",
                    cell_type, ncol(gene_mat) - sum(keep)))
  }
  gene_mat <- gene_mat[, keep, drop = FALSE]
  meta     <- meta[keep, , drop = FALSE]

  meta$condition <- factor(meta$condition, levels = COND_LEVELS)
  meta$sex       <- factor(meta$sex)

  for (v in all.vars(FORMULA_PRIMARY)) {
    if (is.numeric(meta[[v]])) meta[[v]] <- as.numeric(scale(meta[[v]]))
  }
  if (nlevels(droplevels(meta$condition)) < 2 || ncol(gene_mat) < 4) {
    message("  [", cell_type, "] insufficient samples/conditions — skipping")
    return(NULL)
  }

  dds <- DESeqDataSetFromMatrix(countData = gene_mat, colData = meta,
                                design = FORMULA_PRIMARY)
  dds <- tryCatch(DESeq(dds, quiet = TRUE),
                  error = function(e) { message("  [", cell_type, "] DESeq error: ",
                                                conditionMessage(e)); NULL })
  if (is.null(dds)) return(NULL)

  res <- results(dds, name = "condition_AD_vs_Control")
  out <- data.frame(
    gene_id   = rownames(res),
    cell_type = cell_type,
    log2FC    = res$log2FoldChange,
    lfcSE     = res$lfcSE,
    pval      = res$pvalue,
    padj      = res$padj,
    baseMean  = res$baseMean,
    stringsAsFactors = FALSE
  )
  n_sig <- sum(!is.na(out$padj) & out$padj < 0.05)
  message(sprintf("  [%s] %d genes tested, %d DE at padj<0.05",
                  cell_type, nrow(out), n_sig))
  out
}

cat(strrep("=", 60), "\n")
cat("Step 15 (GX) — Gene-level DE (DESeq2) on long-read gene-level pseudobulk\n")
cat(strrep("=", 60), "\n")

all_de <- list()
for (ct in CELL_TYPES) {
  message("\n", ct)
  de <- run_gene_de_one(ct)
  if (!is.null(de)) all_de[[ct]] <- de
}

if (length(all_de) == 0) {
  stop("No gene-level DE results produced.")
}

combined <- do.call(rbind, all_de)
out_path <- file.path(OUT_DIR, "gene_level_de_results_GX.csv")
write.csv(combined, out_path, row.names = FALSE)
cat(sprintf("\nSaved %d gene x cell_type rows to %s\n",
            nrow(combined), out_path))

# ============================================================
# Step 13: stageR two-stage FDR + ΔPSI computation
# ============================================================
# Input:  fit results from step12, raw counts matrix, metadata
# Output: data.frame with stageR-adjusted p-values + ΔPSI
# ============================================================

suppressPackageStartupMessages(library(stageR))


run_stageR <- function(results_df, alpha = ALPHA) {
  # results_df must have columns: transcript_id, gene_id, pval
  res <- results_df[!is.na(results_df$pval), ]

  tx2gene <- data.frame(
    tx   = res$transcript_id,
    gene = res$gene_id,
    stringsAsFactors = FALSE
  )

  # pScreen: gene-level (minimum per-gene transcript p-value)
  # as.numeric() is mandatory — tapply returns a named array which stageR rejects
  pScreen       <- as.numeric(tapply(res$pval, res$gene_id, min))
  names(pScreen) <- names(tapply(res$pval, res$gene_id, min))

  # pConfirmation: transcript-level matrix (one column)
  pConf_vec <- res$pval
  names(pConf_vec) <- res$transcript_id
  pConfirmation <- matrix(pConf_vec, ncol = 1,
                          dimnames = list(names(pConf_vec), "AD_vs_Control"))

  stageObj <- stageRTx(
    pScreen          = pScreen,
    pConfirmation    = pConfirmation,
    pScreenAdjusted  = FALSE,
    tx2gene          = tx2gene
  )
  stageObj <- stageWiseAdjustment(stageObj, method = "dtu", alpha = alpha)

  padj_df <- getAdjustedPValues(stageObj, onlySignificantGenes = FALSE, order = FALSE)
  # padj_df columns: geneID, txID, gene, transcript (adjusted p-values)

  padj_df
}

compute_delta_psi <- function(counts_mat, meta) {
  # counts_mat: transcripts × samples (raw, unfiltered — but already subset to cell type)
  # Returns data.frame: transcript_id, gene_id, psi_AD, psi_ctrl, delta_psi

  ad_samps   <- rownames(meta)[meta$condition == "AD"]
  ctrl_samps <- rownames(meta)[meta$condition == "Control"]

  ad_samps   <- intersect(ad_samps,   colnames(counts_mat))
  ctrl_samps <- intersect(ctrl_samps, colnames(counts_mat))

  tx_names   <- rownames(counts_mat)
  gene_names <- sub("-[^-]+$", "", tx_names)

  # Gene totals per sample
  gene_totals <- rowsum(counts_mat, gene_names)

  # PSI per transcript per sample = count / gene_total (NA if gene_total == 0)
  compute_psi_group <- function(samps) {
    ct  <- counts_mat[, samps, drop = FALSE]
    gt  <- gene_totals[gene_names, samps, drop = FALSE]
    psi <- ct / gt
    # Mean across donors (ignore NA)
    rowMeans(psi, na.rm = TRUE)
  }

  psi_ad   <- compute_psi_group(ad_samps)
  psi_ctrl <- compute_psi_group(ctrl_samps)

  data.frame(
    transcript_id = tx_names,
    gene_id       = gene_names,
    psi_AD        = psi_ad,
    psi_ctrl      = psi_ctrl,
    delta_psi     = psi_ad - psi_ctrl,
    stringsAsFactors = FALSE
  )
}

merge_results <- function(stageR_df, delta_psi_df, raw_results_df) {
  # stageR_df:     geneID, txID, gene, transcript
  # delta_psi_df:  transcript_id, gene_id, psi_AD, psi_ctrl, delta_psi
  # raw_results_df: transcript_id, gene_id, pval, regularizedDispersion, ...

  # Rename stageR columns to match
  colnames(stageR_df)[colnames(stageR_df) == "txID"]       <- "transcript_id"
  colnames(stageR_df)[colnames(stageR_df) == "geneID"]     <- "gene_id"
  colnames(stageR_df)[colnames(stageR_df) == "gene"]       <- "padj_gene"
  colnames(stageR_df)[colnames(stageR_df) == "transcript"] <- "padj_tx"

  merged <- merge(stageR_df, delta_psi_df, by = c("transcript_id", "gene_id"), all.x = TRUE)

  # Add raw satuRn p-values and empirical flag
  if ("empirical_pval" %in% colnames(raw_results_df)) {
    merged <- merge(merged,
                    raw_results_df[, c("transcript_id", "pval", "empirical_pval",
                                       "regular_FDR", "empirical_FDR")],
                    by = "transcript_id", all.x = TRUE)
  } else {
    merged <- merge(merged,
                    raw_results_df[, c("transcript_id", "pval")],
                    by = "transcript_id", all.x = TRUE)
  }

  merged
}

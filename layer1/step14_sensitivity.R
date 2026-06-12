# ============================================================
# Step 14: Active Control descriptive PSI
# ============================================================
# Purely descriptive — no DTU test.
# Computes mean PSI per transcript for Active Control donors
# and compares against AD and Control PSIs from the primary analysis.
# Output: data.frame with psi_active_ctrl, delta_psi_AD_vs_active columns
# ============================================================

source(file.path(dirname(sys.frame(1)$ofile), "config.R"))

compute_active_psi <- function(counts_sensitivity, meta_sensitivity) {
  # counts_sensitivity: transcripts × samples (all conditions in sensitivity file)
  # meta_sensitivity:   data.frame with condition column

  active_samps <- rownames(meta_sensitivity)[
    meta_sensitivity$condition == "Active control"
  ]
  active_samps <- intersect(active_samps, colnames(counts_sensitivity))

  if (length(active_samps) == 0) {
    message("  No Active control samples found — skipping step 14")
    return(NULL)
  }

  tx_names   <- rownames(counts_sensitivity)
  gene_names <- sub("-[^-]+$", "", tx_names)

  ct      <- counts_sensitivity[, active_samps, drop = FALSE]
  gt_mat  <- rowsum(ct, gene_names)[gene_names, , drop = FALSE]
  psi_mat <- ct / gt_mat

  data.frame(
    transcript_id    = tx_names,
    gene_id          = gene_names,
    psi_active_ctrl  = rowMeans(psi_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

merge_active_psi <- function(primary_results, active_psi_df) {
  if (is.null(active_psi_df)) return(primary_results)

  merge(primary_results, active_psi_df[, c("transcript_id", "psi_active_ctrl")],
        by = "transcript_id", all.x = TRUE)
}

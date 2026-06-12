# ============================================================
# Step 11: Count-based transcript filter (dmFilter equivalent)
# ============================================================
# Input:  counts matrix (samples × transcripts) + metadata
# Output: filtered counts matrix + log summary
#
# Filters applied in order:
#   1. Transcript-level: count >= MIN_TX_COUNT in >= floor(n * MIN_SAMPS_FRAC) samples
#   2. Gene-level:       gene total >= MIN_GENE_COUNT in >= floor(n * MIN_SAMPS_FRAC) samples
#   3. Multi-transcript: keep only genes with >= 2 passing transcripts
# ============================================================

filter_counts <- function(counts_mat, min_tx_count, min_gene_count, min_samps_frac, verbose = TRUE) {
  # counts_mat: transcripts × samples (rows = features, cols = samples)
  n_samples <- ncol(counts_mat)
  min_samps <- max(2L, floor(n_samples * min_samps_frac))

  n_tx_before <- nrow(counts_mat)

  # 1. Transcript-level count filter
  tx_pass <- rowSums(counts_mat >= min_tx_count) >= min_samps
  counts_mat <- counts_mat[tx_pass, , drop = FALSE]

  # 2. Gene totals per sample
  tx_names   <- rownames(counts_mat)
  gene_names <- sub("-[^-]+$", "", tx_names)  # GENE-NNN -> GENE (remove last -NNN)

  gene_totals <- rowsum(counts_mat, gene_names)  # genes × samples

  gene_pass_mask <- rowSums(gene_totals >= min_gene_count) >= min_samps
  passing_genes  <- rownames(gene_totals)[gene_pass_mask]
  counts_mat     <- counts_mat[gene_names %in% passing_genes, , drop = FALSE]
  gene_names     <- sub("-[^-]+$", "", rownames(counts_mat))

  # 3. Multi-transcript gene filter
  tx_per_gene  <- table(gene_names)
  multi_genes  <- names(tx_per_gene)[tx_per_gene >= 2]
  counts_mat   <- counts_mat[gene_names %in% multi_genes, , drop = FALSE]

  n_tx_after   <- nrow(counts_mat)
  n_genes_after <- length(unique(sub("-[^-]+$", "", rownames(counts_mat))))

  if (verbose) {
    message(sprintf("  Transcripts: %d -> %d  |  Genes: %d  |  min_samps=%d/%d",
                    n_tx_before, n_tx_after, n_genes_after, min_samps, n_samples))
  }

  list(
    counts      = counts_mat,
    n_tx        = n_tx_after,
    n_tx_before = n_tx_before,
    n_genes     = n_genes_after,
    min_samps   = min_samps
  )
}

# Returns per-stage exclusion counts for the step 11 retention plot.
# Runs the same three stages as filter_counts but records how many transcripts
# each stage removes, so plots can show WHY transcripts were excluded.
filter_counts_breakdown <- function(counts_mat, min_tx_count, min_gene_count, min_samps_frac) {
  n_samples <- ncol(counts_mat)
  min_samps <- max(2L, floor(n_samples * min_samps_frac))
  n_start   <- nrow(counts_mat)

  # Stage 1: transcript count threshold
  tx_pass1 <- rowSums(counts_mat >= min_tx_count) >= min_samps
  n_failed_tx_count <- sum(!tx_pass1)
  mat1 <- counts_mat[tx_pass1, , drop = FALSE]

  # Stage 2: gene total threshold
  gene_ids1   <- sub("-[^-]+$", "", rownames(mat1))
  gene_totals <- rowsum(mat1, gene_ids1)
  gene_pass   <- rowSums(gene_totals >= min_gene_count) >= min_samps
  ok_genes    <- rownames(gene_totals)[gene_pass]
  tx_pass2    <- gene_ids1 %in% ok_genes
  n_failed_gene_count <- sum(!tx_pass2)
  mat2 <- mat1[tx_pass2, , drop = FALSE]

  # Stage 3: multi-transcript gene filter
  gene_ids2   <- sub("-[^-]+$", "", rownames(mat2))
  tx_per_gene <- table(gene_ids2)
  multi_genes <- names(tx_per_gene)[tx_per_gene >= 2]
  tx_pass3    <- gene_ids2 %in% multi_genes
  n_single_isoform <- sum(!tx_pass3)
  n_passed <- sum(tx_pass3)

  list(
    n_start             = n_start,
    n_failed_tx_count   = n_failed_tx_count,
    n_failed_gene_count = n_failed_gene_count,
    n_single_isoform    = n_single_isoform,
    n_passed            = n_passed,
    min_samps           = min_samps,
    n_samples           = n_samples
  )
}

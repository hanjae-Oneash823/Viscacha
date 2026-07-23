# ============================================================
# Step 11: Structural transcript cleanup (no count thresholds)
# ============================================================
# Input:  counts matrix (samples × transcripts) + metadata
# Output: cleaned counts matrix + log summary
#
# This is the set used for gene-level PSI denominators — no minimum-count/
# prevalence filtering here, so totals reflect true total transcript
# expression. Only two structurally-required removals are applied, in order:
#   1. Transcripts with zero counts in every sample (no information to fit;
#      contribute nothing to any gene's total either way).
#   2. Genes left with < 2 transcripts (DTU is undefined for a single isoform).
#
# A separate, lighter expression filter (filter_counts_for_fit, below) is
# applied downstream to pick the subset that's actually fit/tested by
# satuRn+stageR — fitting on every structurally-valid transcript balloons the
# multiple-testing universe with near-zero-count transcripts and crushes
# power, without affecting the denominators computed here.
# ============================================================

filter_counts <- function(counts_mat, verbose = TRUE) {
  # counts_mat: transcripts × samples (rows = features, cols = samples)
  n_tx_before <- nrow(counts_mat)

  # 1. Drop transcripts never detected in any sample
  tx_pass    <- rowSums(counts_mat) > 0
  counts_mat <- counts_mat[tx_pass, , drop = FALSE]

  # 2. Drop genes left with < 2 transcripts
  tx_names    <- rownames(counts_mat)
  gene_names  <- sub("-[^-]+$", "", tx_names)  # GENE-NNN -> GENE (remove last -NNN)
  tx_per_gene <- table(gene_names)
  multi_genes <- names(tx_per_gene)[tx_per_gene >= 2]
  counts_mat  <- counts_mat[gene_names %in% multi_genes, , drop = FALSE]

  n_tx_after    <- nrow(counts_mat)
  n_genes_after <- length(unique(sub("-[^-]+$", "", rownames(counts_mat))))

  if (verbose) {
    message(sprintf("  Transcripts: %d -> %d  |  Genes: %d  (structural cleanup only)",
                    n_tx_before, n_tx_after, n_genes_after))
  }

  list(
    counts  = counts_mat,
    n_tx    = n_tx_after,
    n_genes = n_genes_after
  )
}

# Lighter expression filter, applied on top of filter_counts() output, that
# restricts the set passed into satuRn fitting/testing. Keeps the
# multiple-testing universe (and stageR's gene-level FDR correction) sized to
# transcripts with enough signal to be testable, while filter_counts()'s full
# structural set is still what PSI denominators are computed from.
filter_counts_for_fit <- function(counts_mat, min_count = MIN_TX_COUNT,
                                   min_samps_frac = MIN_SAMPS_FRAC, verbose = TRUE) {
  n_tx_before <- nrow(counts_mat)

  pass_frac  <- rowMeans(counts_mat >= min_count)
  tx_pass    <- pass_frac >= min_samps_frac
  counts_mat <- counts_mat[tx_pass, , drop = FALSE]

  # Re-apply multi-isoform requirement: a gene can drop to 1 isoform here
  # even though it had >= 2 in the structural set.
  tx_names    <- rownames(counts_mat)
  gene_names  <- sub("-[^-]+$", "", tx_names)
  tx_per_gene <- table(gene_names)
  multi_genes <- names(tx_per_gene)[tx_per_gene >= 2]
  counts_mat  <- counts_mat[gene_names %in% multi_genes, , drop = FALSE]

  if (verbose) {
    message(sprintf("  Fit-set filter: %d -> %d transcripts (>=%d counts in >=%.0f%% of samples)",
                    n_tx_before, nrow(counts_mat), min_count, 100 * min_samps_frac))
  }

  counts_mat
}

# Returns per-stage exclusion counts for the step 11 retention plot.
filter_counts_breakdown <- function(counts_mat) {
  n_start <- nrow(counts_mat)

  # Stage 1: never detected
  tx_pass1  <- rowSums(counts_mat) > 0
  n_zero    <- sum(!tx_pass1)
  mat1      <- counts_mat[tx_pass1, , drop = FALSE]

  # Stage 2: single-isoform genes
  gene_ids1        <- sub("-[^-]+$", "", rownames(mat1))
  tx_per_gene      <- table(gene_ids1)
  multi_genes      <- names(tx_per_gene)[tx_per_gene >= 2]
  tx_pass2         <- gene_ids1 %in% multi_genes
  n_single_isoform <- sum(!tx_pass2)
  n_passed         <- sum(tx_pass2)

  list(
    n_start          = n_start,
    n_zero           = n_zero,
    n_single_isoform = n_single_isoform,
    n_passed         = n_passed
  )
}

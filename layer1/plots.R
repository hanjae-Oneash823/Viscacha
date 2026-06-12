# ============================================================
# Viscacha Layer 1 — visualization
# All functions write PNG files to OUT_DIR/plots/
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(pheatmap)
  library(tidyr)
  library(dplyr)
})

# ---- Style constants (match Layer 0 notebook palette) ----
COND_COLORS <- c("AD" = "#d62728", "Control" = "#2ca02c", "Active control" = "#ff7f0e")
COND_ORDER  <- c("Control", "Active control", "AD")
CT_COLORS   <- setNames(
  c("#1f77b4","#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f"),
  c("Excitatory_neuron","Inhibitory_neuron","Oligodendrocyte",
    "OPC","Astrocyte","Microglia","Vascular_cell")
)

THEME_VIS <- theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "#555555"),
    axis.title       = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.4)
  )

save_plot <- function(p, name, width = 8, height = 6) {
  dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(OUT_DIR, "plots", paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  message("  Plot: ", basename(path))
  invisible(path)
}

# ============================================================
# Step 11 — Filter retention with per-stage exclusion breakdown
# filter_stats_all: named list, each element from filter_counts_breakdown()
# ============================================================
plot_step11_filter <- function(filter_stats_all) {
  EXCL_LEVELS <- c("Passed", "Low isoform count", "Low gene total", "Single isoform")
  EXCL_COLORS <- c(
    "Passed"            = "#2171b5",
    "Low isoform count" = "#fd8d3c",
    "Low gene total"    = "#e6550d",
    "Single isoform"    = "#bdbdbd"
  )

  df <- do.call(rbind, lapply(names(filter_stats_all), function(ct) {
    s <- filter_stats_all[[ct]]
    data.frame(
      cell_type = ct,
      category  = factor(EXCL_LEVELS, levels = EXCL_LEVELS),
      n         = c(s$n_passed,
                    s$n_failed_tx_count,
                    s$n_failed_gene_count,
                    s$n_single_isoform),
      stringsAsFactors = FALSE
    )
  }))
  df$cell_type <- factor(df$cell_type, levels = names(filter_stats_all))

  # Stack order: Passed at bottom, exclusions stacked above
  df$category <- factor(df$category,
    levels = rev(EXCL_LEVELS))  # reversed so Passed draws first (bottom)

  # Percent-retained label at top of each bar
  totals <- sapply(filter_stats_all, `[[`, "n_start")
  passed <- sapply(filter_stats_all, `[[`, "n_passed")
  pct_df <- data.frame(
    cell_type = factor(names(totals), levels = names(filter_stats_all)),
    total     = totals,
    pct       = sprintf("%.1f%%\nretained", 100 * passed / totals)
  )

  p <- ggplot(df, aes(x = cell_type, y = n, fill = category)) +
    geom_col(width = 0.65, color = "white", linewidth = 0.3) +
    geom_text(data = pct_df, aes(x = cell_type, y = total, label = pct),
              inherit.aes = FALSE,
              vjust = -0.25, size = 2.8, fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(
      values = EXCL_COLORS,
      breaks = EXCL_LEVELS,   # legend in natural top-to-bottom order
      labels = c(
        "Passed"            = "Passed (kept)",
        "Low isoform count" = sprintf("Excluded: isoform count < %d in ≥ %.0f%% donors",
                                       MIN_TX_COUNT, MIN_SAMPS_FRAC * 100),
        "Low gene total"    = sprintf("Excluded: gene total < %d in ≥ %.0f%% donors",
                                       MIN_GENE_COUNT, MIN_SAMPS_FRAC * 100),
        "Single isoform"    = "Excluded: only 1 isoform left (DTU untestable)"
      ),
      name = NULL
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = "Step 11 — Why transcripts were excluded from DTU testing",
      subtitle = sprintf(
        "Count filter: isoform ≥ %d counts AND gene total ≥ %d counts, both in ≥ %.0f%% of donors",
        MIN_TX_COUNT, MIN_GENE_COUNT, MIN_SAMPS_FRAC * 100),
      x = NULL, y = "Transcripts (n)"
    ) +
    THEME_VIS +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      legend.position = "bottom",
      legend.key.size = unit(0.45, "cm")
    ) +
    guides(fill = guide_legend(nrow = 2))

  save_plot(p, "step11_filter_retention", width = 10, height = 6)
}

# ============================================================
# Step 12 — P-value distributions (refined)
# Y-axis: proportion (comparable across cell types)
# Reference line: uniform expectation
# Coloring: p < 0.05 bins highlighted
# p-type baked into strip label
# ============================================================
plot_step12_pval_hist <- function(results_all_list) {
  N_BINS    <- 20
  BIN_WIDTH <- 1 / N_BINS
  BIN_MIDS  <- seq(BIN_WIDTH / 2, 1 - BIN_WIDTH / 2, by = BIN_WIDTH)

  # Build per-panel binned proportion data
  bin_list <- lapply(names(results_all_list), function(ct) {
    r <- results_all_list[[ct]]
    if (is.null(r)) return(NULL)

    use_empirical <- "empirical_pval" %in% colnames(r) &&
                     sum(!is.na(r$empirical_pval)) > 10
    p_vec <- if (use_empirical) r$empirical_pval else r$pval
    p_vec <- p_vec[!is.na(p_vec)]
    if (length(p_vec) == 0) return(NULL)

    ptype  <- if (use_empirical) "empirical p" else "regular p"
    label  <- paste0(gsub("_", " ", ct), "  [", ptype, ",  n=", length(p_vec), "]")

    counts <- tabulate(findInterval(p_vec,
                                    seq(0, 1, by = BIN_WIDTH),
                                    rightmost.closed = TRUE),
                       nbins = N_BINS)
    data.frame(
      panel_label = label,
      bin_mid     = BIN_MIDS,
      proportion  = counts / length(p_vec),
      stringsAsFactors = FALSE
    )
  })
  bin_df <- do.call(rbind, Filter(Negate(is.null), bin_list))

  # Preserve cell-type order in panels
  panel_order <- unlist(lapply(names(results_all_list), function(ct) {
    r <- results_all_list[[ct]]
    if (is.null(r)) return(NULL)
    use_emp <- "empirical_pval" %in% colnames(r) && sum(!is.na(r$empirical_pval)) > 10
    p_vec   <- if (use_emp) r$empirical_pval else r$pval
    p_vec   <- p_vec[!is.na(p_vec)]
    ptype   <- if (use_emp) "empirical p" else "regular p"
    paste0(gsub("_", " ", ct), "  [", ptype, ",  n=", length(p_vec), "]")
  }))
  bin_df$panel_label <- factor(bin_df$panel_label, levels = panel_order)
  bin_df$sig_bin     <- bin_df$bin_mid < 0.05   # highlight the first bin

  uniform_ref <- BIN_WIDTH   # expected proportion per bin under null

  p <- ggplot(bin_df, aes(x = bin_mid, y = proportion, fill = sig_bin)) +
    geom_col(width = BIN_WIDTH * 0.92, color = "white", linewidth = 0.15) +
    geom_hline(yintercept = uniform_ref,
               linetype = "dashed", color = "#333333", linewidth = 0.55) +
    annotate("text", x = 0.98, y = uniform_ref,
             label = "uniform", vjust = -0.45, hjust = 1,
             size = 2.6, color = "#333333") +
    scale_fill_manual(
      values = c("TRUE" = "#d62728", "FALSE" = "#9ecae1"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
      name   = NULL
    ) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", ".25", ".50", ".75", "1"),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.18))) +
    facet_wrap(~ panel_label, ncol = 3) +
    labs(
      title    = "Step 12 — P-value distributions (satuRn)",
      subtitle = "Proportion of transcripts per bin  |  dashed = uniform null expectation  |  red = p < 0.05 bin",
      x = "p-value", y = "Proportion of transcripts"
    ) +
    THEME_VIS +
    theme(
      legend.position  = "top",
      strip.text       = element_text(face = "bold", size = 9)
    )

  n_ct <- length(unique(bin_df$panel_label))
  save_plot(p, "step12_pval_hist", width = 11, height = ceiling(n_ct / 3) * 3.4 + 2)
}

# ============================================================
# Step 13 — Volcano plots (one panel per cell type)
# ============================================================
plot_step13_volcano <- function(results_all_list) {
  plot_list <- lapply(names(results_all_list), function(ct) {
    res <- results_all_list[[ct]]
    if (is.null(res) || nrow(res) == 0) return(NULL)

    p_col <- if ("empirical_pval" %in% colnames(res) &&
                 sum(!is.na(res$empirical_pval)) > 10) "empirical_pval" else "pval"
    res$plot_p <- pmax(res[[p_col]], 1e-10)
    res <- res[!is.na(res$plot_p) & !is.na(res$delta_psi), ]
    res$neglog10p  <- -log10(res$plot_p)
    res$significant <- !is.na(res$padj_gene) & res$padj_gene < ALPHA &
                       !is.na(res$padj_tx)   & res$padj_tx   < ALPHA
    sig_df <- res[res$significant, ]

    ggplot(res, aes(x = delta_psi, y = neglog10p)) +
      geom_point(data = res[!res$significant, ],
                 color = "#aaaaaa", alpha = 0.4, size = 0.8) +
      geom_point(data = sig_df,
                 color = "#d62728", alpha = 0.9, size = 2) +
      geom_label_repel(data = sig_df,
                       aes(label = transcript_id),
                       size = 2.5, max.overlaps = 20,
                       box.padding = 0.3, label.padding = 0.15,
                       label.size = 0.2, fontface = "bold",
                       fill = "white", color = "#d62728") +
      geom_hline(yintercept = -log10(ALPHA),
                 linetype = "dashed", color = "#666666", linewidth = 0.5) +
      geom_vline(xintercept = 0,
                 linetype = "solid", color = "#dddddd", linewidth = 0.4) +
      annotate("text", x = Inf, y = -log10(ALPHA), vjust = -0.4, hjust = 1.1,
               label = paste0("p=", ALPHA), size = 2.8, color = "#666666") +
      labs(title = gsub("_", " ", ct),
           subtitle = sprintf("%d sig. transcripts | p-type: %s", nrow(sig_df), p_col),
           x = expression(paste(Delta, "PSI (AD − Control)")),
           y = bquote(-log[10](p))) +
      THEME_VIS
  })
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0) return(invisible(NULL))

  chunk <- 4
  n_pages <- ceiling(length(plot_list) / chunk)
  for (i in seq_len(n_pages)) {
    idx <- ((i - 1) * chunk + 1):min(i * chunk, length(plot_list))
    combined <- wrap_plots(plot_list[idx], ncol = 2)
    n_rows <- ceiling(length(idx) / 2)
    save_plot(combined, sprintf("step13_volcano_page%d", i),
              width = 12, height = 5 * n_rows + 0.5)
  }
}

# ============================================================
# Step 13 — Summary: sig results per cell type
# ============================================================
plot_step13_summary <- function(results_all_list) {
  rows <- lapply(names(results_all_list), function(ct) {
    res <- results_all_list[[ct]]
    if (is.null(res)) return(NULL)
    sig <- !is.na(res$padj_gene) & res$padj_gene < ALPHA &
           !is.na(res$padj_tx)   & res$padj_tx   < ALPHA
    n_robust <- if ("robust_to_braak" %in% colnames(res))
                  sum(res$robust_to_braak & sig, na.rm = TRUE) else 0L
    data.frame(cell_type  = ct,
               n_sig_tx   = sum(sig, na.rm = TRUE),
               n_sig_genes = length(unique(res$gene_id[sig])),
               n_robust   = n_robust,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  df$cell_type <- factor(df$cell_type, levels = names(results_all_list))

  # Panel A: genes + transcripts
  df_long <- pivot_longer(df, cols = c("n_sig_genes", "n_sig_tx"),
                          names_to = "type", values_to = "n")
  df_long$type <- factor(df_long$type,
    levels = c("n_sig_genes", "n_sig_tx"),
    labels = c("Genes", "Transcripts"))

  pA <- ggplot(df_long, aes(x = cell_type, y = n, fill = type)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3, fontface = "bold") +
    scale_fill_manual(values = c("Genes" = "#6baed6", "Transcripts" = "#2171b5"), name = NULL) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                       expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Step 13 — Significant DTU results (stageR OFDR < 5%)",
         x = NULL, y = "Count") +
    THEME_VIS +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # Panel B: robust to braak
  pB <- ggplot(df, aes(x = cell_type, y = n_robust,
                        fill = factor(ifelse(n_robust > 0, "robust", "not_robust")))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = ifelse(n_sig_tx > 0,
                                  sprintf("%d/%d", n_robust, n_sig_tx), "")),
              vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("robust" = "#74c476", "not_robust" = "#c7e9c0"),
                      guide = "none") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Robustness: also significant in braak sensitivity model",
         subtitle = "robust / total significant per cell type",
         x = NULL, y = "Robust transcripts (n)") +
    THEME_VIS +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  combined <- pA / pB + plot_layout(heights = c(2, 1.2))
  save_plot(combined, "step13_summary", width = 9, height = 9)
}

# ============================================================
# Step 14 — Dumbbell plot: Control → AD segment, Active Control diamond
# Each row = one significant transcript, ordered by ΔPSI.
# % of switch = (PSI_active − PSI_ctrl) / (PSI_AD − PSI_ctrl)
# ============================================================
plot_step14_active_psi <- function(sig_results) {
  if (!"psi_active_ctrl" %in% colnames(sig_results)) return(invisible(NULL))
  has_active <- !is.na(sig_results$psi_active_ctrl)
  if (sum(has_active) == 0) {
    message("  No Active control PSI data — skipping step14 plot")
    return(invisible(NULL))
  }

  df <- sig_results[has_active, c("transcript_id", "gene_id", "cell_type",
                                   "psi_AD", "psi_ctrl", "psi_active_ctrl",
                                   "delta_psi", "padj_gene")]
  df <- df[order(df$delta_psi), ]

  # Fraction of the Control→AD shift reproduced by Active Control
  # 0% = Active looks like Control; 100% = Active looks like AD
  df$switch_pct <- (df$psi_active_ctrl - df$psi_ctrl) /
                   (df$psi_AD          - df$psi_ctrl) * 100

  # Row label: transcript + cell type + ΔPSI
  df$row_label <- factor(
    paste0(df$transcript_id, "   (",
           gsub("_", " ", df$cell_type), ")"),
    levels = paste0(df$transcript_id, "   (",
                    gsub("_", " ", df$cell_type), ")")
  )

  # Segment midpoint x for placing the % label
  df$seg_mid <- (df$psi_ctrl + df$psi_AD) / 2

  # Alternating row background
  n_tx   <- nrow(df)
  bg_df  <- data.frame(
    ymin  = seq_len(n_tx) - 0.5,
    ymax  = seq_len(n_tx) + 0.5,
    shade = seq_len(n_tx) %% 2 == 0
  )
  bg_df$row_num <- seq_len(n_tx)

  p <- ggplot(df, aes(y = row_label)) +
    # Alternating row shading
    geom_rect(data = bg_df[bg_df$shade, ],
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, fill = "#f5f5f5", color = NA) +
    # Segment: Control → AD
    geom_segment(aes(x    = psi_ctrl, xend = psi_AD,
                     y    = row_label, yend = row_label,
                     color = delta_psi > 0),
                 linewidth = 2.8, lineend = "round") +
    # Control endpoint
    geom_point(aes(x = psi_ctrl), shape = 21, size = 5,
               fill = COND_COLORS["Control"], color = "white", stroke = 1) +
    # AD endpoint
    geom_point(aes(x = psi_AD), shape = 21, size = 5,
               fill = COND_COLORS["AD"], color = "white", stroke = 1) +
    # Active Control — diamond, sits on the segment
    geom_point(aes(x = psi_active_ctrl), shape = 23, size = 5,
               fill = COND_COLORS["Active control"], color = "white", stroke = 1) +
    # % of switch label above each segment
    geom_text(aes(x = seg_mid,
                  label = sprintf("%+.0f%%", switch_pct)),
              vjust = -0.7, size = 3, fontface = "bold", color = "#333333") +
    # Vertical reference at PSI = 0.5
    geom_vline(xintercept = 0.5, linetype = "dotted",
               color = "#bbbbbb", linewidth = 0.5) +
    scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "#2171b5"),
                       guide = "none") +
    scale_x_continuous(limits  = c(0, 1),
                       breaks  = c(0, 0.25, 0.5, 0.75, 1),
                       labels  = c("0", ".25", ".50", ".75", "1"),
                       expand  = expansion(mult = c(0.03, 0.03))) +
    # Manual legend
    annotate("point", x = 0.72, y = n_tx + 0.8, shape = 21, size = 4,
             fill = COND_COLORS["Control"], color = "white", stroke = 1) +
    annotate("text",  x = 0.74, y = n_tx + 0.8, label = "Control",
             hjust = 0, size = 3.2) +
    annotate("point", x = 0.84, y = n_tx + 0.8, shape = 21, size = 4,
             fill = COND_COLORS["AD"], color = "white", stroke = 1) +
    annotate("text",  x = 0.86, y = n_tx + 0.8, label = "AD",
             hjust = 0, size = 3.2) +
    annotate("point", x = 0.92, y = n_tx + 0.8, shape = 23, size = 4,
             fill = COND_COLORS["Active control"], color = "white", stroke = 1) +
    annotate("text",  x = 0.94, y = n_tx + 0.8, label = "Active ctrl",
             hjust = 0, size = 3.2) +
    labs(
      title    = "Step 14 — Active Control PSI relative to AD vs Control switch",
      subtitle = "Segment = Control (●) → AD (●) isoform shift  |  Diamond (◆) = Active control mean PSI  |  % = fraction of AD shift reproduced",
      x = "Mean PSI", y = NULL
    ) +
    THEME_VIS +
    theme(
      legend.position   = "none",
      panel.grid.major.x = element_line(color = "#eeeeee", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      axis.text.y       = element_text(size = 9, face = "bold")
    ) +
    coord_cartesian(clip = "off")

  save_plot(p, "step14_active_psi_comparison",
            width = 9, height = max(4, n_tx * 0.55 + 2.5))
}

# ============================================================
# Viz — Per-donor PSI strip plots for each significant transcript
# ============================================================
compute_donor_psi_df <- function(counts_mat, meta, tx_subset = NULL) {
  if (!is.null(tx_subset))
    counts_mat <- counts_mat[rownames(counts_mat) %in% tx_subset, , drop = FALSE]

  tx_names   <- rownames(counts_mat)
  gene_names <- sub("-[^-]+$", "", tx_names)
  gene_totals <- rowsum(counts_mat, gene_names)

  do.call(rbind, lapply(colnames(counts_mat), function(donor) {
    ct_vec <- counts_mat[, donor]
    gt_vec <- gene_totals[gene_names, donor]
    psi    <- ifelse(gt_vec > 0, ct_vec / gt_vec, NA_real_)
    data.frame(transcript_id = tx_names, gene_id = gene_names,
               donor = donor, condition = meta[donor, "condition"],
               psi = psi, stringsAsFactors = FALSE)
  }))
}

plot_viz_psi_donors <- function(sig_results) {
  if (nrow(sig_results) == 0) return(invisible(NULL))

  for (ct in unique(sig_results$cell_type)) {
    sig_ct <- sig_results[sig_results$cell_type == ct, ]

    counts_raw <- tryCatch(
      read.csv(file.path(IN_DIR, paste0("counts_", ct, ".csv")),
               row.names = 1, check.names = FALSE),
      error = function(e) NULL)
    if (is.null(counts_raw)) next

    meta_raw   <- read.csv(file.path(IN_DIR, paste0("metadata_", ct, ".csv")),
                            row.names = 1, check.names = FALSE)
    counts_mat <- t(as.matrix(counts_raw))

    donor_psi <- compute_donor_psi_df(counts_mat, meta_raw,
                                       tx_subset = sig_ct$transcript_id)
    donor_psi <- donor_psi[donor_psi$condition %in% c("AD", "Control"), ]
    if (nrow(donor_psi) == 0) next

    donor_psi$condition <- factor(donor_psi$condition, levels = c("Control", "AD"))

    # Add ΔPSI facet subtitle
    dpsi_map <- setNames(round(sig_ct$delta_psi, 3), sig_ct$transcript_id)
    donor_psi$tx_label <- paste0(donor_psi$transcript_id, "\nΔPSI = ",
                                   sprintf("%+.3f", dpsi_map[donor_psi$transcript_id]))

    n_tx  <- nrow(sig_ct)
    n_col <- min(4, n_tx)

    p <- ggplot(donor_psi, aes(x = condition, y = psi, color = condition)) +
      geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
      stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "crossbar",
                   width = 0.45, linewidth = 0.7, color = "#2C3E50",
                   middle.linewidth = 1.5) +
      scale_color_manual(values = COND_COLORS, name = NULL) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
      facet_wrap(~ tx_label, ncol = n_col) +
      labs(title = paste0(gsub("_", " ", ct), " — Per-donor isoform PSI"),
           subtitle = "Points = individual donors  |  crossbar = group mean",
           x = NULL, y = "PSI") +
      THEME_VIS +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "none")

    save_plot(p, paste0("viz_psi_donors_", ct),
              width  = max(6, n_col * 3.5),
              height = max(4, ceiling(n_tx / n_col) * 4))
  }
}

# ============================================================
# Viz — ΔPSI heatmap across all significant transcripts × cell types
# ============================================================
plot_viz_heatmap_dpsi <- function(results_all_list) {
  rows <- lapply(names(results_all_list), function(ct) {
    res <- results_all_list[[ct]]
    if (is.null(res)) return(NULL)
    sig <- !is.na(res$padj_gene) & res$padj_gene < ALPHA &
           !is.na(res$padj_tx)   & res$padj_tx   < ALPHA
    if (sum(sig) == 0) return(NULL)
    data.frame(transcript_id = res$transcript_id[sig],
               cell_type     = ct,
               delta_psi     = res$delta_psi[sig],
               stringsAsFactors = FALSE)
  })
  sig_df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(sig_df) || nrow(sig_df) == 0) {
    message("  No significant results — skipping heatmap")
    return(invisible(NULL))
  }

  mat <- pivot_wider(sig_df, id_cols = transcript_id, names_from = cell_type,
                     values_from = delta_psi, values_fill = NA)
  rnames  <- mat$transcript_id
  mat     <- as.matrix(mat[, -1])
  rownames(mat) <- rnames

  n_tx <- nrow(mat)
  n_ct <- ncol(mat)
  if (n_tx < 2) return(invisible(NULL))

  limit  <- max(abs(mat), na.rm = TRUE)
  breaks <- seq(-limit, limit, length.out = 101)
  colors <- colorRampPalette(c("#2171b5", "#f7f7f7", "#d62728"))(100)

  # pheatmap cannot cluster when NAs are present — use 0-imputed matrix for clustering
  mat_for_clust <- mat
  mat_for_clust[is.na(mat_for_clust)] <- 0
  row_clust <- if (n_tx > 1) hclust(dist(mat_for_clust)) else FALSE
  col_clust <- if (n_ct > 1) hclust(dist(t(mat_for_clust))) else FALSE

  dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(OUT_DIR, "plots", "viz_heatmap_dpsi.png")
  png(out_path,
      width  = max(5, n_ct * 1.5 + 3),
      height = max(4, n_tx * 0.4 + 2.5),
      units  = "in", res = 300, bg = "white")
  pheatmap(mat,
           color        = colors,
           breaks       = breaks,
           na_col       = "#eeeeee",
           cluster_rows = row_clust,
           cluster_cols = col_clust,
           border_color = NA,
           cellwidth    = 45,
           fontsize_row = 9,
           fontsize_col = 10,
           main         = "ΔPSI heatmap — significant transcripts across cell types",
           angle_col    = 45,
           legend_breaks = c(-limit, -limit/2, 0, limit/2, limit),
           legend_labels = c(sprintf("%.2f", -limit), "", "0",
                             "", sprintf("%.2f", limit)))
  dev.off()
  message("  Plot: viz_heatmap_dpsi.png")
}

# ============================================================
# Viz — Isoform slope charts for DTU genes
# One panel per gene; Control → AD slope per isoform.
# Significant isoforms: thick solid line + filled dot + ΔPSI label.
# Non-significant isoforms: thin dashed line + open dot (context only).
# ============================================================
plot_viz_isoform_props <- function(sig_results, max_genes_per_ct = 12) {
  if (nrow(sig_results) == 0) return(invisible(NULL))

  # Bold, saturated, colorblind-distinguishable palette (not pastel)
  ISO_PALETTE <- c(
    "#1b7837", "#762a83", "#e08214", "#4575b4",
    "#d73027", "#01665e", "#8073ac", "#bf812d",
    "#de77ae", "#4d4d4d", "#35978f", "#c51b7d"
  )

  for (ct in unique(sig_results$cell_type)) {
    sig_ct    <- sig_results[sig_results$cell_type == ct, ]
    top_genes <- head(unique(sig_ct$gene_id), max_genes_per_ct)

    counts_raw <- tryCatch(
      read.csv(file.path(IN_DIR, paste0("counts_", ct, ".csv")),
               row.names = 1, check.names = FALSE),
      error = function(e) NULL)
    if (is.null(counts_raw)) next

    meta_raw   <- read.csv(file.path(IN_DIR, paste0("metadata_", ct, ".csv")),
                            row.names = 1, check.names = FALSE)
    counts_mat <- t(as.matrix(counts_raw))

    # One ggplot per gene
    gene_plots <- lapply(top_genes, function(gn) {
      tx_in <- rownames(counts_mat)[grepl(paste0("^", gn, "-"), rownames(counts_mat))]
      if (length(tx_in) < 2) return(NULL)

      gene_tot <- colSums(counts_mat[tx_in, , drop = FALSE])
      psi_long <- do.call(rbind, lapply(tx_in, function(tx) {
        psi  <- ifelse(gene_tot > 0, counts_mat[tx, ] / gene_tot, NA_real_)
        cond <- meta_raw[names(psi), "condition"]
        data.frame(transcript = tx, condition = cond, psi = psi,
                   stringsAsFactors = FALSE)
      }))
      psi_long <- psi_long[psi_long$condition %in% c("Control", "AD"), ]
      psi_long$condition <- factor(psi_long$condition, levels = c("Control", "AD"))

      # Mean PSI per transcript × condition
      summ <- aggregate(psi ~ transcript + condition, data = psi_long,
                        FUN = function(x) mean(x, na.rm = TRUE))

      # Sort isoforms by Control PSI descending for consistent color assignment
      ctrl_ord  <- summ$transcript[summ$condition == "Control"]
      ctrl_psi  <- summ$psi[summ$condition == "Control"]
      iso_levels <- ctrl_ord[order(-ctrl_psi)]
      summ$transcript <- factor(summ$transcript, levels = iso_levels)

      # Assign one bold color per isoform
      n_iso      <- nlevels(summ$transcript)
      iso_colors <- setNames(ISO_PALETTE[seq_len(n_iso)], levels(summ$transcript))

      # Mark significance and build short label
      summ$is_sig   <- as.character(summ$transcript) %in% sig_ct$transcript_id
      summ$tx_short <- sub(paste0("^", gn, "-"), "", as.character(summ$transcript))

      # ΔPSI for each isoform
      psi_wide <- reshape(summ[, c("transcript", "condition", "psi")],
                          idvar = "transcript", timevar = "condition", direction = "wide")
      names(psi_wide) <- gsub("psi\\.", "", names(psi_wide))
      psi_wide$dpsi <- psi_wide$AD - psi_wide$Control

      # Label data: AD endpoint for significant isoforms only
      label_df <- merge(
        summ[summ$condition == "AD" & summ$is_sig, ],
        psi_wide[, c("transcript", "dpsi")],
        by = "transcript"
      )
      label_df$label_text <- paste0(
        label_df$tx_short, "\n",
        sprintf("(%+.2f)", label_df$dpsi)
      )

      # Split for layered drawing
      sig_data   <- summ[summ$is_sig, ]
      other_data <- summ[!summ$is_sig, ]

      n_sig_in_gene <- sum(unique(as.character(summ$transcript)) %in%
                           sig_ct$transcript_id)
      panel_title <- sprintf("%s  ·  %d sig. isoform%s",
                             gn, n_sig_in_gene,
                             if (n_sig_in_gene != 1) "s" else "")

      ggplot(summ, aes(x = condition, y = psi,
                        group = transcript, color = transcript)) +
        # Non-significant: thin solid, faded — individual color kept
        (if (nrow(other_data) > 0)
          geom_line(data = other_data, linewidth = 0.5, alpha = 0.30) else NULL) +
        (if (nrow(other_data) > 0)
          geom_point(data = other_data, size = 2.0, shape = 1,
                     alpha = 0.38, stroke = 0.9) else NULL) +
        # Significant: white halo then thick colored line
        (if (nrow(sig_data) > 0)
          geom_line(data = sig_data, color = "white", linewidth = 4.2) else NULL) +
        (if (nrow(sig_data) > 0)
          geom_line(data = sig_data, linewidth = 2.4) else NULL) +
        (if (nrow(sig_data) > 0)
          geom_point(data = sig_data, color = "white", size = 7.5, shape = 16) else NULL) +
        (if (nrow(sig_data) > 0)
          geom_point(data = sig_data, size = 5, shape = 16) else NULL) +
        # ΔPSI labels at AD endpoint
        (if (nrow(label_df) > 0)
          geom_text(data = label_df,
                    aes(x = condition, y = psi, label = label_text),
                    hjust = -0.12, size = 2.9, fontface = "bold",
                    lineheight = 0.85, show.legend = FALSE) else NULL) +
        scale_color_manual(values = iso_colors, guide = "none") +
        scale_y_continuous(limits = c(0, 1),
                           labels = percent_format(accuracy = 1),
                           expand = expansion(mult = c(0.06, 0.06))) +
        scale_x_discrete(expand = expansion(add = c(0.25, 1.1))) +
        labs(title = panel_title, x = NULL, y = "Mean PSI") +
        THEME_VIS +
        theme(
          plot.title         = element_text(face = "bold", size = 10),
          axis.title.y       = element_text(size = 8),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "#e8e8e8", linewidth = 0.4)
        ) +
        coord_cartesian(clip = "off")
    })

    gene_plots <- Filter(Negate(is.null), gene_plots)
    if (length(gene_plots) == 0) next

    n_plots <- length(gene_plots)
    # Pick ncol (2–4) that leaves the fewest empty slots; break ties toward more columns
    cands  <- seq(2L, min(4L, n_plots))
    empties <- n_plots %% cands
    n_cols  <- max(cands[empties == min(empties)])
    n_rows  <- ceiling(n_plots / n_cols)

    combined <- wrap_plots(gene_plots, ncol = n_cols) +
      plot_annotation(
        title    = paste0(gsub("_", " ", ct), " — Isoform usage (DTU genes)"),
        subtitle = "Each line = one isoform  |  thick = significant  |  label = ΔPSI (AD − Control)",
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "#555555")
        )
      )

    save_plot(combined, paste0("viz_isoform_props_", ct),
              width  = max(7, n_cols * 3.8),
              height = max(4, n_rows * 4.5 + 1.2))
  }
}

# ============================================================
# Viz — ΔPSI comparison primary vs braak sensitivity model
# ============================================================
plot_viz_braak_comparison <- function(results_all_list) {
  rows <- lapply(names(results_all_list), function(ct) {
    res <- results_all_list[[ct]]
    if (is.null(res) || !"robust_to_braak" %in% colnames(res)) return(NULL)
    sig <- !is.na(res$padj_gene) & res$padj_gene < ALPHA &
           !is.na(res$padj_tx)   & res$padj_tx   < ALPHA
    if (sum(sig) == 0) return(NULL)
    sig_res <- res[sig, ]
    sig_res$cell_type <- ct
    sig_res$robust_label <- ifelse(sig_res$robust_to_braak,
                                    "Robust to braak", "Primary only")
    sig_res
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  df$tx_label <- paste0(df$transcript_id, "\n(", gsub("_", " ", df$cell_type), ")")

  p <- ggplot(df, aes(x = reorder(tx_label, delta_psi), y = delta_psi,
                       fill = robust_label)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c("Robust to braak" = "#2ca02c", "Primary only" = "#d62728"),
                      name = NULL) +
    scale_y_continuous(labels = function(x) sprintf("%+.2f", x)) +
    labs(title = "Significant DTU transcripts — ΔPSI and braak robustness",
         subtitle = "Green = significant in both primary and braak sensitivity model",
         x = NULL, y = "ΔPSI (AD − Control)") +
    THEME_VIS +
    theme(legend.position = "top")

  n_tx <- nrow(df)
  save_plot(p, "viz_sig_dpsi_robustness",
            width = 8, height = max(4, n_tx * 0.45 + 2))
}

# ============================================================
# Master function: generate all plots
# ============================================================
# ============================================================
# Viz — Per-isoform deep-dive: one image per significant transcript
# Left panel : stacked bar of mean PSI per condition (+ Active Control if available)
# Right panel: gene slope chart with focal isoform emphasised
# ============================================================
plot_viz_per_isoform <- function(sig_results) {
  if (nrow(sig_results) == 0) return(invisible(NULL))

  ISO_PALETTE <- c(
    "#1b7837", "#762a83", "#e08214", "#4575b4",
    "#d73027", "#01665e", "#8073ac", "#bf812d",
    "#de77ae", "#4d4d4d", "#35978f", "#c51b7d"
  )

  for (ct in unique(sig_results$cell_type)) {
    sig_ct <- sig_results[sig_results$cell_type == ct, ]

    counts_raw <- tryCatch(
      read.csv(file.path(IN_DIR, paste0("counts_", ct, ".csv")),
               row.names = 1, check.names = FALSE),
      error = function(e) NULL)
    if (is.null(counts_raw)) next
    meta_raw   <- read.csv(file.path(IN_DIR, paste0("metadata_", ct, ".csv")),
                            row.names = 1, check.names = FALSE)
    counts_mat <- t(as.matrix(counts_raw))

    # Active control (optional)
    counts_raw_ac <- tryCatch(
      read.csv(file.path(IN_DIR, paste0("counts_", ct, "_sensitivity.csv")),
               row.names = 1, check.names = FALSE),
      error = function(e) NULL)
    meta_raw_ac   <- if (!is.null(counts_raw_ac)) tryCatch(
      read.csv(file.path(IN_DIR, paste0("metadata_", ct, "_sensitivity.csv")),
               row.names = 1, check.names = FALSE),
      error = function(e) NULL) else NULL
    counts_mat_ac <- if (!is.null(counts_raw_ac)) t(as.matrix(counts_raw_ac)) else NULL

    for (i in seq_len(nrow(sig_ct))) {
      tx_id <- sig_ct$transcript_id[i]
      gn    <- sig_ct$gene_id[i]
      if (!tx_id %in% rownames(counts_mat)) next

      tx_short  <- sub(paste0("^", gn, "-"), "", tx_id)
      dpsi_val  <- sig_ct$delta_psi[i]
      padj_gene <- sig_ct$padj_gene[i]
      padj_tx   <- sig_ct$padj_tx[i]

      # ── Shared: PSI for all isoforms in this gene ─────────────────
      tx_in <- rownames(counts_mat)[grepl(paste0("^", gn, "-"), rownames(counts_mat))]
      if (length(tx_in) < 2) next

      gene_tot <- colSums(counts_mat[tx_in, , drop = FALSE])
      psi_long <- do.call(rbind, lapply(tx_in, function(tx) {
        psi    <- ifelse(gene_tot > 0, counts_mat[tx, ] / gene_tot, NA_real_)
        cond_v <- meta_raw[colnames(counts_mat), "condition"]
        data.frame(transcript = tx, condition = cond_v, psi = psi,
                   stringsAsFactors = FALSE)
      }))
      psi_long <- psi_long[psi_long$condition %in% c("Control", "AD"), ]
      psi_long$condition <- factor(psi_long$condition, levels = c("Control", "AD"))

      # Mean PSI for slope chart (Control + AD)
      summ <- aggregate(psi ~ transcript + condition, data = psi_long,
                        FUN = function(x) mean(x, na.rm = TRUE))

      # Color assignment: focal isoform first → gets color [1]
      iso_levels <- c(tx_id, setdiff(unique(summ$transcript), tx_id))
      summ$transcript <- factor(summ$transcript, levels = iso_levels)
      n_iso      <- length(iso_levels)
      iso_colors <- setNames(ISO_PALETTE[seq_len(n_iso)], iso_levels)

      # ── Left panel: stacked bar (Control, AD, + Active Control) ───
      bar_df <- aggregate(psi ~ transcript + condition, data = psi_long,
                          FUN = function(x) mean(x, na.rm = TRUE))

      # Add Active Control mean PSI if available
      if (!is.null(counts_mat_ac)) {
        tx_in_ac <- intersect(tx_in, rownames(counts_mat_ac))
        if (length(tx_in_ac) > 0) {
          gene_tot_ac <- colSums(counts_mat_ac[tx_in_ac, , drop = FALSE])
          psi_ac <- do.call(rbind, lapply(tx_in_ac, function(tx) {
            psi    <- ifelse(gene_tot_ac > 0, counts_mat_ac[tx, ] / gene_tot_ac, NA_real_)
            cond_v <- meta_raw_ac[colnames(counts_mat_ac), "condition"]
            data.frame(transcript = tx, condition = cond_v, psi = psi,
                       stringsAsFactors = FALSE)
          }))
          psi_ac <- psi_ac[psi_ac$condition == "Active control", ]
          if (nrow(psi_ac) > 0) {
            bar_ac <- aggregate(psi ~ transcript + condition, data = psi_ac,
                                FUN = function(x) mean(x, na.rm = TRUE))
            bar_df <- rbind(bar_df, bar_ac)
          }
        }
      }

      cond_levs <- c("Control", "AD")
      if ("Active control" %in% bar_df$condition) cond_levs <- c(cond_levs, "Active control")
      bar_df$transcript <- factor(bar_df$transcript, levels = iso_levels)
      bar_df$condition  <- factor(bar_df$condition, levels = cond_levs)
      bar_df$is_focal   <- as.character(bar_df$transcript) == tx_id

      p_bar <- ggplot(bar_df, aes(x = condition, y = psi,
                                   fill = transcript, color = is_focal)) +
        geom_col(position = "stack", linewidth = 0.8, width = 0.65) +
        scale_fill_manual(values = iso_colors, guide = "none") +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"),
                           guide = "none") +
        scale_y_continuous(labels = percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.02)),
                           limits = c(0, 1.01)) +
        labs(title = "Stacked proportions",
             subtitle = paste0("black border = ", tx_short),
             x = NULL, y = "Mean PSI") +
        THEME_VIS +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "#e8e8e8", linewidth = 0.4))

      # ── Right panel: gene slope chart, focal isoform prominent ───

      summ$is_focal <- as.character(summ$transcript) == tx_id
      summ$is_sig   <- as.character(summ$transcript) %in% sig_ct$transcript_id &
                       !summ$is_focal   # other sig isoforms (not focal)
      summ$tx_lbl   <- sub(paste0("^", gn, "-"), "", as.character(summ$transcript))

      # ΔPSI labels only for focal isoform at AD endpoint
      psi_wide <- reshape(summ[, c("transcript", "condition", "psi")],
                          idvar = "transcript", timevar = "condition", direction = "wide")
      names(psi_wide) <- sub("^psi\\.", "", names(psi_wide))
      psi_wide$dpsi <- psi_wide$AD - psi_wide$Control

      focal_label <- merge(
        summ[summ$condition == "AD" & summ$is_focal, ],
        psi_wide[, c("transcript", "dpsi")],
        by = "transcript"
      )
      focal_label$label_text <- paste0(focal_label$tx_lbl, "\n",
                                        sprintf("(%+.2f)", focal_label$dpsi))

      focal_data <- summ[summ$is_focal, ]
      sig_data   <- summ[summ$is_sig, ]
      other_data <- summ[!summ$is_focal & !summ$is_sig, ]

      # ── Legend data: transcript short name + padj_tx ─────────────
      tx_lgd_df <- data.frame(transcript = iso_levels,
                               tx_short   = sub(paste0("^", gn, "-"), "", iso_levels),
                               stringsAsFactors = FALSE)
      gene_sig_df <- sig_ct[sig_ct$gene_id == gn, c("transcript_id", "padj_tx")]
      tx_lgd_df   <- merge(tx_lgd_df, gene_sig_df,
                            by.x = "transcript", by.y = "transcript_id", all.x = TRUE)
      tx_lgd_df   <- tx_lgd_df[match(iso_levels, tx_lgd_df$transcript), ]
      tx_lgd_df$label <- ifelse(
        is.na(tx_lgd_df$padj_tx),
        paste0(tx_lgd_df$tx_short, "  n.s."),
        ifelse(
          tx_lgd_df$padj_tx == 0,
          paste0(tx_lgd_df$tx_short, "  p < 2.2e-16"),
          paste0(tx_lgd_df$tx_short,
                 "  p = ", formatC(tx_lgd_df$padj_tx, format = "e", digits = 1))
        )
      )

      # Line style per isoform (used in both slope chart and legend panel)
      lwd_lgd   <- ifelse(iso_levels == tx_id, 3.2,
                           ifelse(iso_levels %in% sig_ct$transcript_id, 1.6, 0.6))
      ltype_lgd <- rep("solid", length(iso_levels))   # all solid — no dashes
      alpha_lgd <- ifelse(iso_levels == tx_id, 1.0,
                           ifelse(iso_levels %in% sig_ct$transcript_id, 0.75, 0.30))
      is_drawn  <- iso_levels == tx_id | iso_levels %in% sig_ct$transcript_id

      p_slope <- ggplot(summ, aes(x = condition, y = psi,
                                   group = transcript, color = transcript)) +
        # Non-significant: thin solid, faded — individual color kept
        (if (nrow(other_data) > 0)
          geom_line(data = other_data, linewidth = 0.5, alpha = 0.30) else NULL) +
        (if (nrow(other_data) > 0)
          geom_point(data = other_data, size = 1.8, shape = 1,
                     alpha = 0.38, stroke = 0.9) else NULL) +
        # Other significant isoforms: medium weight
        (if (nrow(sig_data) > 0)
          geom_line(data = sig_data, linewidth = 1.8, alpha = 0.78) else NULL) +
        (if (nrow(sig_data) > 0)
          geom_point(data = sig_data, size = 3.5, shape = 16, alpha = 0.85) else NULL) +
        # Focal isoform: white halo first, then colored line on top
        geom_line(data = focal_data, color = "white", linewidth = 5.8) +
        geom_line(data = focal_data, linewidth = 3.2) +
        geom_point(data = focal_data, color = "white", size = 9.5, shape = 16) +
        geom_point(data = focal_data, size = 6.5, shape = 16) +
        # ΔPSI label for focal isoform at AD endpoint
        (if (nrow(focal_label) > 0)
          geom_text(data = focal_label,
                    aes(x = condition, y = psi, label = label_text),
                    hjust = -0.12, size = 3.2, fontface = "bold",
                    lineheight = 0.85, show.legend = FALSE) else NULL) +
        scale_color_manual(values = iso_colors, guide = "none") +
        scale_y_continuous(limits = c(0, 1),
                           labels = percent_format(accuracy = 1),
                           expand = expansion(mult = c(0.06, 0.06))) +
        scale_x_discrete(expand = expansion(add = c(0.25, 0.3))) +
        labs(title = "Isoform proportions",
             subtitle = sprintf("%s  (%d isoforms in gene)", gn, length(tx_in)),
             x = NULL, y = "Mean PSI") +
        THEME_VIS +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "#e8e8e8", linewidth = 0.4)) +
        coord_cartesian(clip = "off")

      # ── Custom legend panel ───────────────────────────────────────
      n_iso  <- length(iso_levels)
      leg_df <- data.frame(
        y      = rev(seq_len(n_iso)),
        col    = iso_colors[iso_levels],
        label  = tx_lgd_df$label,
        lwd    = lwd_lgd,
        ltype  = ltype_lgd,
        alp    = alpha_lgd,
        solid  = is_drawn,
        stringsAsFactors = FALSE
      )

      p_legend <- ggplot(leg_df, aes(y = y)) +
        # Colored square = stacked bar segment
        geom_tile(aes(x = 0.22, fill = I(col)),
                  width = 0.32, height = 0.55, show.legend = FALSE) +
        # Short line = slope chart line style
        geom_segment(aes(x = 0.50, xend = 0.96, yend = y,
                          color = I(col),
                          linewidth = I(pmin(lwd * 0.38, 1.3)),
                          linetype = ltype,
                          alpha = I(alp)),
                     show.legend = FALSE) +
        scale_linetype_identity() +
        # Dot = slope chart point (filled for sig, open for n.s.)
        geom_point(data = leg_df[leg_df$solid, ],
                   aes(x = 0.73, color = I(col), alpha = I(alp)),
                   size = 2.3, shape = 16, show.legend = FALSE) +
        geom_point(data = leg_df[!leg_df$solid, ],
                   aes(x = 0.73, color = I(col), alpha = I(alp)),
                   size = 1.6, shape = 1, stroke = 0.8, show.legend = FALSE) +
        # Isoform label
        geom_text(aes(x = 1.08, label = label),
                  hjust = 0, size = 2.55, family = "mono") +
        scale_x_continuous(limits = c(0.02, 5.5), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0.2, n_iso + 0.8), expand = c(0, 0)) +
        theme_void() +
        theme(plot.margin = margin(4, 2, 4, 2))

      # ── Combine ───────────────────────────────────────────────────
      combined <- (p_bar | p_slope | p_legend) +
        plot_layout(widths = c(1, 1, 0.65)) +
        plot_annotation(
          title    = sprintf("%s  ·  %s  (%s)", gn, tx_short, gsub("_", " ", ct)),
          subtitle = sprintf("ΔPSI = %+.3f  |  padj_gene = %.2e  |  padj_tx = %.2e",
                             dpsi_val, padj_gene, padj_tx),
          theme = theme(
            plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "#444444")
          )
        )

      tx_safe <- gsub("[^A-Za-z0-9_-]", "_", tx_id)
      save_plot(combined,
                paste0("per_isoform_", ct, "_", tx_safe),
                width = 9, height = 4.5)
    }
  }
}

generate_all_plots <- function(results_all_list, filter_stats_all, sig_results) {
  message("\n", strrep("=", 60))
  message("Generating Layer 1 plots...")
  message(strrep("=", 60))

  tryCatch(plot_step11_filter(filter_stats_all),
           error = function(e) message("  step11 plot failed: ", conditionMessage(e)))

  tryCatch(plot_step12_pval_hist(results_all_list),
           error = function(e) message("  step12 plot failed: ", conditionMessage(e)))

  tryCatch(plot_step13_volcano(results_all_list),
           error = function(e) message("  step13 volcano failed: ", conditionMessage(e)))

  tryCatch(plot_step13_summary(results_all_list),
           error = function(e) message("  step13 summary failed: ", conditionMessage(e)))

  if (!is.null(sig_results) && nrow(sig_results) > 0) {
    tryCatch(plot_step14_active_psi(sig_results),
             error = function(e) message("  step14 plot failed: ", conditionMessage(e)))

    tryCatch(plot_viz_psi_donors(sig_results),
             error = function(e) message("  viz psi donors failed: ", conditionMessage(e)))

    tryCatch(plot_viz_isoform_props(sig_results),
             error = function(e) message("  viz isoform props failed: ", conditionMessage(e)))

    tryCatch(plot_viz_per_isoform(sig_results),
             error = function(e) message("  viz per-isoform failed: ", conditionMessage(e)))

    tryCatch(plot_viz_braak_comparison(results_all_list),
             error = function(e) message("  viz braak comparison failed: ", conditionMessage(e)))
  }

  tryCatch(plot_viz_heatmap_dpsi(results_all_list),
           error = function(e) message("  viz heatmap failed: ", conditionMessage(e)))

  message("\nPlots saved to: ", file.path(OUT_DIR, "plots"))
}

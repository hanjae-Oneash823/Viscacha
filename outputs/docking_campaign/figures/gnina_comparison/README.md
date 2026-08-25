# GNINA comparison figures

These figures were rebuilt with R 4.5.2, ggplot2, patchwork, ragg, and svglite. Each is exported as a 3408-pixel-wide PNG for slides and as vector PDF and SVG files for lossless scaling.

- `gnina_matched_comparison.*`: the main findings figure. The left estimation plot shows seed-paired alternate-minus-canonical CNNaffinity shifts with means and 95% t intervals. The right plot shows the RMSD of each CNN-selected pose and the 2 Å pose-valid region.
- `gnina_pose_selection_validation.*`: a single paired canonical-validation chart. Thin lines retain the link between Vina and GNINA selections from each search; large symbols show group means and direct labels report the number of GNINA selections below 2 Å.
- `vina_gnina_rank_agreement.*`: a horizontal agreement dashboard. The left plot reports which original Vina rank GNINA selected; the right plot shows the complete run-level distribution of within-ensemble Spearman correlations.

The source is `02_SURVEYOR/master_surveyor/plot_gnina_comparison.R`. It reads the audited `analysis/gnina/run_summary.csv` table and does not recompute or alter docking results. The design uses direct labeling, restrained color, confidence intervals, and shaded decision/reference regions. Vina energy and CNNaffinity are never placed on a shared numerical axis.

GNINA used `--score_only`; the plotted pose changes are selection changes within the retained Vina ensemble, not coordinate optimization.

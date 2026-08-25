# GNINA comparison figures

Each figure is available as a 360-dpi PNG plus vector PDF and SVG.

- `gnina_pose_selection_validation.*`: paired Vina-versus-GNINA pose selection for the three exact canonical controls. This is the clearest figure for showing BACE1/CHRNA7 preservation and the FYN disagreement.
- `gnina_matched_comparison.*`: primary Vina-rank-1 CNNaffinity rescoring and secondary CNN-selected geometry for BACE1 and CHRNA7/CHRFAM7A.
- `vina_gnina_rank_agreement.*`: original Vina ranks selected by CNNscore and mean within-run rank correlation.

GNINA used `--score_only`; the plotted pose changes are selection changes within the retained Vina ensemble, not coordinate optimization.

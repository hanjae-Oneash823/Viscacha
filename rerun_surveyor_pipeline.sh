#!/usr/bin/env bash
set -euo pipefail

PY=/home/welcome3/anaconda3/envs/oneash_dtu/bin/python
cd /home/welcome3/Viscacha_pipeline

# Cap BLAS/OMP threading to 6 cores across every step (must be set before
# numpy/pandas/anndata are imported, so export at the shell level, not in Python).
export OPENBLAS_NUM_THREADS=6
export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6
export GOTO_NUM_THREADS=6
export NUMEXPR_NUM_THREADS=6

# 1. Rebuild the significant-hits table from the raw permutation results + h5ad
$PY 02_SURVEYOR/export_significant_hits.py

# 2. Classify hits two ways (empirical dominance, then MANE Select role)
$PY classify_hit_scenarios.py
$PY classify_hit_scenarios_mane.py

# 3. Gate down to the two priority candidate groups
$PY 02_SURVEYOR/initial_filter.py

# 4. ASSISTANT_SURVEYOR (L1-L5 enrichment: biotype, AD-prior, Open Targets, UniProt, consequences)
(cd 02_SURVEYOR && $PY -m assistant_surveyor.run_assistant_surveyor)
# add --no-cache above to force-refetch every API response instead of reusing outputs/assistant_surveyor/cache/

# 5. JUNIOR_SURVEYOR (J1-J4: sequence diff, Pfam domains @ 6 CPUs, drug targets, gating)
$PY 02_SURVEYOR/junior_surveyor/run_junior_surveyor.py

# 6. MASTER_SURVEYOR (final plot suite over the J4-gated shortlist)
$PY 02_SURVEYOR/master_surveyor/plot_results.py

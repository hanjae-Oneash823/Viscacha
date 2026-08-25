#!/usr/bin/env bash
set -euo pipefail

cd /home/welcome3/Viscacha_pipeline
conda run -n pocket_dock python 02_SURVEYOR/master_surveyor/run_vina_redock.py \
  --receptor outputs/docking_campaign/CACNA1D_isradipine/prepared/cacna1d_8E59_chain_A.pdbqt \
  --ligand outputs/docking_campaign/CACNA1D_isradipine/prepared/amiodarone_3PE_A2203.pdbqt \
  --outdir outputs/docking_campaign/CACNA1D_isradipine/runs/vina_amiodarone_seed20260824_ex32 \
  --center 137.086 151.020 168.491 \
  --box-size 20 22 30 \
  --seed 20260824 --cpu 16 --exhaustiveness 32 --n-poses 10 \
  2>&1 | tee outputs/docking_campaign/CACNA1D_isradipine/logs/vina_amiodarone_seed20260824_ex32.log

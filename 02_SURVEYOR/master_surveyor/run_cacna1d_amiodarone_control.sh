#!/usr/bin/env bash
set -euo pipefail

cd /home/welcome3/Viscacha_pipeline
conda run -n pocket_dock python 02_SURVEYOR/master_surveyor/run_vina_redock.py \
  --receptor outputs/docking_campaign/systems/CACNA1D_isradipine/prepared/cacna1d_8E59_chain_A.pdbqt \
  --ligand outputs/docking_campaign/systems/CACNA1D_isradipine/prepared/amiodarone_BBI_A2201.pdbqt \
  --outdir outputs/docking_campaign/systems/CACNA1D_isradipine/runs/amiodarone_corrected_seed20260825_ex32 \
  --center 151.334 167.442 149.793 \
  --box-size 18.851 27.189 19.189 \
  --seed 20260825 --cpu 16 --exhaustiveness 32 --n-poses 10 \
  2>&1 | tee outputs/docking_campaign/systems/CACNA1D_isradipine/logs/amiodarone_corrected_seed20260825_ex32.log

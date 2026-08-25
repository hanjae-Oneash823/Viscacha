#!/usr/bin/env bash
# Independent masitinib searches in the validated experimental c-KIT 1T46 pocket.
# Runs serially: never more than 16 Vina CPU cores are allocated in total.
set -euo pipefail

project_dir=${1:?usage: launch_kit_masitinib_template_replicates.sh PROJECT_DIR}
runner="$project_dir/02_SURVEYOR/master_surveyor/run_vina_redock.py"
campaign="$project_dir/outputs/docking_campaign/KIT_masitinib"
receptor="$campaign/prepared/kit_1T46_chain_A.pdbqt"
ligand="$campaign/prepared/masitinib.pdbqt"

run_one() {
    local seed=$1
    local exhaustiveness=$2
    local label="masitinib_1T46_seed${seed}_ex${exhaustiveness}"
    conda run -n pocket_dock python "$runner" \
      --receptor "$receptor" --ligand "$ligand" \
      --outdir "$campaign/runs/$label" \
      --center 26.181 26.061 40.364 --box-size 22 20 28 \
      --seed "$seed" --cpu 16 --exhaustiveness "$exhaustiveness" --n-poses 10 \
      >"$campaign/logs/$label.log" 2>&1
}

for spec in '1103 32' '2207 32' '3301 32' '4409 32' '5519 64' '6619 64' '7723 64' '8831 64'; do
  read -r seed exhaustiveness <<<"$spec"
  run_one "$seed" "$exhaustiveness"
done

#!/usr/bin/env bash
# Eight independent FYN--saracatinib Vina searches: four at exhaustiveness 32
# and four at 64. Four jobs run concurrently, each pinned to 16 Vina CPUs.
set -euo pipefail

project_dir=${1:?usage: launch_fyn_vina_replicates.sh PROJECT_DIR}
runner="$project_dir/02_SURVEYOR/master_surveyor/run_vina_redock.py"
campaign="$project_dir/outputs/docking_campaign/systems/FYN_saracatinib"
receptor="$campaign/prepared/fyn_chain_A.pdbqt"
ligand="$campaign/prepared/saracatinib_H8H_A601.pdbqt"

run_one() {
    local seed=$1
    local exhaustiveness=$2
    local label="vina_seed${seed}_ex${exhaustiveness}"
    conda run -n pocket_dock python "$runner" \
        --receptor "$receptor" --ligand "$ligand" \
        --outdir "$campaign/runs/$label" \
        --center -11.255 14.853 -9.445 --box-size 20 20 26 \
        --seed "$seed" --cpu 16 --exhaustiveness "$exhaustiveness" --n-poses 10 \
        >"$campaign/logs/$label.log" 2>&1
}

active=0
for spec in '1103 32' '2207 32' '3301 32' '4409 32' '5519 64' '6619 64' '7723 64' '8831 64'; do
    read -r seed exhaustiveness <<<"$spec"
    run_one "$seed" "$exhaustiveness" &
    active=$((active + 1))
    if [ "$active" -eq 4 ]; then
        wait
        active=0
    fi
done
wait

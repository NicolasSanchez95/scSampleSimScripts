#!/bin/bash
#SBATCH --job-name=sim_anls
#SBATCH --output=./logs/sim_anls_%A_%a.out
#SBATCH --error=./logs/sim_anls_%A_%a.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
# Array range is set by run_param_creator.sh (default 0-49).

# One simulation replicate. Shared parameters come from the environment
# (set by run_param_creator.sh). id_check = BASE_ID + SLURM_ARRAY_TASK_ID.

set -euo pipefail

# SLURM copies this file to a temp path; do not use dirname(BASH_SOURCE).
load_paper_config() {
  local f
  for f in \
    ${REPO_ROOT:+"${REPO_ROOT}/config.sh"} \
    ${SLURM_SUBMIT_DIR:+"${SLURM_SUBMIT_DIR}/../config.sh"} \
    ${SLURM_SUBMIT_DIR:+"${SLURM_SUBMIT_DIR}/config.sh"} \
    "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../config.sh"
  do
    if [[ -n "$f" && -f "$f" ]]; then
      # shellcheck disable=SC1090
      source "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")"
      return 0
    fi
  done
  echo "Cannot find config.sh (REPO_ROOT=${REPO_ROOT:-unset} SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset})." >&2
  exit 1
}
load_paper_config

SCRIPT_DIR="${REPO_ROOT}/simulation_scripts"

: "${SIM_PREFIX:=augSimLeiden}"
: "${SIM_TYPE:=augData}"
: "${CLUSTER_TYPE:=leiden}"
: "${LEIDEN_RES:=0.1}"
: "${KMEANS_K:=6}"
: "${NUM_GENES_FOR_TVSS_CLUSTERING:=5000}"
: "${RESULTS_SUBDIR:=results}"
: "${NUM_DE_GENES:=100}"
: "${BASE_ID:=1}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "Run via sbatch --array=... (see run_param_creator.sh)." >&2
  exit 1
fi

id_check=$((BASE_ID + SLURM_ARRAY_TASK_ID))

cd "${SCRIPT_DIR}"
mkdir -p logs

echo "Running analysis: num_de_genes=${NUM_DE_GENES}, id_check=${id_check}, sim_prefix=${SIM_PREFIX}, cluster_type=${CLUSTER_TYPE}"
echo "Repo: ${REPO_ROOT}"

"${RSCRIPT}" sim_analysis_pipeline.R \
  --repo-root "${REPO_ROOT}" \
  --sim_type "${SIM_TYPE}" \
  --sim_prefix "${SIM_PREFIX}" \
  --num_de_genes "${NUM_DE_GENES}" \
  --cluster_type "${CLUSTER_TYPE}" \
  --leiden_res "${LEIDEN_RES}" \
  --kmeans_k "${KMEANS_K}" \
  --id_check "${id_check}" \
  --augDataRep "${SIM_PREFIX}" \
  --num_genes_for_TvsS_clustering "${NUM_GENES_FOR_TVSS_CLUSTERING}" \
  --run_PVE_metrics FALSE \
  --results_subdir "${RESULTS_SUBDIR}"

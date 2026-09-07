#!/bin/bash
#SBATCH --job-name=master_conf_pwr
#SBATCH --output=./logs/master_conf_pwr_%A_%a.out
#SBATCH --error=./logs/master_conf_pwr_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-60%15

# Per-phenotype pipeline: array of phentype_rem_conf_table.R tasks, then
# one combine_batch_results.R. Array size is overridden on submit.

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

SCRIPT_DIR="${REPO_ROOT}/master_table_scripts"

args=("$@")
parse_flag() {
  local flag="$1"
  local i
  for i in "${!args[@]}"; do
    if [ "${args[i]}" = "$flag" ] && [ $((i + 1)) -lt ${#args[@]} ]; then
      echo "${args[i+1]}"
      return
    fi
  done
  echo ""
}

if [ "$#" -ge 1 ] && [ "$1" = "--repo-root" ]; then
  REPO_ROOT=$(parse_flag "--repo-root")
  PHEN_TYPE=$(parse_flag "--phen-type-removal")
  COHORT_NAME=$(parse_flag "--cohort-name")
  COHORT_DATA_DIR=$(parse_flag "--cohort-data-dir")
  SIM_RESULTS_DIR=$(parse_flag "--sim-results-dir")
  RESULTS_BASE=$(parse_flag "--results-base")
else
  PHEN_TYPE="${1:?phen_type required}"
  COHORT_NAME="${2:?cohort_name required}"
  COHORT_DATA_DIR="${COHORT_DATA_DIR:-}"
  SIM_RESULTS_DIR="${SIM_RESULTS_DIR:-}"
  RESULTS_BASE="${RESULTS_BASE:-}"
fi

[ -z "${COHORT_DATA_DIR:-}" ] && COHORT_DATA_DIR="${REPO_ROOT}/master_table_scripts/data"
[ -z "${SIM_RESULTS_DIR:-}" ] && SIM_RESULTS_DIR="${REPO_ROOT}/simulation_scripts/results"
[ -z "${RESULTS_BASE:-}" ] && RESULTS_BASE="${REPO_ROOT}/master_table_scripts/results"

mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

if [ -z "${SLURM_JOB_ID:-}" ]; then
  COHORT_IDS_RDS="${COHORT_DATA_DIR}/cohort_ids_${COHORT_NAME}.rds"
  if [ ! -f "$COHORT_IDS_RDS" ]; then
    echo "Error: cohort_ids file not found: $COHORT_IDS_RDS" >&2
    exit 1
  fi
  N=$("${RSCRIPT}" -e "cat(ceiling(length(readRDS('$COHORT_IDS_RDS'))/10))")
  if [ -z "$N" ] || [ "$N" -lt 1 ]; then
    echo "Error: could not determine array size from $COHORT_IDS_RDS" >&2
    exit 1
  fi
  echo "Submitting array job for phenotype: $PHEN_TYPE, cohort: $COHORT_NAME (array 1-$N)"
  mapfile -t SBATCH_EXTRA < <(sbatch_extra_args)
  ARRAY_JOB=$(sbatch "${SBATCH_EXTRA[@]}" --chdir="${SCRIPT_DIR}" --array=1-${N}%15 \
    --export=ALL,REPO_ROOT="$REPO_ROOT",COHORT_DATA_DIR="$COHORT_DATA_DIR",SIM_RESULTS_DIR="$SIM_RESULTS_DIR",RESULTS_BASE="$RESULTS_BASE" \
    "${REPO_ROOT}/master_table_scripts/run_single_phentype_batch.sh" "$PHEN_TYPE" "$COHORT_NAME" | awk '{print $4}')
  echo "Array job ID: $ARRAY_JOB"
  echo "Submitting combination job (runs after array)..."
  sbatch "${SBATCH_EXTRA[@]}" --chdir="${SCRIPT_DIR}" --dependency=afterok:"$ARRAY_JOB" --array=0 \
    --export=ALL,MASTER_ARRAY_JOB_ID="$ARRAY_JOB",REPO_ROOT="$REPO_ROOT",COHORT_DATA_DIR="$COHORT_DATA_DIR",SIM_RESULTS_DIR="$SIM_RESULTS_DIR",RESULTS_BASE="$RESULTS_BASE" \
    "${REPO_ROOT}/master_table_scripts/run_single_phentype_batch.sh" "$PHEN_TYPE" "$COHORT_NAME"
  echo "Done. Combination runs after all array tasks complete."
  exit 0
fi

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && [ "$SLURM_ARRAY_TASK_ID" -gt 0 ]; then
  echo "Array task $SLURM_ARRAY_TASK_ID: phentype_rem_conf_table.R $PHEN_TYPE $COHORT_NAME"
  "${RSCRIPT}" --vanilla "${REPO_ROOT}/master_table_scripts/phentype_rem_conf_table.R" \
    --repo-root "$REPO_ROOT" \
    --phen-type-removal "$PHEN_TYPE" \
    --cohort-name "$COHORT_NAME" \
    --cohort-data-dir "$COHORT_DATA_DIR" \
    --sim-results-dir "$SIM_RESULTS_DIR" \
    --results-base "$RESULTS_BASE" \
    --array-job-id "${SLURM_ARRAY_JOB_ID}" \
    --array-task-id "${SLURM_ARRAY_TASK_ID}" \
    2>&1 | tee -a "logs/master_conf_pwr_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
  echo "Array task $SLURM_ARRAY_TASK_ID completed"
else
  echo "Running combine_batch_results.R for $PHEN_TYPE (cohort: $COHORT_NAME, array job: ${MASTER_ARRAY_JOB_ID:-})"
  "${RSCRIPT}" "${REPO_ROOT}/master_table_scripts/combine_batch_results.R" \
    "$REPO_ROOT" "$PHEN_TYPE" "$COHORT_NAME" "${MASTER_ARRAY_JOB_ID:-}" "$RESULTS_BASE"
  echo "Combination completed"
fi

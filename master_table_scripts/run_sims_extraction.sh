#!/bin/bash
# Build a cohort from simulation outputs, then submit per-phenotype extraction jobs.
#
# Usage:
#   ./run_sims_extraction.sh --cohort-name demo50 --id-pattern "augSimLeiden_" \
#     [--repo-root PATH] [--sim-results-dir DIR] [--cohort-data-dir DIR] [--results-base DIR]

set -euo pipefail

# Same locator as the sbatch wrappers (safe if this script is also submitted).
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

usage() {
  cat <<EOF
Usage:
  ./run_sims_extraction.sh --cohort-name <name> --id-pattern <pattern> \\
    [--repo-root <path>] [--sim-results-dir <dir>] [--cohort-data-dir <dir>] [--results-base <dir>]

Example:
  ./run_sims_extraction.sh --cohort-name demo50 --id-pattern "augSimLeiden_"

Notes:
  --sim-results-dir must contain anls_info/ (and de_final/ for extraction).
  Default repo root: ${REPO_ROOT}
EOF
}

COHORT_NAME=""
ID_PATTERN=""
COHORT_DATA_DIR=""
SIM_RESULTS_DIR=""
RESULTS_BASE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-root) REPO_ROOT="${2:-}"; shift 2 ;;
    --cohort-name) COHORT_NAME="${2:-}"; shift 2 ;;
    --id-pattern) ID_PATTERN="${2:-}"; shift 2 ;;
    --cohort-data-dir) COHORT_DATA_DIR="${2:-}"; shift 2 ;;
    --sim-results-dir) SIM_RESULTS_DIR="${2:-}"; shift 2 ;;
    --results-base) RESULTS_BASE="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$COHORT_NAME" || -z "$ID_PATTERN" ]]; then
  usage
  exit 2
fi

if [[ -z "$SIM_RESULTS_DIR" ]]; then
  SIM_RESULTS_DIR="${REPO_ROOT}/simulation_scripts/results"
fi
if [[ -z "$COHORT_DATA_DIR" ]]; then
  COHORT_DATA_DIR="${REPO_ROOT}/master_table_scripts/data"
fi
if [[ -z "$RESULTS_BASE" ]]; then
  RESULTS_BASE="${REPO_ROOT}/master_table_scripts/results"
fi

mkdir -p "${COHORT_DATA_DIR}" "${RESULTS_BASE}" "${SCRIPT_DIR}/logs"

export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

"${RSCRIPT}" --vanilla "${REPO_ROOT}/master_table_scripts/cohort_builder.R" \
  --repo-root "$REPO_ROOT" \
  --cohort-name "$COHORT_NAME" \
  --id-pattern "$ID_PATTERN" \
  --sim-results-dir "$SIM_RESULTS_DIR" \
  --cohort-data-dir "$COHORT_DATA_DIR" \
  --results-base "$RESULTS_BASE"

PHENOTYPES=(
  "noharm_"
  "noharm_null_"
  "harmony_null_"
  "harmony_"
  "TvsS_"
  "fmnn_"
)

echo "Phenotypes: ${PHENOTYPES[*]}"
echo "Submitting ${#PHENOTYPES[@]} jobs..."

cd "${SCRIPT_DIR}"
for phen in "${PHENOTYPES[@]}"; do
  echo "Submitting master_conf_pwr for phenotype: $phen (cohort: $COHORT_NAME)"
  ./run_single_phentype_batch.sh \
    --repo-root "$REPO_ROOT" \
    --phen-type-removal "$phen" \
    --cohort-name "$COHORT_NAME" \
    --cohort-data-dir "$COHORT_DATA_DIR" \
    --sim-results-dir "$SIM_RESULTS_DIR" \
    --results-base "$RESULTS_BASE"
done

echo "Done. Each phenotype runs as an array job plus a combination job."

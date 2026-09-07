#!/bin/bash
#SBATCH --job-name=allHyp_metrics
#SBATCH --output=./logs/allHyp_metrics.out
#SBATCH --error=./logs/allHyp_metrics.err
#SBATCH --mem=10G

# Compute all-hypothesis conf metrics (overall and by cluster) for one cohort.
#
# Usage:
#   ./run_allHyp_metrics.sh --cohort-name demo50 \
#     [--repo-root PATH] [--overlap-def union] [--cut-off-true 0] [--cut-off-false 0] \
#     [--sig-threshold 0.1] [--p-val-col adj_p_val]

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

SCRIPT_DIR="${REPO_ROOT}/decision_parameters"

COHORT_NAME=""
OVERLAP_DEF="union"
CUT_OFF_TRUE="0.0"
CUT_OFF_FALSE="0.0"
SIG_THRESHOLD="0.1"
P_VAL_COL="adj_p_val"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-root) REPO_ROOT="$2"; shift 2 ;;
    --cohort-name) COHORT_NAME="$2"; shift 2 ;;
    --overlap-def) OVERLAP_DEF="$2"; shift 2 ;;
    --cut-off-true) CUT_OFF_TRUE="$2"; shift 2 ;;
    --cut-off-false) CUT_OFF_FALSE="$2"; shift 2 ;;
    --sig-threshold) SIG_THRESHOLD="$2"; shift 2 ;;
    --p-val-col) P_VAL_COL="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --cohort-name NAME [--repo-root PATH] [--overlap-def union] ..."
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$COHORT_NAME" ]]; then
  echo "Required: --cohort-name" >&2
  exit 2
fi

mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

echo "Computing all-hypothesis metrics for cohort ${COHORT_NAME}"
"${RSCRIPT}" compute_allHyp_metrics.R \
  "${REPO_ROOT}" "${COHORT_NAME}" "${OVERLAP_DEF}" "${CUT_OFF_TRUE}" "${CUT_OFF_FALSE}" \
  "${SIG_THRESHOLD}" "${P_VAL_COL}"

"${RSCRIPT}" compute_allHyp_metrics_by_cluster.R \
  "${REPO_ROOT}" "${COHORT_NAME}" "${OVERLAP_DEF}" "${CUT_OFF_TRUE}" "${CUT_OFF_FALSE}" \
  "${SIG_THRESHOLD}" "${P_VAL_COL}"

echo "Done."

#!/bin/bash
#SBATCH --job-name=param_creator
#SBATCH --output=logs/param_creator_%j.out
#SBATCH --error=logs/param_creator_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00

# Orchestrator: run a param creator, then submit simulation array jobs
# (one num_de_genes). Safe to sbatch or run on a login/compute node.
#
# Usage (from this directory; pass -p/--partition if your site requires it):
#   sbatch -p PARTITION run_param_creator.sh --sim_prefix augSimLeiden \
#     --sim_type augData --cluster_type leiden --n_sims 50 --num_de_genes 100
#
# mkdir -p logs first so SLURM can open the log files.

set -euo pipefail

# SLURM copies this file to a temp path, so BASH_SOURCE cannot find config.sh.
# Fall back to the directory from which sbatch was invoked.
resolve_script_dir() {
  local src="${BASH_SOURCE[0]}"
  local dir
  dir="$(cd "$(dirname "$src")" && pwd)"
  if [[ -f "${dir}/../config.sh" ]]; then
    echo "$dir"
    return 0
  fi
  local submit="${SLURM_SUBMIT_DIR:-}"
  if [[ -n "$submit" && -f "${submit}/../config.sh" ]]; then
    echo "$(cd "$submit" && pwd)"
    return 0
  fi
  if [[ -n "$submit" && -f "${submit}/config.sh" && -d "${submit}/simulation_scripts" ]]; then
    echo "$(cd "${submit}/simulation_scripts" && pwd)"
    return 0
  fi
  echo "Cannot find config.sh (script dir=${dir}, SLURM_SUBMIT_DIR=${submit:-unset})." >&2
  echo "Submit from simulation_scripts/ or the repo root." >&2
  exit 1
}

SCRIPT_DIR="$(resolve_script_dir)"
# shellcheck source=../config.sh
source "${SCRIPT_DIR}/../config.sh"

# Reuse this job's partition/account for the simulation array if the user
# passed sbatch -p / --account rather than exporting SLURM_PARTITION.
if [[ -z "${SLURM_PARTITION:-}" && -n "${SLURM_JOB_PARTITION:-}" ]]; then
  SLURM_PARTITION="${SLURM_JOB_PARTITION}"
fi
if [[ -z "${SLURM_ACCOUNT:-}" && -n "${SLURM_JOB_ACCOUNT:-}" ]]; then
  SLURM_ACCOUNT="${SLURM_JOB_ACCOUNT}"
fi

SIM_PREFIX="${SIM_PREFIX:-augSimLeiden}"
SIM_TYPE="${SIM_TYPE:-augData}"
CLUSTER_TYPE="${CLUSTER_TYPE:-leiden}"
LEIDEN_RES="${LEIDEN_RES:-0.1}"
KMEANS_K="${KMEANS_K:-6}"
NUM_GENES_FOR_TVSS_CLUSTERING="${NUM_GENES_FOR_TVSS_CLUSTERING:-5000}"
RESULTS_SUBDIR="${RESULTS_SUBDIR:-results}"
N_SIMS="${N_SIMS:-50}"
NUM_DE_GENES="${NUM_DE_GENES:-100}"
BASE_ID="${BASE_ID:-1}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sim_prefix) SIM_PREFIX="$2"; shift 2 ;;
    --sim_type) SIM_TYPE="$2"; shift 2 ;;
    --cluster_type) CLUSTER_TYPE="$2"; shift 2 ;;
    --leiden_res) LEIDEN_RES="$2"; shift 2 ;;
    --kmeans_k) KMEANS_K="$2"; shift 2 ;;
    --num_genes_for_TvsS_clustering) NUM_GENES_FOR_TVSS_CLUSTERING="$2"; shift 2 ;;
    --results_subdir) RESULTS_SUBDIR="$2"; shift 2 ;;
    --n_sims) N_SIMS="$2"; shift 2 ;;
    --num_de_genes) NUM_DE_GENES="$2"; shift 2 ;;
    --base_id) BASE_ID="$2"; shift 2 ;;
    --repo-root) REPO_ROOT="$2"; shift 2 ;;
    *)
      if [[ "$1" != --* ]]; then
        SIM_PREFIX="$1"
      else
        echo "Unknown argument: $1" >&2
        exit 2
      fi
      shift
      ;;
  esac
done

PARAM_CREATOR="${SCRIPT_DIR}/param_creator_${SIM_PREFIX}.R"
if [[ ! -f "${PARAM_CREATOR}" ]]; then
  echo "Missing param creator: ${PARAM_CREATOR}" >&2
  echo "Expected param_creator_<sim_prefix>.R (e.g. param_creator_augSimLeiden.R)." >&2
  exit 1
fi

mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

echo "Using: sim_prefix=${SIM_PREFIX} sim_type=${SIM_TYPE} cluster_type=${CLUSTER_TYPE}"
echo "       n_sims=${N_SIMS} num_de_genes=${NUM_DE_GENES} base_id=${BASE_ID} results_subdir=${RESULTS_SUBDIR}"
echo "Repo: ${REPO_ROOT}"
echo "Partition: ${SLURM_PARTITION:-<default>}"

echo "Running ${PARAM_CREATOR} ..."
"${RSCRIPT}" "${PARAM_CREATOR}" \
  --sim_prefix "${SIM_PREFIX}" \
  --leiden_res "${LEIDEN_RES}" \
  --results_subdir "${RESULTS_SUBDIR}" \
  --repo-root "${REPO_ROOT}"

echo "Submitting run_sim_batch.sh as array 0-$((N_SIMS - 1)) ..."
ARRAY_MAX=$((N_SIMS - 1))
if [[ ${ARRAY_MAX} -lt 0 ]]; then
  echo "N_SIMS must be >= 1" >&2
  exit 2
fi

mapfile -t SBATCH_EXTRA < <(sbatch_extra_args)
sbatch "${SBATCH_EXTRA[@]}" \
  --chdir="${SCRIPT_DIR}" \
  --array="0-${ARRAY_MAX}%25" \
  --export=ALL,SIM_PREFIX="${SIM_PREFIX}",SIM_TYPE="${SIM_TYPE}",CLUSTER_TYPE="${CLUSTER_TYPE}",LEIDEN_RES="${LEIDEN_RES}",KMEANS_K="${KMEANS_K}",NUM_GENES_FOR_TVSS_CLUSTERING="${NUM_GENES_FOR_TVSS_CLUSTERING}",RESULTS_SUBDIR="${RESULTS_SUBDIR}",NUM_DE_GENES="${NUM_DE_GENES}",BASE_ID="${BASE_ID}",REPO_ROOT="${REPO_ROOT}" \
  "${SCRIPT_DIR}/run_sim_batch.sh"

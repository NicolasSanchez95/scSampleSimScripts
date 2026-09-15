#!/usr/bin/env bash
# Shared paths and SLURM/R settings for the paper-facing scripts.
# Source from other bash wrappers:  source "$(dirname "$0")/../config.sh"
# or, from this directory:          source ./config.sh
#
# SLURM copies #SBATCH scripts to a temp path, so dirname(BASH_SOURCE) is not
# the repo. Sbatch'd wrappers should locate this file via REPO_ROOT or
# SLURM_SUBMIT_DIR (see load_paper_config in the batch shells).

_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-${_CONFIG_DIR}}"

# Pinned scSampleSim release used by all R entrypoints (must match R/ensure_scSampleSim.R).
# Do not reuse tested-backup-20260904 / ba61526 (that commit is still Package: postselect).
# After the Package: scSampleSim rename is committed, tag it scSampleSim-20260910 and
# set SCSAMPLESIM_SHA to that commit's full SHA.
SCSAMPLESIM_REPO="${SCSAMPLESIM_REPO:-epurdom/scSampleSim}"
SCSAMPLESIM_REF="${SCSAMPLESIM_REF:-scSampleSim-20260910}"
SCSAMPLESIM_SHA="${SCSAMPLESIM_SHA:-PENDING_AFTER_SCSAMPLESIM_RENAME_TAG}"

# Rscript: R_SCRIPT env override, else Rscript on PATH.
# Optional: module load R/4.5.0 on clusters that provide Environment Modules.
if [[ -n "${R_SCRIPT:-}" && -x "${R_SCRIPT}" ]]; then
  RSCRIPT="${R_SCRIPT}"
else
  RSCRIPT="$(command -v Rscript || true)"
fi
if [[ -z "${RSCRIPT}" ]]; then
  echo "Rscript not found on PATH. Set R_SCRIPT to an Rscript executable." >&2
  exit 1
fi

# Do not override an existing user library; only export if the caller set it.
# (Leave R_LIBS_USER unset to use the default user library.)

sbatch_extra_args() {
  local args=()
  if [[ -n "${SLURM_PARTITION:-}" ]]; then
    args+=(--partition="${SLURM_PARTITION}")
  fi
  if [[ -n "${SLURM_ACCOUNT:-}" ]]; then
    args+=(--account="${SLURM_ACCOUNT}")
  fi
  printf '%s\n' "${args[@]}"
}

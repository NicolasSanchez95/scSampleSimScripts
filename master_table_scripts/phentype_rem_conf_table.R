#!/usr/bin/env Rscript
# Phenotype-removal confidence table extraction (SLURM array).
# Each array task processes a chunk of id_checks and saves a temp RDS.

library(pbapply)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage:",
  "  Rscript phentype_rem_conf_table.R --repo-root <path> --phen-type-removal <name>",
  "    --cohort-name <name> --cohort-data-dir <dir>",
  "",
  "Optional:",
  "  --sim-results-dir <dir>   Default: <repo-root>/simulation_scripts/results",
  "  --results-base <dir>      Default: <repo-root>/master_table_scripts/results",
  sep = "\n"
)

parse_flag <- function(flag, args) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) return(NULL)
  args[i + 1]
}

.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .script_args, value = TRUE)
.this_script <- if (length(.file_arg) >= 1L) {
  normalizePath(sub("^--file=", "", .file_arg[[1]]), mustWork = TRUE)
} else {
  stop("Run this file with Rscript so --file= is set.")
}

repo_root <- parse_flag("--repo-root", args)
if (is.null(repo_root) || !nzchar(repo_root)) {
  repo_root <- normalizePath(file.path(dirname(.this_script), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
}
source(file.path(repo_root, "R", "ensure_postselect.R"))
ensure_postselect()

phen_type_removal <- parse_flag("--phen-type-removal", args)
cohort_name <- parse_flag("--cohort-name", args)
cohort_data_dir <- parse_flag("--cohort-data-dir", args)
sim_results_dir <- parse_flag("--sim-results-dir", args)
results_base <- parse_flag("--results-base", args)
array_job_id_arg <- parse_flag("--array-job-id", args)
array_task_id_arg <- parse_flag("--array-task-id", args)
if (is.null(repo_root) || is.null(phen_type_removal) || is.null(cohort_name) || is.null(cohort_data_dir)) {
  stop(usage)
}
cat("Running for phen_type_removal:", phen_type_removal, ", cohort_name:", cohort_name, "\n")

data_dir <- cohort_data_dir
if (is.null(sim_results_dir)) sim_results_dir <- file.path(repo_root, "simulation_scripts", "results")
if (is.null(results_base)) results_base <- file.path(repo_root, "master_table_scripts", "results")
if (!dir.exists(sim_results_dir)) stop(paste0("Directory does not exist: ", sim_results_dir))
dir.create(results_base, recursive = TRUE, showWarnings = FALSE)
analysis_results_dir <- file.path(sim_results_dir, "anls_info")
de_outputs_dir <- file.path(sim_results_dir, "de_final")
# per_sim_data_extraction concatenates prefixes with paste0; keep trailing slashes.
analysis_results_dir <- paste0(analysis_results_dir, "/")
de_outputs_dir <- paste0(de_outputs_dir, "/")

overlap_matrix_prefix <- "de_overlap_info_"
de_pvals_by_cluster_prefix <- "res_de_"
pa_de_save_prfx <- "pa_de_int_"
mv_PVE_metrics_save_prfx <- "mv_PVE_embed_"

chunk_size <- 10

array_job_id <- array_job_id_arg
if (is.null(array_job_id) || !nzchar(array_job_id)) {
  array_job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")
}
array_task_id <- if (!is.null(array_task_id_arg) && nzchar(array_task_id_arg)) {
  as.integer(array_task_id_arg)
} else {
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}
if (!nzchar(array_job_id) || is.na(array_task_id)) {
  stop("This script must be run as a SLURM array job. Pass --array-job-id and --array-task-id (or set SLURM_ARRAY_JOB_ID / SLURM_ARRAY_TASK_ID).")
}
cat(sprintf("Running as SLURM array job: Job ID %s, Task ID %d\n", array_job_id, array_task_id))

output_dir <- file.path(results_base, cohort_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

temp_output_dir <- file.path(output_dir, "temp_pwr_table", array_job_id)
if (!dir.exists(data_dir)) {
  stop(paste0("Directory does not exist: ", data_dir))
}
cohort_ids_path <- file.path(data_dir, paste0("cohort_ids_", cohort_name, ".rds"))
if (!file.exists(cohort_ids_path)) {
  stop(paste0("cohort_ids file not found: ", cohort_ids_path))
}
cat("Loading cohort ids from:", cohort_ids_path, "\n")
cohort_ids <- readRDS(cohort_ids_path)
if (is.list(cohort_ids) && length(cohort_ids) == 1L && is.character(cohort_ids[[1]])) {
  cohort_ids <- cohort_ids[[1]]
}
if (is.data.frame(cohort_ids)) {
  if ("id_check" %in% names(cohort_ids)) {
    cohort_ids <- cohort_ids[["id_check"]]
  } else {
    cohort_ids <- cohort_ids[[1]]
  }
}
cat("cohort_ids class:", paste(class(cohort_ids), collapse = ","), "\n")
cat("cohort_ids length:", length(cohort_ids), "\n")

options(pbapply.progress.type = "txt")
options(pbapply.nout = 1000L)

dir.create(temp_output_dir, recursive = TRUE, showWarnings = FALSE)

total_cohort_ids <- length(cohort_ids)
n_chunks <- ceiling(total_cohort_ids / chunk_size)
chunk_start <- (array_task_id - 1) * chunk_size + 1
chunk_end <- min(array_task_id * chunk_size, total_cohort_ids)

if (chunk_start > total_cohort_ids) {
  cat(sprintf("Array task ID %d exceeds number of chunks (%d). Exiting.\n", array_task_id, n_chunks))
  quit(status = 0)
}

cohort_ids_to_process <- cohort_ids[chunk_start:chunk_end]
cat(sprintf(
  "Array task %d: Processing id_checks %d-%d (chunk %d of %d, %d items)\n",
  array_task_id, chunk_start, chunk_end, array_task_id, n_chunks, length(cohort_ids_to_process)
))
flush.console()

start_time <- Sys.time()
all_results <- pbapply::pblapply(
  cohort_ids_to_process,
  function(x) {
    result <- postselect::per_sim_data_extraction(
      id_check_single = x,
      phen_type_removal = phen_type_removal,
      analysis_results_dir = analysis_results_dir,
      de_outputs_dir = de_outputs_dir,
      overlap_matrix_prefix = overlap_matrix_prefix,
      de_pvals_by_cluster_prefix = de_pvals_by_cluster_prefix,
      pa_de_save_prfx = pa_de_save_prfx,
      mv_PVE_metrics_save_prfx = mv_PVE_metrics_save_prfx
    )
    flush.console()
    result
  },
  cl = NULL
)

chunk_combined <- do.call(rbind, all_results)
output_file <- file.path(
  temp_output_dir,
  paste0("temp_result_", phen_type_removal, "_cohort_", cohort_name, "_task_", array_task_id, ".rds")
)
saveRDS(chunk_combined, output_file)

cat(sprintf("Saved chunk result (%d rows) to %s\n", nrow(chunk_combined), output_file))
cat("Time taken for chunk:", format(Sys.time() - start_time), "\n")
cat(sprintf("Array task %d completed successfully.\n", array_task_id))

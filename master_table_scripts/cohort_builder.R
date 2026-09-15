#!/usr/bin/env Rscript
# Build a cohort ID vector from simulation analysis RDS filenames.

library(stringr)

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage:",
  "  Rscript cohort_builder.R --cohort-name <name> --id-pattern <pattern>",
  "    [--repo-root <path>] [--sim-results-dir <dir>] [--cohort-data-dir <dir>]",
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
source(file.path(repo_root, "R", "ensure_scSampleSim.R"))
ensure_scSampleSim()

cohort_name <- parse_flag("--cohort-name", args)
id_pattern <- parse_flag("--id-pattern", args)
sim_results_dir <- parse_flag("--sim-results-dir", args)
cohort_data_dir <- parse_flag("--cohort-data-dir", args)
results_base <- parse_flag("--results-base", args)

print("Building cohort for the following parameters:")
print(paste0("cohort_name: ", cohort_name))
print(paste0("id_pattern: ", id_pattern))
print(paste0("sim_results_dir: ", sim_results_dir))
print(paste0("cohort_data_dir: ", cohort_data_dir))
print(paste0("repo_root: ", repo_root))
print(paste0("results_base (unused here): ", results_base))

if (is.null(cohort_name) || is.null(id_pattern)) {
  stop(usage)
}

if (is.null(sim_results_dir) || is.null(cohort_data_dir)) {
  if (is.null(sim_results_dir)) {
    sim_results_dir <- file.path(repo_root, "simulation_scripts", "results")
  }
  if (is.null(cohort_data_dir)) {
    cohort_data_dir <- file.path(repo_root, "master_table_scripts", "data")
  }
}
dir.create(cohort_data_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(sim_results_dir)) {
  stop(paste0("Directory does not exist: ", sim_results_dir))
}

sims_anls_info_dir <- file.path(sim_results_dir, "anls_info")
if (!dir.exists(sims_anls_info_dir)) {
  stop(paste0("Directory does not exist: ", sims_anls_info_dir))
}

cohort_prefix <- "cohort_ids_"

sims_pa_de_int_files <- list.files(
  path = sims_anls_info_dir,
  pattern = id_pattern,
  full.names = TRUE
)
sims_pa_de_int_files_df <- data.frame(
  file_name = sims_pa_de_int_files,
  id_check = stringr::str_extract(sims_pa_de_int_files, "(?<=pa_de_int_).*(?=\\.rds)")
)
cohort_ids <- unique(sims_pa_de_int_files_df$id_check)
cohort_ids <- cohort_ids[!is.na(cohort_ids)]

print(paste0("Found ", length(cohort_ids), " cohort ids"))
cohort_ids_path <- file.path(cohort_data_dir, paste0(cohort_prefix, cohort_name, ".rds"))
print(paste0("Saving cohort ids to file: ", cohort_ids_path))
saveRDS(cohort_ids, file = cohort_ids_path)

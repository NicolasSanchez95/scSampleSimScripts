#!/usr/bin/env Rscript
# All-hypothesis-level FDR (FDP) and power metrics.
#
# Reads master_pwr_table_<phen_type>.rds from master_table_scripts/results/<cohort>/
# and writes one RDS (list by phenotype method) into decision_parameters/results/<cohort>/.
#
# Usage: Rscript compute_allHyp_metrics.R <repo_root> <cohort_name> <overlap_def>
#          <cut_off_true> <cut_off_false> <sig_threshold_fixed> <p_val_col_name>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript compute_allHyp_metrics.R <repo_root> <cohort_name> <overlap_def> <cut_off_true> <cut_off_false> <sig_threshold_fixed> <p_val_col_name>")
}
repo_root <- normalizePath(args[[1]], mustWork = TRUE)
cohort_name <- args[[2]]
overlap_def <- args[[3]]
cut_off_true <- as.numeric(args[[4]])
cut_off_false <- as.numeric(args[[5]])
sig_threshold_fixed <- as.numeric(args[[6]])
p_val_col_name <- args[[7]]

source(file.path(repo_root, "R", "ensure_scSampleSim.R"))
ensure_scSampleSim()
library(dplyr)

phen_type_removal_set <- c("fmnn_", "harmony_", "noharm_", "noharm_null_", "harmony_null_", "TvsS_")

master_table_results_dir <- file.path(repo_root, "master_table_scripts", "results", cohort_name)
decision_parameters_results_dir <- file.path(repo_root, "decision_parameters", "results", cohort_name)
dir.create(decision_parameters_results_dir, recursive = TRUE, showWarnings = FALSE)

list_of_conf_metrics_tested <- list()
saved_with_params <- paste0(
  "conf_metrics_tested_", overlap_def, "_", cut_off_true, "_", cut_off_false, "_",
  sig_threshold_fixed, "_", p_val_col_name, ".rds"
)
for (phen_type_removal in phen_type_removal_set) {
  start_time <- Sys.time()
  print(paste0("Processing phen_type_removal: ", phen_type_removal))
  master_pwr_table_fn <- file.path(master_table_results_dir, paste0("master_pwr_table_", phen_type_removal, ".rds"))
  if (!file.exists(master_pwr_table_fn)) {
    warning("Missing master table, skipping: ", master_pwr_table_fn)
    next
  }
  master_pwr_table <- readRDS(master_pwr_table_fn)

  conf_metrics_tested <- scSampleSim::get_conf_metrics_tested_all_hypotheses(
    master_pwr_table_used = master_pwr_table,
    sig_threshold_fixed = sig_threshold_fixed,
    p_val_col_name = p_val_col_name,
    overlap_def = overlap_def,
    cut_off_true = cut_off_true,
    cut_off_false = cut_off_false
  )

  print(paste0("Number of NAs in conf_metrics_tested: ", paste(colSums(is.na(conf_metrics_tested)), collapse = ", ")))
  list_of_conf_metrics_tested[[phen_type_removal]] <- conf_metrics_tested
  rm(master_pwr_table)
  gc()
  end_time <- Sys.time()
  print(paste0("Time taken for phen_type_removal: ", phen_type_removal, " ", end_time - start_time))
}
saveRDS(list_of_conf_metrics_tested, file.path(decision_parameters_results_dir, saved_with_params))
print(paste0("Saved file ", file.path(decision_parameters_results_dir, saved_with_params)))

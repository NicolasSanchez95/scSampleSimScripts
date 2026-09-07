#!/usr/bin/env Rscript
# Combine temp RDS chunks from phentype_rem_conf_table.R into one master table.
# Usage: Rscript combine_batch_results.R <repo_root> <phen_type_removal> <cohort_name> [array_job_id] [results_base]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript combine_batch_results.R <repo_root> <phen_type_removal> <cohort_name> [array_job_id] [results_base]")
}
repo_root <- args[[1]]
phen_type_removal <- args[[2]]
cohort_name <- args[[3]]
array_job_id <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else ""
results_base <- if (length(args) >= 5 && nzchar(args[[5]])) {
  args[[5]]
} else {
  file.path(repo_root, "master_table_scripts", "results")
}

output_dir <- file.path(results_base, cohort_name)

cat("phen_type_removal:", phen_type_removal, "\n")
cat("cohort_name:", cohort_name, "\n")
cat("array_job_id:", array_job_id, "\n")
cat("results_base:", results_base, "\n")
cat("output_dir:", output_dir, "\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

temp_output_dir <- if (nzchar(array_job_id)) {
  file.path(output_dir, "temp_pwr_table", array_job_id)
} else {
  file.path(output_dir, "temp_pwr_table")
}
cat("temp_output_dir:", temp_output_dir, "\n")

temp_pattern <- paste0("temp_result_", phen_type_removal)
temp_files <- list.files(temp_output_dir, pattern = temp_pattern, full.names = TRUE)

cat(sprintf("Found %d temporary result files\n", length(temp_files)))

all_results <- list()
for (file in temp_files) {
  cat(sprintf("Reading %s...\n", basename(file)))
  result <- readRDS(file)
  if (!is.null(result)) {
    all_results[[length(all_results) + 1]] <- result
  }
}

if (length(all_results) == 0) {
  stop("No valid results to combine")
}

master_pwr_combined_sim_res_DE_all <- do.call(rbind, all_results)

final_output_file <- file.path(output_dir, paste0("master_pwr_table_", phen_type_removal, ".rds"))
saveRDS(master_pwr_combined_sim_res_DE_all, final_output_file)
cat(sprintf("Saved combined results to %s\n", final_output_file))
cat(sprintf("Total rows: %d\n", nrow(master_pwr_combined_sim_res_DE_all)))

cat("Cleaning up temporary files...\n")
for (file in temp_files) {
  file.remove(file)
}
cat("Done!\n")

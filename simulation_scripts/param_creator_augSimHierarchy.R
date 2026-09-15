#' Example parameter setup for k-means hierarchy augmented simulations.
#'
#' Builds the base Seurat object and `pa_preseurat_<sim_prefix>.Rds` used by
#' `sim_analysis_pipeline.R`. Place the input SCE at
#' `simulation_scripts/data/filtered_sce_data.Rda` with colData columns
#' `group`, `sample`, `batch`, and `cell_type`.
#'
#' Hierarchy-specific settings: `clustering = "kmeans"`, `cut_at = c(3, 10, 20)`,
#' and a vector of mean log-fold-changes.
#'
#' Usage:
#'   Rscript param_creator_augSimHierarchy.R [--sim_prefix NAME] [--leiden_res 0.1]
#'     [--results_subdir results] [--repo-root PATH]

.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .script_args, value = TRUE)
if (length(.file_arg) < 1L) {
  stop("Run this file with Rscript so --file= is set.")
}
.this_script <- normalizePath(sub("^--file=", "", .file_arg[[1]]), mustWork = TRUE)
.repo_root_default <- normalizePath(file.path(dirname(.this_script), ".."), mustWork = TRUE)

library(mclust)
library(Seurat)
library(DESeq2)
library(scuttle)
library(ggplot2)
library(tidyr)
library(muscat)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
sim_prefix <- "augSimHierarchy"
leiden_res <- 0.1
results_subdir <- "results"
repo_root <- .repo_root_default
i <- 1
while (i <= length(args)) {
  if (args[i] == "--sim_prefix" && i < length(args)) { sim_prefix <- args[i + 1]; i <- i + 2; next }
  if (args[i] == "--leiden_res" && i < length(args)) { leiden_res <- as.numeric(args[i + 1]); i <- i + 2; next }
  if (args[i] == "--results_subdir" && i < length(args)) { results_subdir <- args[i + 1]; i <- i + 2; next }
  if (args[i] == "--repo-root" && i < length(args)) { repo_root <- args[i + 1]; i <- i + 2; next }
  i <- i + 1
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)
source(file.path(repo_root, "R", "ensure_scSampleSim.R"))
ensure_scSampleSim()

working_dir <- file.path(repo_root, "simulation_scripts")
raw_data_path <- file.path(working_dir, "data", "filtered_sce_data.Rda")
if (!file.exists(raw_data_path)) {
  stop("Input SCE not found: ", raw_data_path)
}
param_file_prefix <- "pa_preseurat_"

results_dir <- file.path(working_dir, results_subdir)
configs_dir <- file.path(results_dir, "de_configs")
sims_info_dir <- file.path(results_dir, "sims_info")
dir.create(configs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sims_info_dir, recursive = TRUE, showWarnings = FALSE)
target_base_seurat_fn <- file.path(configs_dir, "seurat_null_genes_augSimHierarchy.rds")

print(paste0("param_creator args: sim_prefix=", sim_prefix, " leiden_res=", leiden_res, " results_subdir=", results_subdir))

print("Creating params list for data augmentation")
params <- list()
params$id_int <- "kmeans_all_lfc_1"
params$datapath <- raw_data_path
params$dest_file <- target_base_seurat_fn
params$seed <- NULL
params$condition <- NULL
params$randomization <- "samples"
params$clustering <- "kmeans"
params$cut_at <- c(3, 10, 20)
params$use_harmony <- TRUE
params$n_de_genes <- 750
params$lfc_mean <- c(1.125, 2.25, 3.375)
params$working_dir <- working_dir
params$result_id <- "augSimHierarchy_base.Rda"
params$n_hvgs <- 5000
params$clust_to_sample_from <- NULL
params$leiden_resolution <- leiden_res

main_covariate <- "group"
sample_covariate <- "sample"
batch_covariate <- "batch"
assay_continuous <- "logcounts"
cell_type_column <- "cell_type"

print("Creating args list pa")
pa <- list(
  datapath = params$datapath,
  dest_file = params$dest_file,
  condition = params$condition,
  seed = params$seed,
  randomization = params$randomization,
  n_hvgs = params$n_hvgs,
  cut_at = params$cut_at,
  n_de_genes = params$n_de_genes,
  lfc_mean = params$lfc_mean,
  clustering = params$clustering,
  working_dir = params$working_dir,
  result_id = params$result_id,
  use_harmony = params$use_harmony,
  clust_to_sample_from = params$clust_to_sample_from,
  leiden_resolution = params$leiden_resolution,
  main_covariate = main_covariate,
  sample_covariate = sample_covariate,
  assay_continuous = assay_continuous,
  cell_type_column = cell_type_column,
  batch_covariate = batch_covariate
)

verbose <- TRUE
print("Preparing augmented data")
args_prep <- pa[names(pa) %in% formalArgs(scSampleSim::prep_aug_sim)]
args_prep$verbose <- verbose
pa <- do.call(scSampleSim::prep_aug_sim, args_prep)

print("Saving pa before processing for use in analysis")
saveRDS(pa, file.path(sims_info_dir, paste0(param_file_prefix, sim_prefix, ".Rds")))
print(paste0("Wrote ", file.path(sims_info_dir, paste0(param_file_prefix, sim_prefix, ".Rds"))))

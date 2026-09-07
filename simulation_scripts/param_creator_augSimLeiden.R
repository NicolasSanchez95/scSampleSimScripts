#' Example parameter setup for Leiden-based augmented simulations.
#'
#' Builds the base Seurat object and `pa_preseurat_<sim_prefix>.Rds` used by
#' `sim_analysis_pipeline.R`. Place the input SCE at
#' `simulation_scripts/data/filtered_sce_data.Rda` with colData columns
#' `group`, `sample`, `batch`, and `cell_type`.
#'
#' Usage:
#'   Rscript param_creator_augSimLeiden.R [--sim_prefix NAME] [--leiden_res 0.1]
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
sim_prefix <- "augSimLeiden"
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
source(file.path(repo_root, "R", "ensure_postselect.R"))
ensure_postselect()

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
target_base_seurat_fn <- file.path(configs_dir, "seurat_null_genes_augSimLeiden.rds")

print(paste0("param_creator args: sim_prefix=", sim_prefix, " leiden_res=", leiden_res, " results_subdir=", results_subdir))

print("Creating params list for data augmentation")
params <- list()
params$id_int <- "leiden_single_lfc_Q"
params$datapath <- raw_data_path
params$dest_file <- target_base_seurat_fn
params$seed <- NULL
params$condition <- NULL
params$randomization <- "samples"
params$clustering <- "leiden"
params$cut_at <- 1
params$use_harmony <- TRUE
params$n_de_genes <- 50 # will be overwritten by the number of DE genes specified in the simulation script
params$lfc_mean <- 3
params$working_dir <- working_dir
params$result_id <- "augSimLeiden_base.Rda"
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
args_prep <- pa[names(pa) %in% formalArgs(postselect::prep_aug_sim)]
args_prep$verbose <- verbose
pa <- do.call(postselect::prep_aug_sim, args_prep)

print("Saving pa before processing for use in analysis")
pa_before_processing <- pa
saveRDS(pa_before_processing, file.path(sims_info_dir, paste0(param_file_prefix, sim_prefix, ".Rds")))

new_id_check <- "param_setup_nsim0"
cluster_type <- "leiden"
num_pcs <- 50
num_sim_genes <- 0

print("Simulating 0 genes (oracle cluster labels on null genes)")
sim_data <- postselect::create_full_augmented_data_w_num_sim_genes(
  sim_prefix = sim_prefix,
  new_id_check = new_id_check,
  num_sim_genes = num_sim_genes,
  base_seurat_fn = target_base_seurat_fn,
  sims_info_dir = paste0(sims_info_dir, "/"),
  param_file_prefix = param_file_prefix,
  de_args_save_prfx = "de_args_",
  condition_save_prfx = "fake_condition_",
  sample_ids_save_prfx = "sample_ids_",
  gt_clusterings_save_prfx = "gt_clusterings_",
  verbose = TRUE,
  save_down = FALSE
)
used_sce <- sim_data$used_sce
pa_de <- sim_data$pa_de

print("Getting null genes")
null_genes <- setdiff(rownames(used_sce), pa_de$set_de_genes)
null_sce <- used_sce[null_genes, ]
print("Running PCA and Harmony clustering for null genes")
clust_results_null <- postselect::run_pca_harmony_leiden(
  null_sce,
  cluster_type = cluster_type,
  leiden_res = leiden_res,
  num_pcs = num_pcs,
  include_batch = pa_de$include_batch
)
pa_before_processing$leiden_clusters <- clust_results_null$cluster_assignment_list_harmony
saveRDS(pa_before_processing, file.path(sims_info_dir, paste0(param_file_prefix, sim_prefix, ".Rds")))
print(paste0("Wrote ", file.path(sims_info_dir, paste0(param_file_prefix, sim_prefix, ".Rds"))))

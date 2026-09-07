library(SingleCellExperiment)
library(Seurat)
library(harmony)
library(tibble)
library(dplyr)
library(stringr)
library(variancePartition)
library(limma)
library(pbapply)
library(edgeR)
library(scuttle)
library(SummarizedExperiment)

print("defining paths")
###################### DEFINE PATHS AND PREFIXES ######################

# define file prefixes for intermediate saves
de_args_save_prfx <- "de_args_"
condition_save_prfx <- "fake_condition_"
sample_ids_save_prfx <- "sample_ids_"
de_overlap_info_save_prfx <- "de_overlap_info_"
pa_de_save_prfx <- "pa_de_int_"
fdps_storage_save_prfx <- "control_store_all_"
res_de_save_prfx <- "res_de_"
clustering_prefix <- "cluster_assignment_"
var_fracs_save_prfx <- "var_fracs_de_"
cca_anls_save_prfx <- "cca_anls_de_"
PVE_metrics_save_prfx <- "PVE_metrics_"
mv_PVE_metrics_save_prfx <- "mv_PVE_embed_"
param_file_prefix <-  "pa_preseurat_"
file_naming_utils_prfx <- "file_naming_utils"
gt_clusterings_save_prfx <- "gt_clusterings_"


# Define the analysis patterns and control types
anls_patterns <- list(
  "no_harmony" = "noharm_",
  "null_noharm" = "noharm_null_",
  "null_harmony" = "harmony_null_",
  "null_fmnn" = "fmnn_null_",
  "full_harmony" = "harmony_",
  "TvsS" = "TvsS_",
  "fmnn" = "fmnn_",
  "oracle" = "oracle_"
)

sim_prefixes_w_oracle <- c("augSimLeiden", "augSimLeidenNoHarm", "augSimLeidenLow")

ctrl_types <- c("ctrl_clusts_loc", "ctrl_glob", "ctrl_clusts_glob")
ctrl_type <- "ctrl_glob"

###################### PARSE COMMAND LINE ARGUMENTS ######################
args <- commandArgs(trailingOnly = TRUE)
parser <- argparse::ArgumentParser()
parser$add_argument("--repo-root", type="character", default="",
                    help="Root of the paper-facing scripts repo (default: parent of this script)")
parser$add_argument("--results_subdir", type="character", default="results",
                    help="Subdirectory under simulation_scripts for simulation outputs")
parser$add_argument("--sim_prefix", type="character", help="Sim prefix")
parser$add_argument("--num_de_genes", type="integer", help="Number of de genes",default=0)
parser$add_argument("--sim_type", type="character", help="Sim type")
parser$add_argument("--cluster_type", type="character", help="Cluster type")
parser$add_argument("--leiden_res", type="numeric", help="Leiden resolution",default=0.1)
parser$add_argument("--kmeans_k", type="numeric", help="Kmeans k",default=-1)
parser$add_argument("--id_check", type="integer", help="ID check number")
parser$add_argument("--lfc", type="numeric", help="LFC")
parser$add_argument("--pde", type="numeric", help="PDE")
parser$add_argument("--augDataRep", type="character", help="AugData rep", default="augSim31020")
parser$add_argument("--num_genes_for_TvsS_clustering", type="integer", help="Number of genes for TvsS clustering", default=2500)
parser$add_argument("--use_pseudobulk_TvsS", type="logical", help="Use pseudobulk for TvsS clustering", default=TRUE)
parser$add_argument("--sig_threshold", type="numeric", help="Significance threshold", default=0.1)
parser$add_argument("--num_pcs", type="integer", help="Number of PCs to use for PCA embedding", default=50)
parser$add_argument("--run_PVE_metrics", type="logical", help="Run PVE metrics", default=FALSE)
parser$add_argument(
  "--anls_patterns",
  type = "character",
  default = "",
  help = "Comma-separated list of analysis pattern names to run (subset of: no_harmony,null_noharm,null_harmony,null_fmnn,full_harmony,TvsS,fmnn,oracle). Empty means run all."
)


args <- parser$parse_args(args)
.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .script_args, value = TRUE)
.this_script <- if (length(.file_arg) >= 1L) {
  normalizePath(sub("^--file=", "", .file_arg[[1]]), mustWork = TRUE)
} else {
  stop("Run this file with Rscript so --file= is set.")
}
repo_root <- args$repo_root
if (is.null(repo_root) || !nzchar(repo_root)) {
  repo_root <- normalizePath(file.path(dirname(.this_script), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
}
source(file.path(repo_root, "R", "ensure_postselect.R"))
ensure_postselect()
working_dir <- file.path(repo_root, "simulation_scripts")
raw_data_path <- file.path(working_dir, "data", "filtered_sce_data.Rda")

results_subdir <- args$results_subdir
sim_prefix <- args$sim_prefix
num_de_genes <- args$num_de_genes
id_check <- args$id_check
lfc <- args$lfc
pde <- args$pde
sim_type <- args$sim_type
cluster_type <- args$cluster_type
leiden_res <- args$leiden_res
kmeans_k <- args$kmeans_k
sig_threshold <- args$sig_threshold
num_pcs <- args$num_pcs
augDataRep <- args$augDataRep
num_genes_for_TvsS_clustering <- args$num_genes_for_TvsS_clustering
use_pseudobulk_TvsS <- args$use_pseudobulk_TvsS
run_PVE_metrics <- args$run_PVE_metrics
anls_patterns_arg <- args$anls_patterns

if (is.null(sim_prefix)) {
  stop("Invalid sim_prefix")
}
if (is.null(id_check)) {
  stop("Invalid id_check")
}

# Decide which analyses to run
valid_anls <- names(anls_patterns)
selected_anls <- valid_anls
if (!is.null(anls_patterns_arg) && nzchar(trimws(anls_patterns_arg))) {
  selected_anls <- unlist(strsplit(anls_patterns_arg, ",", fixed = TRUE))
  selected_anls <- trimws(selected_anls)
  selected_anls <- selected_anls[nzchar(selected_anls)]
  unknown <- setdiff(selected_anls, valid_anls)
  if (length(unknown) > 0) {
    stop(
      paste0(
        "Unknown --anls_patterns entries: ",
        paste(unknown, collapse = ", "),
        ". Valid: ",
        paste(valid_anls, collapse = ", ")
      )
    )
  }
}

# define paths to folders for intermediate and final saves (after --results_subdir)
results_dir <- file.path(working_dir, results_subdir)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
configs_dir <- paste0(results_dir, "/de_configs/")
if (!dir.exists(configs_dir)) {
  dir.create(configs_dir, recursive = TRUE)
}
sims_info_dir <- paste0(results_dir, "/sims_info/")
if (!dir.exists(sims_info_dir)) {
  dir.create(sims_info_dir, recursive = TRUE)
}
anls_info_dir <- paste0(results_dir, "/anls_info/")
if (!dir.exists(anls_info_dir)) {
  dir.create(anls_info_dir, recursive = TRUE)
}
des_dir <- paste0(results_dir, "/de_final/")
if (!dir.exists(des_dir)) {
  dir.create(des_dir, recursive = TRUE)
}

file_naming_utils <- list(
  results_dir = results_dir,
  sims_info_dir = sims_info_dir,
  anls_info_dir = anls_info_dir,
  raw_data_path = raw_data_path,
  configs_dir = configs_dir,
  des_dir = des_dir,
  param_file_prefix = param_file_prefix,
  de_args_save_prfx = de_args_save_prfx,
  condition_save_prfx = condition_save_prfx,
  sample_ids_save_prfx = sample_ids_save_prfx,
  de_overlap_info_save_prfx = de_overlap_info_save_prfx,
  pa_de_save_prfx = pa_de_save_prfx,
  fdps_storage_save_prfx = fdps_storage_save_prfx,
  res_de_save_prfx = res_de_save_prfx,
  clustering_prefix = clustering_prefix,
  var_fracs_save_prfx = var_fracs_save_prfx,
  cca_anls_save_prfx = cca_anls_save_prfx,
  PVE_metrics_save_prfx = PVE_metrics_save_prfx,
  mv_PVE_metrics_save_prfx = mv_PVE_metrics_save_prfx,
  anls_patterns = anls_patterns,
  ctrl_types = ctrl_types,
  file_naming_utils_prfx = file_naming_utils_prfx,
  gt_clusterings_save_prfx = gt_clusterings_save_prfx
)

pa_rds_path <- file.path(sims_info_dir, paste0(param_file_prefix, sim_prefix, ".Rds"))
if (!file.exists(pa_rds_path)) {
  stop(paste0("Missing param RDS (run param creator first): ", pa_rds_path))
}
pa_pre <- readRDS(pa_rds_path)
if (is.null(pa_pre$dest_file) || !nzchar(as.character(pa_pre$dest_file))) {
  stop(paste0("Param RDS has no dest_file: ", pa_rds_path))
}
base_seurat_fn <- as.character(pa_pre$dest_file)
if (!file.exists(base_seurat_fn)) {
  stop(paste0("Base Seurat not found at dest_file path: ", base_seurat_fn))
}
file_naming_utils$base_seurat_fn <- base_seurat_fn

saveRDS(file_naming_utils, paste0(anls_info_dir, "/", file_naming_utils_prfx, ".rds"))

###################### DEFINE DISPLAY CONTROLPARAMETERS ######################
cut_off_true <- 0.00
cut_off_false <- 0.00
overlap_type <- "union"
sig_threshold <- 0.1

print(paste0("Running analysis with sim_type: ", sim_type, " and sim_prefix: ", sim_prefix, " and id_check: ", id_check))
#assert that sim_type is either mSim or augData



# if no arguments are provided, use default values for debuggin
# if(is.null(args$sim_type)){
#   sim_type <- "augData"
#   sim_prefix <- "augSimHierarchy" 
#   num_de_genes <- 500 # 500# augData "mSim"
#   cluster_type <- "leiden" # leiden"kmeans" 
#   id_check <- 11111111 # 90710
# }

new_id_check <- paste0(sim_prefix, "_", id_check, "_nsim_", num_de_genes)
last_anls <- if (length(selected_anls) > 0) tail(selected_anls, 1) else "fmnn"
pve_results_fn <- file.path(file_naming_utils$anls_info_dir, paste0(file_naming_utils$PVE_metrics_save_prfx, file_naming_utils$anls_patterns[[last_anls]], new_id_check, ".rds"))

print(file.exists(pve_results_fn))

if (file.exists(pve_results_fn)){
  print(paste0("PVE results already exist for ", new_id_check))
  pve_results <- readRDS(pve_results_fn)
  print(paste0("PVE results: ", pve_results))
  stop("PVE results already exist, skipping analysis")
}

print(paste0("id_check: ", id_check))
print(paste0("num_de_genes: ", num_de_genes))
print(paste0("lfc: ", lfc))
print(paste0("sim_type: ", sim_type))
print(paste0("sim_prefix: ", sim_prefix))
print(paste0("cluster_type: ", cluster_type))
print(paste0("leiden_res: ", leiden_res))
print(paste0("kmeans_k: ", kmeans_k))
print(paste0("augDataRep: ", augDataRep))

print(paste0("creating the simulated data according to the sim_type: ", sim_type))
start_time <- Sys.time()
args_sim <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::create_simulated_data)]
args_sim$sim_type <- sim_type
args_sim$sim_prefix <- sim_prefix
args_sim$new_id_check <- new_id_check
args_sim$num_de_genes <- num_de_genes
sim_data <- do.call(postselect::create_simulated_data, args_sim)
used_sce <- sim_data$used_sce
end_time <- Sys.time()
print(paste0("Time taken to create the simulated data: ", end_time - start_time))

print("saving parameters relevant to the DE analysis")
pa_de <- sim_data$pa_de
pa_de$cluster_type <- cluster_type
pa_de$kmeans_k <- kmeans_k
pa_de$leiden_res <- leiden_res
pa_de$sig_threshold <- sig_threshold
pa_de$cut_off_true <- cut_off_true
pa_de$cut_off_false <- cut_off_false
pa_de$sim_type <- sim_type
pa_de$num_pcs <- num_pcs
pa_de$num_genes_for_TvsS_clustering <- num_genes_for_TvsS_clustering
pa_de$use_pseudobulk_TvsS <- use_pseudobulk_TvsS

new_id_check <- pa_de$new_id_check
saveRDS(pa_de, file.path(file_naming_utils$anls_info_dir, paste0(file_naming_utils$pa_de_save_prfx, new_id_check, ".rds")))


############### PCA AND HARMONY CLUSTERING ON ALL GENES ###############
start_time <- Sys.time()
print(paste0("Running PCA and Harmony clustering according to the clustering type: ", cluster_type))
clust_results <- postselect::run_pca_harmony_leiden(used_sce, cluster_type = cluster_type, leiden_res = leiden_res,num_pcs = num_pcs, include_batch = pa_de$include_batch)
end_time <- Sys.time()

############### ORACLE CLUSTERING ###############
if ("oracle" %in% selected_anls && sim_prefix %in% sim_prefixes_w_oracle){
  print("Running oracle analysis")
  start_time <- Sys.time()
  curr_anls <- "oracle"
  cluster_assignment_curr_anls <- pa_de$leiden_clusters
  lblnorm_counts <- clust_results$lblnorm_counts
  args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
  args_anls$clustering_assignment_curr_anls <- cluster_assignment_curr_anls
  args_anls$lblnorm_counts <- lblnorm_counts
  args_anls$used_sce <- used_sce
  args_anls$pa_de <- pa_de
  args_anls$sim_type <- sim_type
  args_anls$curr_anls <- curr_anls
  args_anls$new_id_check <- new_id_check
  args_anls$cut_off_true <- cut_off_true
  args_anls$cut_off_false <- cut_off_false
  args_anls$sig_threshold <- sig_threshold
  args_anls$overlap_type <- overlap_type
  args_anls$run_PVE <- run_PVE_metrics
  args_anls$cell_embeddings <- clust_results$pca_embeds
  args_anls$embedding_pc_weights <- (clust_results$pca_stdev)^2
  anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
  end_time <- Sys.time()
  print(paste0("Time taken to run oracle analysis: ", end_time - start_time))
  postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
}




print(paste0("Time taken to run PCA and Harmony clustering: ", end_time - start_time))
if ("no_harmony" %in% selected_anls) {
  print("Running analysis on all of the data with no harmony")
  start_time <- Sys.time()
  curr_anls <- "no_harmony"
  cluster_label_name <- "cluster_assignment_list"
  embedding_name <- "pca_embeds"
  cluster_assignment_curr_anls <- clust_results$cluster_assignment_list
  lblnorm_counts <- clust_results$lblnorm_counts
  embedding_name <- "pca_embeds"
  args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
  args_anls$clustering_assignment_curr_anls <- cluster_assignment_curr_anls
  args_anls$lblnorm_counts <- lblnorm_counts
  args_anls$used_sce <- used_sce
  args_anls$pa_de <- pa_de
  args_anls$sim_type <- sim_type
  args_anls$curr_anls <- curr_anls
  args_anls$new_id_check <- new_id_check
  args_anls$cut_off_true <- cut_off_true
  args_anls$cut_off_false <- cut_off_false
  args_anls$sig_threshold <- sig_threshold
  args_anls$overlap_type <- overlap_type
  args_anls$run_PVE <- run_PVE_metrics
  args_anls$cell_embeddings <- clust_results$pca_embeds
  args_anls$embedding_pc_weights <- (clust_results$pca_stdev)^2
  anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
  end_time <- Sys.time()
  print(paste0("Time taken to run analysis on all of the data with no harmony: ", end_time - start_time))
  postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
}

if ("full_harmony" %in% selected_anls) {
  print("Running analysis on all of the data with harmony")
  start_time <- Sys.time()
  curr_anls <- "full_harmony"
  cluster_assignment_curr_anls <- clust_results$cluster_assignment_list_harmony
  args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
  args_anls$clustering_assignment_curr_anls <- cluster_assignment_curr_anls
  args_anls$lblnorm_counts <- lblnorm_counts
  args_anls$used_sce <- used_sce
  args_anls$pa_de <- pa_de
  args_anls$sim_type <- sim_type
  args_anls$curr_anls <- curr_anls
  args_anls$new_id_check <- new_id_check
  args_anls$cut_off_true <- cut_off_true
  args_anls$cut_off_false <- cut_off_false
  args_anls$sig_threshold <- sig_threshold
  args_anls$overlap_type <- overlap_type
  args_anls$run_PVE <- run_PVE_metrics
  args_anls$cell_embeddings <- clust_results$harmony_embeddings
  args_anls$embedding_pc_weights <- NULL
  anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
  end_time <- Sys.time()
  print(paste0("Time taken to run analysis on all of the data with harmony: ", end_time - start_time))
  postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
}



############### PCA AND HARMONY CLUSTERING ON NULL GENES ###############
if (any(c("null_noharm", "null_harmony", "null_fmnn") %in% selected_anls)) {
  print("Running analysis on null genes (subset requested)")
  start_time <- Sys.time()
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
  lblnorm_counts <- clust_results_null$lblnorm_counts

  if ("null_noharm" %in% selected_anls) {
    print("Running analysis on null genes without harmony")
    start_time <- Sys.time()
    curr_anls <- "null_noharm"
    cluster_assignment_curr_anls <- clust_results_null$cluster_assignment_list
    embedding_name <- "pca_embeds"
    args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
    args_anls$clustering_assignment_curr_anls <- cluster_assignment_curr_anls
    args_anls$lblnorm_counts <- lblnorm_counts
    args_anls$used_sce <- used_sce
    args_anls$pa_de <- pa_de
    args_anls$sim_type <- sim_type
    args_anls$curr_anls <- curr_anls
    args_anls$new_id_check <- new_id_check
    args_anls$cut_off_true <- cut_off_true
    args_anls$cut_off_false <- cut_off_false
    args_anls$sig_threshold <- sig_threshold
    args_anls$overlap_type <- overlap_type
    args_anls$run_PVE <- run_PVE_metrics
    args_anls$cell_embeddings <- clust_results_null$pca_embeds
    args_anls$embedding_pc_weights <- (clust_results_null$pca_stdev)^2
    anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
    end_time <- Sys.time()
    print(paste0("Time taken to run analysis on null genes without harmony: ", end_time - start_time))
    postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
  }

  if ("null_harmony" %in% selected_anls) {
    print("Running analysis on null genes with harmony")
    start_time <- Sys.time()
    curr_anls <- "null_harmony"
    cluster_assignment_curr_anls <- clust_results_null$cluster_assignment_list_harmony
    lblnorm_counts <- clust_results_null$lblnorm_counts
    args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
    args_anls$clustering_assignment_curr_anls <- cluster_assignment_curr_anls
    args_anls$lblnorm_counts <- lblnorm_counts
    args_anls$used_sce <- used_sce
    args_anls$pa_de <- pa_de
    args_anls$sim_type <- sim_type
    args_anls$curr_anls <- curr_anls
    args_anls$new_id_check <- new_id_check
    args_anls$cut_off_true <- cut_off_true
    args_anls$cut_off_false <- cut_off_false
    args_anls$sig_threshold <- sig_threshold
    args_anls$overlap_type <- overlap_type
    args_anls$run_PVE <- run_PVE_metrics
    args_anls$cell_embeddings <- clust_results_null$harmony_embeddings
    args_anls$embedding_pc_weights <- NULL
    anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
    end_time <- Sys.time()
    postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
  }

  # run fastMNN analysis on null genes
  if ("null_fmnn" %in% selected_anls) {
    print("Running fastMNN analysis on null genes")
    start_time <- Sys.time()
    curr_anls <- "null_fmnn"
    assay(null_sce, "logcounts", withDimnames = FALSE) <- clust_results_null$lblnorm_counts
    fastMNN_corrected <- batchelor::fastMNN(null_sce, batch = null_sce$sample_id, BSPARAM = BiocSingular::RandomParam())
    assay(fastMNN_corrected, "counts", withDimnames = FALSE) <- assay(null_sce, "counts")
    assay(fastMNN_corrected, "logcounts", withDimnames = FALSE) <- assay(null_sce, "logcounts")
    seurat_obj_fastMNN <- as.Seurat(fastMNN_corrected, counts = "counts", data = "logcounts")
    if (cluster_type == "leiden") {
      cluster_assignment_list_fastMNN <- postselect::cluster_using_leiden(seurat_obj_fastMNN, leiden_res = leiden_res, reduction_type = "corrected")
    } else if (cluster_type == "kmeans") {
      cluster_assignment_list_fastMNN <- postselect::cluster_using_kmeans(seurat_obj_fastMNN, kmeans_k = kmeans_k, reduction_type = "corrected")
    }
    args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
    args_anls$clustering_assignment_curr_anls <- cluster_assignment_list_fastMNN
    args_anls$lblnorm_counts <- clust_results_null$lblnorm_counts
    args_anls$used_sce <- used_sce
    args_anls$pa_de <- pa_de
    args_anls$sim_type <- sim_type
    args_anls$curr_anls <- curr_anls
    args_anls$new_id_check <- new_id_check
    args_anls$cut_off_true <- cut_off_true
    args_anls$cut_off_false <- cut_off_false
    args_anls$sig_threshold <- sig_threshold
    args_anls$overlap_type <- overlap_type
    args_anls$run_PVE <- run_PVE_metrics
    args_anls$cell_embeddings <- Seurat::Embeddings(seurat_obj_fastMNN, reduction = "corrected")
    args_anls$embedding_pc_weights <- NULL
    anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
    end_time <- Sys.time()
    print(paste0("Time taken to run fastMNN analysis on null genes: ", end_time - start_time))
    postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
  }
}



############### TvsS CLUSTERING ###############
if ("TvsS" %in% selected_anls) {
  print("Running TvsS analysis")
  start_time <- Sys.time()
  curr_anls <- "TvsS"
  args_tvs <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::find_TvsS_genes)]
  args_tvs$used_sce <- used_sce
  args_tvs$clust_results <- clust_results
  args_tvs$num_genes_for_TvsS_clustering <- num_genes_for_TvsS_clustering
  args_tvs$new_id_check <- new_id_check
  args_tvs$use_pseudobulk_TvsS <- use_pseudobulk_TvsS
  TvsS_genes <- do.call(postselect::find_TvsS_genes, args_tvs)
  used_sce_TvsS <- used_sce[TvsS_genes, ]
  clust_results_TvsS <- postselect::run_pca_harmony_leiden(
    used_sce_TvsS,
    cluster_type = cluster_type,
    leiden_res = leiden_res,
    num_pcs = num_pcs,
    include_batch = pa_de$include_batch
  )

  cluster_assignment_curr_anls <- clust_results_TvsS$cluster_assignment_list

  args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
  args_anls$clustering_assignment_curr_anls <- cluster_assignment_curr_anls
  args_anls$lblnorm_counts <- lblnorm_counts
  args_anls$used_sce <- used_sce
  args_anls$pa_de <- pa_de
  args_anls$sim_type <- sim_type
  args_anls$curr_anls <- curr_anls
  args_anls$new_id_check <- new_id_check
  args_anls$cut_off_true <- cut_off_true
  args_anls$cut_off_false <- cut_off_false
  args_anls$sig_threshold <- sig_threshold
  args_anls$overlap_type <- overlap_type
  args_anls$run_PVE <- run_PVE_metrics
  args_anls$cell_embeddings <- clust_results_TvsS$pca_embeds
  args_anls$embedding_pc_weights <- (clust_results_TvsS$pca_stdev)^2
  anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
  end_time <- Sys.time()
  print(paste0("Time taken to run TvsS analysis: ", end_time - start_time))
  postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
}



############### FMNN CLUSTERING ###############
if ("fmnn" %in% selected_anls) {
  print("Running FMNN analysis")
  curr_anls <- "fmnn"
  start_time <- Sys.time()
  assay(used_sce, "logcounts", withDimnames = FALSE) <- clust_results$lblnorm_counts
  fmnn_corrected <- batchelor::fastMNN(used_sce, batch = used_sce$sample_id, BSPARAM = BiocSingular::RandomParam())
  assay(fmnn_corrected, "counts", withDimnames = FALSE) <- assay(used_sce, "counts")
  assay(fmnn_corrected, "logcounts", withDimnames = FALSE) <- assay(used_sce, "logcounts")
  seurat_obj_fmnn <- as.Seurat(fmnn_corrected, counts = "counts", data = "logcounts")

  if (cluster_type == "leiden") {
    cluster_assignment_list_fmnn <- postselect::cluster_using_leiden(seurat_obj_fmnn, leiden_res = leiden_res, reduction_type = "corrected")
  } else if (cluster_type == "kmeans") {
    cluster_assignment_list_fmnn <- postselect::cluster_using_kmeans(seurat_obj_fmnn, kmeans_k = kmeans_k, reduction_type = "corrected")
  }

  args_anls <- file_naming_utils[names(file_naming_utils) %in% formalArgs(postselect::run_analysis_for_clustering)]
  args_anls$clustering_assignment_curr_anls <- cluster_assignment_list_fmnn
  args_anls$lblnorm_counts <- clust_results$lblnorm_counts
  args_anls$used_sce <- used_sce
  args_anls$pa_de <- pa_de
  args_anls$sim_type <- sim_type
  args_anls$curr_anls <- curr_anls
  args_anls$new_id_check <- new_id_check
  args_anls$cut_off_true <- cut_off_true
  args_anls$cut_off_false <- cut_off_false
  args_anls$sig_threshold <- sig_threshold
  args_anls$overlap_type <- overlap_type
  args_anls$run_PVE <- run_PVE_metrics
  args_anls$cell_embeddings <- Seurat::Embeddings(seurat_obj_fmnn, reduction = "corrected")
  args_anls$embedding_pc_weights <- NULL
  anls_res <- do.call(postselect::run_analysis_for_clustering, args_anls)
  end_time <- Sys.time()
  print(paste0("Time taken to run FMNN analysis: ", end_time - start_time))
  postselect::print_quick_glob_ctrl(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type)
}




print("All analyses completed, explicitly removing all objects from memory")
rm(list = ls())
gc()


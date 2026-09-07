# Install and load the pinned postselectPaper build (tested backup tag).
# Sourced by every R entrypoint. Keep POSTSELECT_REF / POSTSELECT_SHA in sync with config.sh.

POSTSELECT_REPO <- "epurdom/postselectPaper"
POSTSELECT_REF <- "tested-backup-20260904"
POSTSELECT_SHA <- "ba615262f3447073e8e67018e95d5240bfe2893e"

#' Directory containing this file's parent (the paper-facing repo root).
find_repo_root <- function(cli_repo_root = NULL) {
  if (!is.null(cli_repo_root) && nzchar(cli_repo_root)) {
    return(normalizePath(cli_repo_root, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) >= 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    d <- dirname(script_path)
    if (basename(d) %in% c("simulation_scripts", "master_table_scripts", "decision_parameters", "R")) {
      return(normalizePath(dirname(d), mustWork = TRUE))
    }
    return(d)
  }
  normalizePath(getwd(), mustWork = FALSE)
}

# Bioconductor Imports of postselect that remotes will not find on CRAN alone.
POSTSELECT_BIOC_IMPORTS <- c(
  "batchelor", "BiocSingular", "BiocNeighbors", "bluster",
  "SingleCellExperiment", "SummarizedExperiment", "MatrixGenerics",
  "S4Vectors", "edgeR", "limma", "muscat", "scuttle", "lemur"
)

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

ensure_postselect <- function() {
  need_install <- TRUE
  if (requireNamespace("postselect", quietly = TRUE)) {
    sha <- utils::packageDescription("postselect")$RemoteSha
    if (!is.null(sha) && startsWith(tolower(sha), substr(POSTSELECT_SHA, 1L, 7L))) {
      need_install <- FALSE
    }
  }
  if (need_install) {
    ensure_pkg("remotes")
    ensure_pkg("BiocManager")
    missing_bioc <- POSTSELECT_BIOC_IMPORTS[
      !vapply(POSTSELECT_BIOC_IMPORTS, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_bioc) > 0L) {
      message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
      BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
    }
    message(
      "Installing postselect from GitHub ", POSTSELECT_REPO,
      " @ ", POSTSELECT_REF, " (", substr(POSTSELECT_SHA, 1L, 7L), ")"
    )
    remotes::install_github(
      POSTSELECT_REPO,
      ref = POSTSELECT_REF,
      upgrade = "never",
      repos = BiocManager::repositories()
    )
  }
  library(postselect)
  invisible(TRUE)
}

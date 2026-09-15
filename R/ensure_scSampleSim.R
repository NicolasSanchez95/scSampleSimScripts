# Install and load the pinned scSampleSim build.
# Sourced by every R entrypoint. Keep SCSAMPLESIM_REF / SCSAMPLESIM_SHA in sync with config.sh.
#
# Pin must be a commit whose DESCRIPTION says Package: scSampleSim (not the
# old tested-backup-20260904 / ba61526 postselect tag). Update REF and SHA
# after that rename is committed and tagged on GitHub.

SCSAMPLESIM_REPO <- "epurdom/scSampleSim"
SCSAMPLESIM_REF <- "scSampleSim-20260910"
SCSAMPLESIM_SHA <- "PENDING_AFTER_SCSAMPLESIM_RENAME_TAG"

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

# Bioconductor Imports of scSampleSim that remotes will not find on CRAN alone.
SCSAMPLESIM_BIOC_IMPORTS <- c(
  "batchelor", "BiocSingular", "BiocNeighbors", "bluster",
  "SingleCellExperiment", "SummarizedExperiment", "MatrixGenerics",
  "S4Vectors", "edgeR", "limma", "muscat", "scuttle", "lemur"
)

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# BiocManager::install(ask = FALSE) will not create a personal library. On
# clusters the site library is often not writable, so make sure one exists.
ensure_writable_lib <- function() {
  can_write <- function(p) {
    nzchar(p) && dir.exists(p) && file.access(p, 2L) == 0L
  }
  if (any(vapply(.libPaths(), can_write, logical(1)))) {
    return(invisible(NULL))
  }
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib) || identical(user_lib, "NULL")) {
    majmin <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
    user_lib <- file.path(
      path.expand("~"), "R", paste0(R.version$platform, "-library"), majmin
    )
  } else {
    user_lib <- path.expand(user_lib)
  }
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  message("Using user library: ", user_lib)
  invisible(NULL)
}

ensure_scSampleSim <- function() {
  ensure_writable_lib()
  need_install <- TRUE
  if (requireNamespace("scSampleSim", quietly = TRUE)) {
    sha <- utils::packageDescription("scSampleSim")$RemoteSha
    if (!is.null(sha) && startsWith(tolower(sha), substr(SCSAMPLESIM_SHA, 1L, 7L)) &&
        nchar(SCSAMPLESIM_SHA) >= 7L && !startsWith(SCSAMPLESIM_SHA, "PENDING")) {
      need_install <- FALSE
    }
  }
  if (need_install) {
    ensure_pkg("remotes")
    ensure_pkg("BiocManager")
    missing_bioc <- SCSAMPLESIM_BIOC_IMPORTS[
      !vapply(SCSAMPLESIM_BIOC_IMPORTS, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_bioc) > 0L) {
      message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
      BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
    }
    message(
      "Installing scSampleSim from GitHub ", SCSAMPLESIM_REPO,
      " @ ", SCSAMPLESIM_REF, " (", substr(SCSAMPLESIM_SHA, 1L, 7L), ")"
    )
    remotes::install_github(
      SCSAMPLESIM_REPO,
      ref = SCSAMPLESIM_REF,
      upgrade = "never",
      repos = BiocManager::repositories()
    )
  }
  library(scSampleSim)
  invisible(TRUE)
}

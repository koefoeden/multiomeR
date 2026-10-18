get_project_root <- function() {
  rprojroot::find_root(rprojroot::has_file("pixi.toml"))
}

#' Stop when the BSgenome data packages are missing
#'
#' Pixi installs these Bioconda data packages only when post-link scripts run.
#' `_targets.R` calls this so a pipeline run names the install command instead
#' of failing at the first target that reads a genome sequence; sessions that
#' only test, document or install packages do not need the genomes.
#'
#' @return Invisibly returns `TRUE`.
#' @keywords internal
assert_genome_packages_installed <- function() {
  genome_packages <- c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Mmusculus.UCSC.mm39")
  missing_genome_packages <- Filter(\(package) !nzchar(system.file(package = package)), genome_packages)
  if (length(missing_genome_packages)) {
    stop(
      paste(missing_genome_packages, collapse = ", "),
      " is not installed. Run `pixi install --locked --run-post-link-scripts` from the repository root.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop when R's build configuration is empty
#'
#' r-base's Conda activation script reruns `R CMD javareconf` at every uncached
#' Pixi activation, and a failed run can leave `Makeconf` and `ldpaths` empty.
#' R still starts, but it cannot compile the native kernels that targets build
#' at run time, so stop at startup with the repair instead.
#'
#' @return Invisibly returns `TRUE`.
#' @keywords internal
assert_R_build_configuration <- function() {
  config_files <- file.path(R.home("etc"), c("Makeconf", "ldpaths"))
  empty_files <- config_files[file.size(config_files) %in% 0]
  if (length(empty_files)) {
    stop(
      paste(empty_files, collapse = ", "),
      " is empty, so R cannot compile native code. Run `pixi reinstall r-base` from the repository root, ",
      "and pass `--use-environment-activation-cache` to `pixi run` so activation scripts do not run again.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Load the project runtime
#'
#' Attach the core packages and conflict preferences, source the project
#' helpers, and install the targets options and crew controllers. `.Rprofile`
#' calls this at startup; call it again to reload edited helpers or controllers
#' in a running session.
#'
#' @return Invisibly returns `TRUE`.
#' @keywords internal
load_project_runtime <- function() {
  assert_R_build_configuration()
  suppressPackageStartupMessages({
    library(Matrix)
    library(purrr)
    library(tibble)
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(stringr)
    library(dplyr)
    library(assertthat)
    library(magrittr)
    library(targets)
    library(tarchetypes)
    library(conflicted)
  })

  conflicted::conflicts_prefer(
    dplyr::filter,
    dplyr::select,
    dplyr::slice,
    dplyr::desc,
    dplyr::rename,
    purrr::set_names,
    purrr::reduce,
    magrittr::extract,
    magrittr::subtract,
    tidyr::expand,
    base::intersect,
    base::setdiff,
    base::unname,
    MatrixGenerics::rowMedians,
    base::as.factor,
    .quiet = TRUE
  )

  targets::tar_source(c("packages/multiomeRCore/R", "R"))

  ggplot2::theme_set(ggplot2::theme_bw())
  ggplot2::theme_update(legend.position = "bottom")
  Sys.setenv("R_MSG_PKG_START_MSG" = "FALSE")
  # Seurat's C++ kernels run single-threaded unless told otherwise; give them the
  # cores of the worker's Slurm job.
  options(Seurat.nthreads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  targets::tar_option_set(
    error = "trim",
    iteration = "list",
    format = "qs",
    # A full gc() before every target costs about 0.2 s per task; R still collects garbage as it needs to.
    garbage_collection = 0L
  )

  apply_crew_controller_options(
    source(configuration_path("crew_controllers.R"))$value
  )
  invisible(TRUE)
}

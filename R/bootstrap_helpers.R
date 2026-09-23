get_project_root <- function() {
  rprojroot::find_root(rprojroot::has_file("pixi.toml"))
}

#' Load the project runtime
#'
#' Attach the core packages and conflict preferences, source the project
#' helpers, and install the targets options and crew controllers. `.Rprofile` calls this at startup; call it again to
#' reload edited helpers or controllers in a running session.
#'
#' @return Invisibly returns `TRUE`.
#' @keywords internal
load_project_runtime <- function() {
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
  targets::tar_option_set(
    error = "trim",
    iteration = "list",
    format = "qs",
    garbage_collection = 1L
  )

  apply_crew_controller_options(
    source(file.path(get_project_root(), "crew_controllers.R"), chdir = TRUE)$value
  )
  invisible(TRUE)
}

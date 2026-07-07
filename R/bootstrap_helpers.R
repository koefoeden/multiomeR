if (!exists("bootstrap_state_env", inherits = FALSE)) {
  bootstrap_state_env <- new.env(parent = emptyenv())
}

if (!exists("helper_packages_loaded", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$helper_packages_loaded <- FALSE
}

if (!exists("shared_helpers_loaded", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$shared_helpers_loaded <- FALSE
}

if (!exists("runtime_options_applied", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$runtime_options_applied <- FALSE
}

if (!exists("targets_patches_applied", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$targets_patches_applied <- FALSE
}

if (!exists("controllers_loaded", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$controllers_loaded <- FALSE
}

if (!exists("controller_resources_tibble", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$controller_resources_tibble <- NULL
}

if (!exists("project_root", envir = bootstrap_state_env, inherits = FALSE)) {
  bootstrap_state_env$project_root <- NA_character_
}

get_project_root <- function(force = FALSE) {
  if (!force && !is.na(bootstrap_state_env$project_root)) {
    return(bootstrap_state_env$project_root)
  }

  bootstrap_state_env$project_root <- rprojroot::find_root(
    criterion = rprojroot::has_file("pixi.toml"),
    path = base::getwd()
  )
  bootstrap_state_env$project_root
}

#' Load project packages
#'
#' Load core packages and install conflict preferences for the targets workflow.
#'
#' @param force Logical; TRUE reruns the setup step even when bootstrap state says it has already completed.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

load_project_packages <- function(force = FALSE) {
  if (isTRUE(bootstrap_state_env$helper_packages_loaded) && !force) {
    return(invisible(FALSE))
  }

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

  bootstrap_state_env$helper_packages_loaded <- TRUE
  invisible(TRUE)
}

source_shared_helpers <- function(force = FALSE) {
  if (isTRUE(bootstrap_state_env$shared_helpers_loaded) && !force) {
    return(invisible(FALSE))
  }

  load_project_packages()
  targets::tar_source("R")
  bootstrap_state_env$shared_helpers_loaded <- TRUE
  invisible(TRUE)
}

apply_runtime_options <- function(force = FALSE) {
  if (isTRUE(bootstrap_state_env$runtime_options_applied) && !force) {
    return(invisible(FALSE))
  }

  ggplot2::theme_set(ggplot2::theme_bw())
  ggplot2::theme_update(legend.position = "bottom")
  Sys.setenv("R_MSG_PKG_START_MSG" = "FALSE")

  if (is.null(getOption("multiomeR.save_serialized_plot_objects"))) {
    options(multiomeR.save_serialized_plot_objects = TRUE)
  }

  targets::tar_option_set(
    error = "trim",
    iteration = "list",
    format = "qs"
  )

  bootstrap_state_env$runtime_options_applied <- TRUE
  invisible(TRUE)
}

#' Apply targets monkey patches
#'
#' Patch selected targets internals used by the current workflow bootstrap.
#'
#' @param force Logical; TRUE reruns the setup step even when bootstrap state says it has already completed.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

apply_targets_monkey_patches <- function(force = FALSE) {
  if (isTRUE(bootstrap_state_env$targets_patches_applied) && !force) {
    return(invisible(FALSE))
  }

  tar_format_trace_chr <- deparse(targets:::tar_format_trace) |>
    stringr::str_c(collapse = "\n") |>
    stringr::str_replace(
      stringr::fixed("width <- min(getOption(\"width\"), 79L) - 4L"),
      "width <- 1000L"
    )

  tar_format_trace_modified <- parse(text = tar_format_trace_chr) |>
    eval()

  utils::assignInNamespace(
    x = "tar_format_trace",
    value = tar_format_trace_modified,
    ns = "targets"
  )

  hash_import_object_modified <- function(value, name, hashes, graph) {
    if (is.function(value) && (grepl("_noinval$", name) || isTRUE(attr(value, "skip_invalidate")))) {
      meta <- targets:::meta_init(path_store = targets::tar_config_get("store"))
      meta$database$ensure_storage()
      data <- meta$database$read_condensed_data()
      old_row <- data[data$name == name, , drop = FALSE]
      meta$database$close()

      if (nrow(old_row) == 1L && !is.na(old_row$data)) {
        assign(x = name, value = old_row$data, envir = hashes)
        return(invisible())
      }
    }

    UseMethod("hash_import_object", value)
  }

  utils::assignInNamespace(
    x = "hash_import_object",
    value = hash_import_object_modified,
    ns = "targets"
  )

  bootstrap_state_env$targets_patches_applied <- TRUE
  invisible(TRUE)
}

source_crew_controllers <- function(force = FALSE, envir = .GlobalEnv) {
  if (isTRUE(bootstrap_state_env$controllers_loaded) && !force && identical(envir, .GlobalEnv)) {
    return(invisible(FALSE))
  }

  controller_setup <- source(file.path(get_project_root(force = force), "crew_controllers.R"), local = envir, chdir = TRUE)$value
  apply_crew_controller_options(controller_setup)

  if (identical(envir, .GlobalEnv)) {
    bootstrap_state_env$controllers_loaded <- TRUE
  }
  invisible(TRUE)
}

load_project_runtime <- function(force = FALSE) {
  load_project_packages(force = force)
  source_shared_helpers(force = force)
  apply_runtime_options(force = force)
  apply_targets_monkey_patches(force = force)
  source_crew_controllers(force = force)
  invisible(TRUE)
}

load_interactive_helpers <- function(full = FALSE, force = FALSE) {
  load_project_packages(force = force)
  source_shared_helpers(force = force)

  if (isTRUE(full)) {
    load_project_runtime(force = force)
  }

  invisible(TRUE)
}

# TODO: Not cool to overwrite standard targets::tar_make() function - possible renaming?
tar_make_v2 <- function(names = NULL, ...) {
  load_project_runtime(force = TRUE)
  if (base::missing(names)) {
    tar_make_w_cfg_skip_patterns(...)
  } else {
    names_expr <- substitute(names)
    base::eval(rlang::expr(tar_make_w_cfg_skip_patterns(names = !!names_expr, !!!list(...))))
  }
  # targets::tar_make(...)
}

source_Rrofile <- function() {
  source(file.path(get_project_root(), ".Rprofile"))
}

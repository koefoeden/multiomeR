multiomeR_core_runtime_state_env <- new.env(parent = emptyenv())
multiomeR_core_runtime_state_env$controller_resources_tibble <- NULL

get_controller_names <- function(controller_list) {
  vapply(
    controller_list,
    \(controller) controller$launcher$name,
    character(1)
  )
}

#' Validate controller setup
#'
#' Validate the controller/resource table returned by the crew setup file.
#'
#' @param controller_setup Named list returned by `crew_controllers.R`, containing `controller_resources_tibble` and `controller_list`.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

validate_controller_setup <- function(controller_setup) {
  if (!is.list(controller_setup) || is.null(names(controller_setup))) {
    stop("crew_controllers.R must return a named list.", call. = FALSE)
  }

  required_setup_names <- c("controller_resources_tibble", "controller_list")
  missing_setup_names <- setdiff(required_setup_names, names(controller_setup))
  if (length(missing_setup_names) > 0) {
    stop(
      "crew_controllers.R return value is missing: ",
      paste(missing_setup_names, collapse = ", "),
      call. = FALSE
    )
  }

  controller_resources_tibble <- controller_setup$controller_resources_tibble
  controller_list <- controller_setup$controller_list
  required_resource_columns <- c("controller_name", "cores", "RAM_GB", "gpus")

  if (!is.data.frame(controller_resources_tibble)) {
    stop("controller_resources_tibble must be a data frame or tibble.", call. = FALSE)
  }

  if (!identical(names(controller_resources_tibble), required_resource_columns)) {
    stop(
      "controller_resources_tibble must have exactly these columns, in order: ",
      paste(required_resource_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(controller_resources_tibble) == 0) {
    stop("controller_resources_tibble must contain at least one controller row.", call. = FALSE)
  }

  if (anyDuplicated(controller_resources_tibble$controller_name) > 0) {
    stop("controller_resources_tibble$controller_name must be unique.", call. = FALSE)
  }

  numeric_resource_columns <- c("cores", "RAM_GB", "gpus")
  numeric_resource_cols <- controller_resources_tibble[numeric_resource_columns]
  if (!all(vapply(numeric_resource_cols, is.numeric, logical(1)))) {
    stop("controller_resources_tibble$cores, $RAM_GB, and $gpus must be numeric.", call. = FALSE)
  }

  if (any(is.na(unlist(numeric_resource_cols)))) {
    stop("controller_resources_tibble$cores, $RAM_GB, and $gpus cannot contain NA values.", call. = FALSE)
  }

  if (!is.list(controller_list) || length(controller_list) == 0) {
    stop("controller_list must be a non-empty list of crew controllers.", call. = FALSE)
  }

  controller_names <- get_controller_names(controller_list)
  if (anyDuplicated(controller_names) > 0) {
    stop("controller_list names must be unique.", call. = FALSE)
  }

  missing_controller_names <- setdiff(controller_resources_tibble$controller_name, controller_names)
  if (length(missing_controller_names) > 0) {
    stop(
      "controller_resources_tibble contains controller names absent from controller_list: ",
      paste(missing_controller_names, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Apply crew controller options
#'
#' Install crew controller and default resource options into `targets`.
#'
#' @param controller_setup Named list returned by `crew_controllers.R`, containing `controller_resources_tibble` and `controller_list`.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

apply_crew_controller_options <- function(controller_setup) {
  validate_controller_setup(controller_setup)

  controller_resources_tibble <- controller_setup$controller_resources_tibble
  controller_list <- controller_setup$controller_list
  default_controller_name <- controller_resources_tibble$controller_name[[1]]

  targets::tar_option_set(
    retrieval = "worker",
    controller = crew::crew_controller_group(controller_list),
    resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = default_controller_name)
    )
  )

  multiomeR_core_runtime_state_env$controller_resources_tibble <-
    controller_resources_tibble

  if (exists("bootstrap_state_env", inherits = TRUE)) {
    bootstrap_state_env$controller_resources_tibble <- controller_resources_tibble
  }

  invisible(controller_setup)
}

#' Get tar resources
#'
#' Select the smallest crew controller tier satisfying a target's resource needs.
#'
#' @param cores_req Minimum CPU cores required by the target.
#' @param RAM_GB_req Minimum base RAM in GB required by the target.
#' @param RAM_GB_per_extra_core Additional RAM in GB required per CPU core
#'   beyond the first requested core.
#' @param gpus_req Minimum GPU count required by the target; use 0 for CPU-only work.
#' @param controller_resources Tibble of available controller tiers with CPU, RAM, GPU, controller name, and default flags.
#' @return A `targets::tar_resources()` object pointing at the selected crew
#'   controller. Requests exceeding all configured tiers error with the request
#'   details.
#' @keywords internal

get_tar_resources <- function(
  cores_req = NULL,
  RAM_GB_req = NULL,
  RAM_GB_per_extra_core = NULL,
  gpus_req = 0,
  controller_resources = NULL
) {
  if (is.null(controller_resources)) {
    controller_resources <-
      multiomeR_core_runtime_state_env$controller_resources_tibble
    if (is.null(controller_resources)) {
      stop(
        "Controller resources are not loaded. Call ",
        "apply_crew_controller_options() first or provide ",
        "`controller_resources`.",
        call. = FALSE
      )
    }
  }

  default_controller_row <- controller_resources[1, , drop = FALSE]

  if (is.null(cores_req)) {
    cores_req <- default_controller_row$cores[[1]]
  }

  if (!is.numeric(cores_req) || length(cores_req) != 1 || is.na(cores_req)) {
    stop("cores_req must be a single numeric value.", call. = FALSE)
  }
  if (is.null(RAM_GB_req)) {
    RAM_GB_req <- default_controller_row$RAM_GB[[1]]
  }
  if (!is.numeric(RAM_GB_req) || length(RAM_GB_req) != 1 || is.na(RAM_GB_req)) {
    stop("RAM_GB_req must be a single numeric value in GB.", call. = FALSE)
  }
  if (!is.null(RAM_GB_per_extra_core)) {
    if (!is.numeric(RAM_GB_per_extra_core) || length(RAM_GB_per_extra_core) != 1 || is.na(RAM_GB_per_extra_core)) {
      stop("RAM_GB_per_extra_core must be a single numeric value in GB.", call. = FALSE)
    }
    if (RAM_GB_per_extra_core < 0) {
      stop("RAM_GB_per_extra_core cannot be negative.", call. = FALSE)
    }
    RAM_GB_req <- RAM_GB_req + max(cores_req - 1, 0) * RAM_GB_per_extra_core
  }
  if (!is.numeric(gpus_req) || length(gpus_req) != 1 || is.na(gpus_req)) {
    stop("gpus_req must be a single numeric value.", call. = FALSE)
  }
  if (gpus_req < 0) {
    stop("gpus_req cannot be negative.", call. = FALSE)
  }

  gpus_match <- if (gpus_req == 0) {
    controller_resources$gpus == 0
  } else {
    controller_resources$gpus >= gpus_req
  }

  eligible_controllers <- controller_resources[
    controller_resources$cores >= cores_req &
      controller_resources$RAM_GB >= RAM_GB_req &
      gpus_match,
    ,
    drop = FALSE
  ]

  if (nrow(eligible_controllers) == 0) {
    stop(
      paste0(
        "No controller found for cores_req = ", cores_req,
        " and RAM_GB_req = ", RAM_GB_req,
        " and gpus_req = ", gpus_req,
        ". Available controllers: ",
        paste(
          sprintf(
            "%s (%s cores, %s GB, %s GPUs)",
            controller_resources$controller_name,
            controller_resources$cores,
            controller_resources$RAM_GB,
            controller_resources$gpus
          ),
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }

  selected_controller_name <- eligible_controllers$controller_name[[1]]
  targets::tar_resources(crew = targets::tar_resources_crew(controller = selected_controller_name))
}

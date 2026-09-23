multiomeR_core_runtime_state_env <- new.env(parent = emptyenv())
multiomeR_core_runtime_state_env$controller_resources_tibble <- NULL

get_controller_names <- function(controller_list) {
  vapply(
    controller_list,
    \(controller) controller$launcher$name,
    character(1)
  )
}

order_controller_resources <- function(controller_resources) {
  controller_resources[
    order(
      controller_resources$gpus,
      controller_resources$cores,
      controller_resources$RAM_GB,
      controller_resources$controller_name
    ),
    ,
    drop = FALSE
  ]
}

#' Validate controller setup
#'
#' Validate the controller/resource table returned by the crew setup file.
#'
#' @param controller_setup Named list returned by `crew_controllers.R`, containing `controller_resources_tibble` and `controller_list`.
#' @return Invisibly returns after the validation succeeds.
#' @keywords internal

validate_controller_setup <- function(controller_setup) {
  resources <- controller_setup$controller_resources_tibble
  controller_names <- get_controller_names(controller_setup$controller_list)
  resource_values <- resources[c("cores", "RAM_GB", "gpus")]
  if (
    !identical(names(resources), c("controller_name", "cores", "RAM_GB", "gpus")) ||
      nrow(resources) == 0L ||
      anyDuplicated(resources$controller_name) ||
      anyDuplicated(controller_names) ||
      !all(vapply(resource_values, is.numeric, logical(1))) ||
      anyNA(resource_values) ||
      !all(resources$controller_name %in% controller_names)
  ) {
    stop(
      "crew_controllers.R must return `controller_resources_tibble` with columns ",
      "controller_name, cores, RAM_GB, and gpus (unique names, complete numeric ",
      "resources) and a `controller_list` with unique names covering every ",
      "controller_name.",
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

  controller_resources_tibble <- order_controller_resources(
    controller_setup$controller_resources_tibble
  )
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

  invisible(controller_setup)
}

#' Get tar resources
#'
#' Select the smallest crew controller tier satisfying a target's resource needs.
#'
#' @param cores_req Minimum CPU cores required by the target.
#' @param RAM_GB_req Minimum RAM in GB required by the target.
#' @param gpus_req Minimum GPU count required by the target; use 0 for CPU-only work.
#' @return A `targets::tar_resources()` object pointing at the selected crew
#'   controller. Requests exceeding all configured tiers error with the request
#'   details.
#' @keywords internal

get_tar_resources <- function(cores_req = 1, RAM_GB_req = 0, gpus_req = 0) {
  controller_resources <- multiomeR_core_runtime_state_env$controller_resources_tibble
  if (is.null(controller_resources)) {
    stop("Controller resources are not loaded. Call apply_crew_controller_options() first.", call. = FALSE)
  }
  gpus_match <- if (gpus_req == 0) controller_resources$gpus == 0 else controller_resources$gpus >= gpus_req
  eligible <- which(
    controller_resources$cores >= cores_req & controller_resources$RAM_GB >= RAM_GB_req & gpus_match
  )
  if (!length(eligible)) {
    stop(
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
      ),
      call. = FALSE
    )
  }
  targets::tar_resources(
    crew = targets::tar_resources_crew(controller = controller_resources$controller_name[[eligible[[1]]]])
  )
}

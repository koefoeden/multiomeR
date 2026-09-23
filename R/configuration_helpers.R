#' Locate a file in this checkout's selected configuration directory
#'
#' The optional, untracked configuration.local contains one directory path.
#' Relative directory paths are anchored at the repository root. Paths inside
#' configuration files retain their existing interpretation.
#' @param filename Configuration basename, including its extension.
#' @param project_root Repository root; explicit roots are useful for validation.
#' @param must_exist Whether the requested file must exist. Disabled modules
#'   may omit their configuration file.
#' @return Absolute path to the selected configuration file.
configuration_path <- function(filename, project_root = get_project_root(), must_exist = TRUE) {
  selection_file <- file.path(project_root, "configuration.local")
  directory <- "configuration"
  if (file.exists(selection_file)) {
    directory <- trimws(readLines(selection_file, warn = FALSE))
    if (length(directory) != 1L || !nzchar(directory)) {
      stop(selection_file, " must contain exactly one non-empty directory path.", call. = FALSE)
    }
  }
  directory <- path.expand(directory)
  if (!startsWith(directory, "/")) directory <- file.path(project_root, directory)
  if (!dir.exists(directory)) {
    stop("Selected configuration directory does not exist: ", directory, call. = FALSE)
  }
  path <- file.path(normalizePath(directory, winslash = "/", mustWork = TRUE), filename)
  if (must_exist && !file.exists(path)) {
    stop("Missing file in selected configuration directory: ", path, call. = FALSE)
  }
  path
}

read_config_parameter_manifest <- function(manifest_file, scope = NULL) {
  manifest_tibble <- readr::read_tsv(
    manifest_file,
    col_types = readr::cols(
      .default = readr::col_character(),
      allow_missing_after_inheritance = readr::col_logical()
    ),
    na = character(),
    show_col_types = FALSE
  )
  duplicated_params <- with(manifest_tibble, paste(scope, param_name, sep = "/"))
  duplicated_params <- unique(duplicated_params[duplicated(duplicated_params)])
  if (length(duplicated_params) > 0) {
    stop(manifest_file, " contains duplicated scope/param_name value(s): ", paste(duplicated_params, collapse = ", "), call. = FALSE)
  }

  if (!is.null(scope)) {
    manifest_tibble <- dplyr::filter(manifest_tibble, .data$scope == .env$scope)
    if (nrow(manifest_tibble) == 0) {
      stop(manifest_file, " does not define any parameters for scope: ", scope, call. = FALSE)
    }
  }

  dplyr::mutate(
    manifest_tibble,
    default_value = purrr::map(default_value, \(x) yaml::read_yaml(text = paste0("value: ", x), eval.expr = TRUE)$value),
    allowed_values = purrr::map(allowed_values, \(x) if (x == "") character() else stringr::str_split_1(x, ","))
  )
}

read_aggregation_config_tibble <- function(
  config_file = configuration_path("cfg_aggregations.yaml"),
  manifest_file = "cfg_pipeline_parameters.tsv"
) {
  config <- read_manifest_config_tibble(config_file, manifest_file, scope = "aggregation", key_col = "aggregation")
  selected <- validation_aggregations()
  if (length(selected)) {
    missing <- setdiff(selected, config$aggregation)
    if (length(missing)) stop("Missing validation aggregations: ", paste(missing, collapse = ", "))
    config$is_active <- as.list(config$aggregation %in% selected)
  }
  config
}

#' Read a configuration YAML against the parameter manifest
#'
#' Each entry starts from the manifest defaults, applies its `inherits` parents
#' in order, then its own values, and is validated against the manifest.
#' @return One row per configuration entry, keyed by `key_col`.
read_manifest_config_tibble <- function(config_file, manifest_file, scope, key_col) {
  manifest_tibble <- read_config_parameter_manifest(manifest_file, scope = scope)
  raw_cfg <- suppressWarnings(yaml::read_yaml(config_file, eval.expr = TRUE))
  defaults <- purrr::set_names(manifest_tibble$default_value, manifest_tibble$param_name)

  for (config_key in names(raw_cfg)) {
    unknown_params <- setdiff(names(raw_cfg[[config_key]]), c("inherits", manifest_tibble$param_name))
    if (length(unknown_params) > 0) {
      stop(
        "Config ", key_col, " '", config_key, "' in ", config_file,
        " defines parameter(s) absent from the manifest: ", paste(unknown_params, collapse = ", "),
        call. = FALSE
      )
    }
  }

  resolve_values <- function(config_key, seen_configs = character()) {
    if (!config_key %in% names(raw_cfg)) {
      stop("Config ", key_col, " '", config_key, "' is not defined in ", config_file, ".", call. = FALSE)
    }
    if (config_key %in% seen_configs) {
      stop(
        "Circular config inheritance found in ", config_file, ": ",
        paste(c(seen_configs, config_key), collapse = " -> "), ".",
        call. = FALSE
      )
    }
    cfg_list <- raw_cfg[[config_key]]
    values <- defaults
    for (parent_config in cfg_list$inherits) {
      parent_values <- resolve_values(parent_config, c(seen_configs, config_key))
      values[names(parent_values)] <- parent_values
    }
    child_values <- cfg_list[setdiff(names(cfg_list), "inherits")]
    values[names(child_values)] <- child_values
    values
  }

  validate_values <- function(values, config_key) {
    purrr::pwalk(
      manifest_tibble,
      \(param_name, data_type, cardinality, allow_missing_after_inheritance, allowed_values, ...) {
        value <- values[[param_name]]
        problem <- if (is_manifest_missing_value(value)) {
          if (!isTRUE(allow_missing_after_inheritance)) {
            "must not resolve to NULL/NA/empty; set it or provide a non-missing default or inherited value"
          }
        } else if (!isTRUE(switch(cardinality,
          scalar = length(value) == 1 && !is.list(value),
          vector = is.atomic(value),
          named_list = is.list(value) && !is.null(names(value)) && all(nzchar(names(value))),
          list = is.list(value)
        ))) {
          paste("has invalid cardinality; expected", cardinality)
        } else if (!isTRUE(switch(data_type,
          character = ,
          path = ,
          regex = is.character(value),
          logical = is.logical(value),
          numeric = is.numeric(value),
          integer = is.numeric(value) && all(value == as.integer(value)),
          named_list = ,
          list = is.list(value)
        ))) {
          paste("has invalid type; expected", data_type)
        } else if (length(allowed_values) > 0 && length(invalid_values <- setdiff(as.character(value), allowed_values)) > 0) {
          paste0(
            "has invalid value(s) ", paste(invalid_values, collapse = ", "),
            "; allowed values: ", paste(allowed_values, collapse = ", ")
          )
        }
        if (!is.null(problem)) {
          stop("Parameter '", param_name, "' ", problem, " for config ", key_col, ": ", config_key, " in ", config_file, ".", call. = FALSE)
        }
      }
    )
  }

  purrr::map(purrr::set_names(names(raw_cfg)), \(config_key) {
    values <- resolve_values(config_key)
    validate_values(values, config_key)
    tidyr::pivot_wider(tibble::enframe(values))
  }) |>
    dplyr::bind_rows(.id = key_col)
}

is_manifest_missing_value <- function(value) {
  if (is.null(value)) {
    return(TRUE)
  }
  if (!is.atomic(value) || length(value) == 0) {
    return(FALSE)
  }
  all(is.na(value) | trimws(as.character(value)) == "")
}

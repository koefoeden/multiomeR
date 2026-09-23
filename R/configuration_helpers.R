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

  duplicate_param_tibble <- manifest_tibble |>
    dplyr::count(scope, param_name) |>
    dplyr::filter(n > 1)
  if (nrow(duplicate_param_tibble) > 0) {
    stop(
      manifest_file,
      " contains duplicated scope/param_name value(s): ",
      paste(
        stringr::str_c(duplicate_param_tibble$scope, duplicate_param_tibble$param_name, sep = "/"),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (!is.null(scope)) {
    manifest_tibble <- manifest_tibble |>
      dplyr::filter(.data$scope == .env$scope)

    if (nrow(manifest_tibble) == 0) {
      stop(manifest_file, " does not define any parameters for scope: ", scope, call. = FALSE)
    }
  }

  manifest_tibble |>
    dplyr::mutate(
      default_value = purrr::map(default_value, parse_manifest_default_value),
      allowed_values = purrr::map(allowed_values, parse_manifest_allowed_values)
    )
}

parse_manifest_default_value <- function(default_value_chr) {
  yaml::read_yaml(text = paste0("value: ", default_value_chr), eval.expr = TRUE)$value
}

parse_manifest_allowed_values <- function(allowed_values_chr) {
  if (is.na(allowed_values_chr) || allowed_values_chr == "") {
    return(character())
  }
  stringr::str_split_1(allowed_values_chr, ",")
}

manifest_defaults <- function(manifest_tibble) {
  defaults <- manifest_tibble$default_value
  names(defaults) <- manifest_tibble$param_name
  defaults
}

read_aggregation_config_tibble <- function(
  config_file = configuration_path("cfg_aggregations.yaml"),
  manifest_file = "cfg_pipeline_parameters.tsv",
  verbose = FALSE
) {
  config <- read_manifest_config_tibble(
    config_file = config_file,
    manifest_file = manifest_file,
    scope = "aggregation",
    key_col = "aggregation",
    verbose = verbose
  )
  selected <- validation_aggregations()
  if (length(selected)) {
    missing <- setdiff(selected, config$aggregation)
    if (length(missing)) stop("Missing validation aggregations: ", paste(missing, collapse = ", "))
    config$is_active <- as.list(config$aggregation %in% selected)
  }
  config
}

read_manifest_config_tibble <- function(config_file, manifest_file, scope, key_col, verbose = FALSE) {
  manifest_tibble <- read_config_parameter_manifest(manifest_file, scope = scope)
  raw_cfg <- suppressWarnings(yaml::read_yaml(config_file, eval.expr = TRUE))
  config_keys <- names(raw_cfg)

  if (verbose) {
    cat("Reading configs for ", key_col, ":\n", paste("-", config_keys, collapse = "\n"), "\n\n", sep = "")
  }

  validate_manifest_config_names(raw_cfg, config_keys, manifest_tibble, config_file, key_col = key_col)

  values_per_key <- config_keys |>
    purrr::set_names() |>
    purrr::map(
      ~ resolve_manifest_config_values(
        raw_cfg = raw_cfg,
        config_key = .x,
        config_file = config_file,
        manifest_tibble = manifest_tibble,
        key_col = key_col
      )
  )

  per_key_tibbles <- purrr::imap(
    values_per_key,
    \(values, config_key) {
      validate_manifest_config_values(
        values = values,
        config_key = config_key,
        config_file = config_file,
        manifest_tibble = manifest_tibble,
        key_col = key_col
      )
      values |>
        tibble::enframe() |>
        tidyr::pivot_wider()
    }
  )

  dplyr::bind_rows(per_key_tibbles, .id = key_col)
}

validate_manifest_config_names <- function(raw_cfg, config_keys, manifest_tibble, config_file, key_col) {
  purrr::walk(config_keys, \(config_key) {
    entry_params <- setdiff(names(raw_cfg[[config_key]]), "inherits")
    unknown_params <- setdiff(entry_params, manifest_tibble$param_name)
    if (length(unknown_params) > 0) {
      stop(
        "Config ",
        key_col,
        " '",
        config_key,
        "' in ",
        config_file,
        " defines parameter(s) absent from the manifest: ",
        paste(unknown_params, collapse = ", "),
        call. = FALSE
      )
    }
  })
}

resolve_manifest_config_values <- function(
  raw_cfg,
  config_key,
  config_file,
  manifest_tibble,
  key_col,
  seen_configs = character()
) {
  if (!config_key %in% names(raw_cfg)) {
    stop(stringr::str_glue("Config {key_col} '{config_key}' is not defined in {config_file}."), call. = FALSE)
  }

  if (config_key %in% seen_configs) {
    stop(stringr::str_glue(
      "Circular config inheritance found in {config_file}: ",
      "{paste(c(seen_configs, config_key), collapse = ' -> ')}."
    ), call. = FALSE)
  }

  cfg_list <- raw_cfg[[config_key]]
  parent_configs <- cfg_list$inherits %||% character()
  values <- manifest_defaults(manifest_tibble)

  for (parent_config in parent_configs) {
    values <- merge_manifest_config_values(
      parent_values = values,
      child_values = resolve_manifest_config_values(
        raw_cfg = raw_cfg,
        config_key = parent_config,
        config_file = config_file,
        manifest_tibble = manifest_tibble,
        key_col = key_col,
        seen_configs = c(seen_configs, config_key)
      )
    )
  }

  merge_manifest_config_values(
    parent_values = values,
    child_values = cfg_list[setdiff(names(cfg_list), "inherits")]
  )
}

merge_manifest_config_values <- function(parent_values, child_values) {
  parent_values[names(child_values)] <- child_values
  parent_values
}

validate_manifest_config_values <- function(values, config_key, config_file, manifest_tibble, key_col) {
  purrr::pwalk(
    manifest_tibble,
    \(scope, param_name, data_type, cardinality, default_value, allow_missing_after_inheritance, allowed_values, ...) {
      value <- values[[param_name]]
      if (is_manifest_missing_value(value)) {
        if (!isTRUE(allow_missing_after_inheritance)) {
          stop(
            "Parameter '",
            param_name,
            "' is not allowed to resolve to NULL/NA/empty for config ",
            key_col,
            ": ",
            config_key,
            ". Please set this parameter in ",
            config_file,
            " or provide a non-missing default/inherited value.",
            call. = FALSE
          )
        }
        return(NULL)
      }

      validate_manifest_cardinality(value, cardinality, param_name, config_key, config_file, key_col)
      validate_manifest_type(value, data_type, param_name, config_key, config_file, key_col)
      validate_manifest_allowed_values(value, allowed_values, param_name, config_key, config_file, key_col)
    }
  )

  invisible(values)
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

validate_manifest_cardinality <- function(value, cardinality, param_name, config_key, config_file, key_col) {
  valid <- switch(cardinality,
    scalar = length(value) == 1 && !is.list(value),
    vector = is.atomic(value),
    named_list = is.list(value) && !is.null(names(value)) && all(nzchar(names(value))),
    list = is.list(value),
    stop("Unknown manifest cardinality '", cardinality, "' for parameter '", param_name, "'.", call. = FALSE)
  )

  if (!isTRUE(valid)) {
    stop(
      "Parameter '",
      param_name,
      "' has invalid cardinality for config ",
      key_col,
      ": ",
      config_key,
      " in ",
      config_file,
      ". Expected ",
      cardinality,
      ".",
      call. = FALSE
    )
  }
}

validate_manifest_type <- function(value, data_type, param_name, config_key, config_file, key_col) {
  valid <- switch(data_type,
    character = is.character(value),
    path = is.character(value),
    regex = is.character(value),
    logical = is.logical(value),
    numeric = is.numeric(value),
    integer = is.numeric(value) && all(value == as.integer(value)),
    named_list = is.list(value),
    list = is.list(value),
    stop("Unknown manifest data_type '", data_type, "' for parameter '", param_name, "'.", call. = FALSE)
  )

  if (!isTRUE(valid)) {
    stop(
      "Parameter '",
      param_name,
      "' has invalid type for config ",
      key_col,
      ": ",
      config_key,
      " in ",
      config_file,
      ". Expected ",
      data_type,
      ".",
      call. = FALSE
    )
  }
}

validate_manifest_allowed_values <- function(value, allowed_values, param_name, config_key, config_file, key_col) {
  if (length(allowed_values) == 0) {
    return(invisible(NULL))
  }

  invalid_values <- setdiff(as.character(value), allowed_values)
  if (length(invalid_values) > 0) {
    stop(
      "Parameter '",
      param_name,
      "' has invalid value(s) for config ",
      key_col,
      ": ",
      config_key,
      " in ",
      config_file,
      ": ",
      paste(invalid_values, collapse = ", "),
      ". Allowed values: ",
      paste(allowed_values, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(NULL)
}

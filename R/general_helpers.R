`%||%` <- function(x, y) if (is.null(x)) y else x

skip_invalidate <- function(fn) {
  attr(fn, "skip_invalidate") <- TRUE
  fn
}
# assign("skip_invalidate", skip_invalidate, envir = .GlobalEnv) # check if this is really needed.

#' Collapse duplicate names
#'
#' Make duplicated list names unique by appending counters while preserving original ordering.
#'
#' @param input_list_object A list, possibly with duplicated names at one or more
#'   nested levels. Non-list objects and unnamed lists are returned unchanged.
#' @return A list with elements sharing the same name collapsed under that name;
#'   nested list values are processed recursively.
#' @keywords internal

collapse_duplicate_names <- function(input_list_object) {
  # Return immediately if input is not a list or has no names to collapse
  if (!is.list(input_list_object) || is.null(names(input_list_object))) {
    return(input_list_object)
  }

  input_list_object %>%
    # Split the list by names; this groups items sharing the same key
    split(names(.)) %>%
    purrr::map(function(grouped_elements_list) {
      # unlist with recursive = FALSE concatenates the elements at the top level
      # If elements are vectors: list(c(1), c(2)) -> c(1, 2)
      # If elements are lists: list(list(a=1), list(a=2)) -> list(a=1, a=2)
      merged_content_object <- unlist(unname(grouped_elements_list), recursive = FALSE)

      # Recursively apply logic if the result is a list (to handle nested duplicates)
      if (is.list(merged_content_object)) {
        collapse_duplicate_names(merged_content_object)
      } else {
        merged_content_object
      }
    })
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

read_dataset_config_tibble <- function(
  config_file = "cfg_datasets.yaml",
  manifest_file = "cfg_pipeline_parameters.tsv",
  verbose = FALSE
) {
  read_manifest_config_tibble(
    config_file = config_file,
    manifest_file = manifest_file,
    scope = "dataset",
    key_col = "dataset",
    verbose = verbose
  )
}

read_aggregation_config_tibble <- function(
  config_file = "cfg_aggregations.yaml",
  manifest_file = "cfg_pipeline_parameters.tsv",
  verbose = FALSE
) {
  read_manifest_config_tibble(
    config_file = config_file,
    manifest_file = manifest_file,
    scope = "aggregation",
    key_col = "aggregation",
    verbose = verbose
  )
}

read_manifest_config_tibble <- function(config_file, manifest_file, scope, key_col, verbose = FALSE) {
  if (requireNamespace("googlesheets4", quietly = TRUE)) {
    googlesheets4::gs4_deauth()
  }

  manifest_tibble <- read_config_parameter_manifest(manifest_file, scope = scope)
  raw_cfg <- suppressWarnings(yaml::read_yaml(config_file, eval.expr = TRUE))
  config_keys <- names(raw_cfg) |>
    purrr::discard(~ .x %in% c("default"))

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
    entry_params <- names(flatten_manifest_config_entry(
      cfg_list = raw_cfg[[config_key]],
      config_key = config_key,
      config_file = config_file,
      key_col = key_col
    ))
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
  if (!config_key %in% names(raw_cfg) || config_key == "default") {
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
    child_values = flatten_manifest_config_entry(
      cfg_list = cfg_list,
      config_key = config_key,
      config_file = config_file,
      key_col = key_col
    )
  )
}

flatten_manifest_config_entry <- function(cfg_list, config_key, config_file, key_col) {
  cfg_names <- names(cfg_list)
  cfg_list[setdiff(cfg_names, "inherits")]
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


get_cfg_per_dataset_tibble <- function(config_file = "cfg_datasets.yaml", verbose = TRUE) {
  cfg_tibble <- read_dataset_config_tibble(config_file = config_file, verbose = verbose)

  return(cfg_tibble)
}


load_CFG <- function(name) {
  get_cfg_per_dataset_tibble(verbose = FALSE) %>%
    dplyr::filter(dataset == name) %>%
    as.list() %>%
    purrr::flatten() %>%
    list2env(envir = .GlobalEnv)
}

cfg_is_set <- function(x) {
  !is.null(x) ## , msg = null_cfg_msg)
}

null_cfg_msg <- " - required cfg-variable not specified in config. Remaining targets are skipped."

assertthat::on_failure(cfg_is_set) <- function(call, env) {
  paste0(deparse(call$x), null_cfg_msg)
}

assert_cfg_is_set <- function(x, cfg_var_name = NULL) {
  # This
  cfg_var_name <- rlang::enexpr(x) %||% cfg_var_name
  if (is.null(cfg_var_name)) {
    info_string <- "Unexpected NULL value found in target {tar_name_wo_suffixes()} for some configured parameter. Please rectify this."
  } else if (is.null(x)) {
    info_string <- "Required parameter '{cfg_var_name}' for target {tar_name_wo_suffixes()} is not specified. Please rectify this."
  } else {
    return(x)
  }
  message <- stringr::str_glue(
    info_string,
    "If you do not wish to run this target, specify this in tar_make(names = ...), or add it to cfg_tar_make_skip_regex_patterns and use tar_make_w_cfg_skip_patterns().",
    .sep = "\n"
  )
  stop(message)
}


ask_to_continue <- function(type = c("stop", "bool", "break")[1]) {
  answer <- readline(prompt = "Do you want to continue? (y/n): ")

  if (tolower(answer) %in% c("y", "yes")) {
    return(invisible(TRUE))
  } else {
    switch(
      type,
      "stop" = stop("Exiting R session..."),
      "bool" = return(invisible(FALSE)),
      "break" = break
    )
  }
}


#' Suppress warnings matching
#'
#' Evaluate an expression while suppressing only warnings whose messages match the requested pattern.
#'
#' @param expr Expression evaluated in the caller environment, for example while holding a file lock or suppressing matching warnings.
#' @param pattern Character string or vector of warning-message patterns to muffle.
#' @param fixed Logical passed to `grepl()`; use `TRUE` when `pattern` should be
#'   matched literally rather than as a regular expression.
#' @param ignore_case Logical passed to `grepl()` for case-insensitive matching.
#' @return The value of `expr`. Warnings whose messages do not match `pattern`
#'   are allowed to propagate normally.
#' @keywords internal

suppress_warnings_matching <- function(expr, pattern, fixed = FALSE, ignore_case = FALSE) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      is_match <- if (length(pattern) == 1L) {
        grepl(pattern, msg, fixed = fixed, ignore.case = ignore_case)
      } else {
        any(vapply(pattern, function(p) grepl(p, msg, fixed = fixed, ignore.case = ignore_case), logical(1)))
      }
      if (is_match) invokeRestart("muffleWarning")
      # else: fall through and show the warning
    }
  )
}

'%!in%' <- function(x, y) !('%in%'(x, y))


#' Run shell with glue
#'
#' Interpolate and execute an external command with optional clean-environment handling.
#'
#' @param command_string Shell command template. `{{...}}` expressions are
#'   interpolated with `glue::glue()` in the caller environment before running.
#' @param clean_env Logical; when `TRUE`, run through `clean_env_run.sh` to avoid
#'   leaking the current R session environment into external tools.
#' @return The processx result object, or the final command argument when called through `run_w_error_check()` wrappers.
#' @keywords internal

run_shell_with_glue <- function(command_string, clean_env = FALSE) {
  command_string_formatted <- glue::glue(
    stringr::str_replace_all(command_string, "\n", " "),
    .open = "{{",
    .close = "}}",
    .envir = parent.frame()
  )

  if (clean_env) {
    # Run in clean environment like run_w_error_check did
    full_command <- stringr::str_glue("env -i HOME=$HOME bash -lc '{command_string_formatted}'")
  } else {
    full_command <- as.character(command_string_formatted)
  }

  result_list <- processx::run(
    command = "bash",
    args = c("-c", full_command),
    error_on_status = FALSE,
    stdout = "|",
    stderr = "|"
  )

  # Check status manually
  if (result_list$status != 0) {
    base::stop(
      "Command failed (exit code ",
      result_list$status,
      "):\n",
      result_list$stderr,
      "\nstdout:\n",
      result_list$stdout
    )
  }

  return(result_list)
}

get_embedding_matrix_from_metadata <- function(metadata_tibble, umap_cols) {
  embedding_matrix <- metadata_tibble |>
    dplyr::select(dplyr::all_of(umap_cols)) |>
    as.matrix()
  rownames(embedding_matrix) <- NULL
  colnames(embedding_matrix) <- umap_cols
  embedding_matrix
}

#' Plot UMAP from metadata
#'
#' Draw a UMAP overlay for one metadata column or feature-expression row.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param variable Metadata column or feature row to plot.
#' @param value_source Whether `variable` should be read from `metadata_tibble`
#'   or from `feature_matrix`.
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param umap_cols Two metadata columns used as embedding x/y coordinates.
#' @param labels_discrete Logical passed to `BPCells::plot_embedding()`; when
#'   `TRUE`, discrete labels are drawn on the embedding.
#' @param legend_continuous Continuous legend mode passed to `BPCells::plot_embedding()`,
#'   for example `quantile`.
#' @param rasterize Logical; when `TRUE`, rasterize point layers for large cell
#'   embeddings.
#' @param raster_pixels Pixel width/height used for rasterized embedding layers.
#' @param randomize_order Logical; when `TRUE`, randomize plotting order to avoid
#'   systematic overplotting by cell order.
#' @param ... Additional arguments passed to `BPCells::plot_embedding()`.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_UMAP_from_metadata <- function(
  metadata_tibble,
  variable,
  value_source = c("metadata", "feature"),
  feature_matrix = NULL,
  umap_cols = c("LSI_UMAP_1", "LSI_UMAP_2"),
  labels_discrete = FALSE,
  legend_continuous = "quantile",
  rasterize = TRUE,
  raster_pixels = 1024,
  randomize_order = TRUE,
  ...
) {
  value_source <- match.arg(value_source)
  if (!is.character(variable) || length(variable) != 1 || is.na(variable) || !nzchar(variable)) {
    stop("`variable` must be a non-empty length-1 character vector.", call. = FALSE)
  }

  missing_umap_cols <- setdiff(umap_cols, colnames(metadata_tibble))
  if (length(missing_umap_cols) > 0) {
    stop("Required UMAP column(s) not found: ", paste(missing_umap_cols, collapse = ", "))
  }

  plot_col_data <- switch(
    value_source,
    metadata = {
      if (!variable %in% colnames(metadata_tibble)) {
        stop("Required metadata column not found: ", variable, call. = FALSE)
      }
      metadata_tibble |>
        dplyr::select(dplyr::all_of(variable)) |>
        as.data.frame()
    },
    feature = {
      if (is.null(feature_matrix)) {
        stop("`feature_matrix` is required when `value_source = 'feature'`.", call. = FALSE)
      }
      if (!variable %in% rownames(feature_matrix)) {
        stop("Required feature row not found: ", variable, call. = FALSE)
      }
      metadata_tibble |>
        dplyr::select(barcode_w_prefix) |>
        add_feature_matrix_to_metadata(
          feature_matrix = feature_matrix,
          features = variable
        ) |>
        dplyr::select(dplyr::all_of(variable)) |>
        as.data.frame()
    }
  )
  rownames(plot_col_data) <- NULL

  alpha_size_list <- get_BPCells_plot_embedding_aesthetics(metadata_tibble, rasterize = rasterize)
  embedding_matrix <- get_embedding_matrix_from_metadata(metadata_tibble, umap_cols)
  plot_embedding_matrix <- embedding_matrix

  if (is.numeric(plot_col_data[[1]])) {
    keep_rows <- is.finite(plot_col_data[[1]])
    if (!any(keep_rows)) {
      return(structure(list(), class = c("empty_plot_list", "list")))
    }
    plot_col_data <- plot_col_data[keep_rows, , drop = FALSE]
    plot_embedding_matrix <- plot_embedding_matrix[keep_rows, , drop = FALSE]
  }

  BPCells::plot_embedding(
    source = plot_col_data,
    embedding = plot_embedding_matrix,
    features = variable,
    size = alpha_size_list$size,
    rasterize = rasterize,
    raster_pixels = raster_pixels,
    randomize_order = randomize_order,
    labels_discrete = labels_discrete,
    legend_continuous = legend_continuous,
    return_plot_list = TRUE,
    apply_styling = TRUE,
    ...
  ) |>
    apply_alpha_to_plot_embedding(alpha_size_list$alpha) |>
    (\(plot) plot + ggplot2::labs(title = variable))() # feature var name is dropped for one-column dataframes, and replaced with generic "value" so we add it back as title. PR candidate for BPCells.
}

#' Plot 3 by 3 clusters and reduction UMAPs from metadata
#'
#' Compare RNA, WNN, and ATAC cluster labels across their UMAP reductions.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param cluster_col_suffix Which cluster annotation suffix to plot across RNA,
#'   WNN, and ATAC reductions. Must be `named` or `cell_type`.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_3_by_3_clusters_and_reduction_UMAPs_from_metadata <- function(metadata_tibble, cluster_col_suffix = c("named", "cell_type")) {
  cluster_col_suffix <- match.arg(cluster_col_suffix)

  dim_plot_list <- purrr::pmap(
    tibble::tribble(
      ~cluster_name    , ~cluster_col                                          ,
      "RNA clusters"   , stringr::str_c("PCA_harmony_SNN_cluster_", cluster_col_suffix) ,
      "Joint clusters" , stringr::str_c("WNN_harmony_SNN_cluster_", cluster_col_suffix) ,
      "ATAC clusters"  , stringr::str_c("LSI_harmony_SNN_cluster_", cluster_col_suffix)
    ),
    \(cluster_name, cluster_col) {
      tri_plot_list <- purrr::pmap(
        tibble::tribble(
          ~umap_cols                    , ~modality_title ,
          c("GEX_UMAP_1", "GEX_UMAP_2") , "RNA"           ,
          c("WNN_UMAP_1", "WNN_UMAP_2") , "Joint"         ,
          c("LSI_UMAP_1", "LSI_UMAP_2") , "ATAC"
        ),
        \(umap_cols, modality_title) {
          metadata_tibble |>
            plot_UMAP_from_metadata(
              variable = cluster_col,
              umap_cols = umap_cols
            ) +
            CONST_UMAP_ggplot2_theme +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(title = modality_title, x = NULL, y = cluster_name, color = cluster_col)
        }
      )

      patchwork::wrap_plots(tri_plot_list, guides = "collect")
    }
  )

  patchwork::wrap_plots(dim_plot_list, ncol = 1)
}


# Low-level POST to Open Targets GraphQL

get_data_from_exec_query <- function(query_file, variables_list, graph_ql_client = NULL) {
  query_str <- query_file %>% readr::read_lines() %>% paste0(collapse = "\n")

  if (is.null(graph_ql_client)) {
    graph_ql_client <- ghql::GraphqlClient$new(url = "https://api.platform.opentargets.org/api/v4/graphql")
  }
  graph_ql_query <- ghql::Query$new()
  graph_ql_query$query(name = "query", x = query_str)

  query_result <- graph_ql_client$exec(graph_ql_query$queries$query, variables_list, flatten = TRUE) %>% jsonlite::fromJSON()

  return(query_result$data)
}


assert_with_info <- function(..., glue_info = NULL, env = parent.frame()) {
  # Check validity without throwing error immediately
  validation_result_obj <- assertthat::validate_that(..., env = env)

  # If TRUE, proceed
  if (isTRUE(validation_result_obj)) {
    return(TRUE)
  }

  # If not TRUE, validation_result_obj is the default error string.
  # Combine default message with custom context.
  final_message_chr <- validation_result_obj

  if (!is.null(glue_info)) {
    final_message_chr <- paste0(final_message_chr, ". \nInfo: ", stringr::str_glue(glue_info, .envir = env))
  }

  stop(final_message_chr, call. = FALSE)
}


assert_command_on_path <- function(command) {
  if (file.exists(command) || nzchar(Sys.which(command))) {
    return(invisible(TRUE))
  }

  stop(
    "Command '", command, "' was not found on PATH. ",
    "Run through `pixi run ...` or add the tool to pixi.toml.",
    call. = FALSE
  )
}

run_w_error_check <- skip_invalidate(
  function(
    command_string,
    arguments_chr = character(),
    in_shell = FALSE,
    conda_env_path = NULL,
    conda_bin = NULL,
    in_clean_env = FALSE,
    in_minimal_env = FALSE,
    modules = character(),
    inherited_unnamed_env_vars = c("HOME", "TMPDIR", "MODULEPATH"),
    new_named_env_vars = character(),
    ...
  ) {
    processx_inherited_env_vars <- Sys.getenv(inherited_unnamed_env_vars)
    names(processx_inherited_env_vars) <- inherited_unnamed_env_vars
    modules <- modules %||% character()

    if (!is.null(conda_env_path)) {
      arguments_chr <- c("run", "-p", conda_env_path, "--no-capture-output", command_string, arguments_chr)
      command_string <- if (is.null(conda_bin)) "conda" else conda_bin
    }

    use_shell_wrapper <- in_shell || in_minimal_env || length(modules) > 0
    if (!use_shell_wrapper) {
      assert_command_on_path(command_string)
    }

    if (use_shell_wrapper) {
      wrapped_command_string <- if (in_shell) {
        stringr::str_flatten(c(command_string, arguments_chr), " ")
      } else {
        c(command_string, arguments_chr) %>% shQuote() %>% stringr::str_flatten(" ")
      }
      shell_steps <- c(
        "source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/Modules/init/bash 2>/dev/null",
        if (in_minimal_env && length(modules) > 0) "module purge",
        if (length(modules) > 0) paste("module load --auto", modules %>% shQuote() %>% stringr::str_flatten(" ")),
        paste("exec", wrapped_command_string)
      ) %>%
        purrr::discard(is.null)
    }

    if (use_shell_wrapper && in_minimal_env) {
      bash_inherited_env_vars <- processx_inherited_env_vars %>% purrr::imap_vec(~ stringr::str_c(.y, "=", .x))
      bash_new_env_vars <- stringr::str_c(names(new_named_env_vars), "=", new_named_env_vars)
      arguments_chr <- c(
        "-i",
        bash_inherited_env_vars,
        bash_new_env_vars,
        "bash",
        "--noprofile",
        "--norc",
        "-lc",
        stringr::str_flatten(shell_steps, " && ")
      )
      command_string <- stringr::str_c("env")
    } else if (use_shell_wrapper) {
      arguments_chr <- c("-lc", stringr::str_flatten(shell_steps, " && "))
      command_string <- Sys.getenv("SHELL") %||% "bash"
    }
    # else if (in_clean_env) {
    #   arguments_chr <- c("-i", evaluated_env_vars, command_string, arguments_chr)
    #   command_string <- str_c("env")
    # }

    result_list <- processx::run(
      command = command_string,
      args = arguments_chr,
      error_on_status = FALSE,
      stdout_line_callback = function(line, proc) {
        cat(line, "\n", file = stdout())
        flush(stdout())
      },
      stderr_line_callback = function(line, proc) {
        cat(line, "\n", file = stderr())
        flush(stderr())
      },
      echo_cmd = TRUE,
      env = if (in_clean_env && !use_shell_wrapper) c(processx_inherited_env_vars, new_named_env_vars) else NULL, # evaluated_env_vars e
      ...
    )

    # Check status manually
    if (result_list$status != 0) {
      format_process_stream <- function(x) {
        if (is.null(x) || !nzchar(x)) "<empty>" else x
      }

      base::stop(
        "Command failed (exit code ", result_list$status, "):\n",
        "Command:\n",
        stringr::str_flatten(c(command_string, shQuote(arguments_chr)), " "),
        "\n\nstderr:\n",
        format_process_stream(result_list$stderr),
        "\n\nstdout:\n",
        format_process_stream(result_list$stdout),
        call. = FALSE
      )
    } else {
      return(utils::tail(arguments_chr, 1)) # return last argument, which is often the output file)
    }
  }
)

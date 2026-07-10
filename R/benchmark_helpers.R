BENCHMARK_OPTIONAL_PREPROCESSING_REGEX <- c(
  "^amulet_metrics_tibble(\\.|$)",
  "^(cellsnp_dir|vireo_donor_ids_tibble)(\\.|$)"
)


benchmark_targets_match_regex <- function(target_names, regex) {
  regex <- regex[!is.na(regex) & nzchar(regex)]
  if (!length(regex)) {
    return(rep(FALSE, length(target_names)))
  }

  Reduce(`|`, purrr::map(regex, \(pattern) grepl(pattern, target_names, perl = TRUE)))
}


#' Build wall-time weights from recorded targets runtimes
#'
#' @param network Output from `targets::tar_network(targets_only = TRUE)`.
#' @param meta Output from `targets::tar_meta()`.
#' @param branch_parallel If `TRUE`, replace pattern-target runtime sums with
#'   the slowest recorded dynamic branch runtime.
#' @param exclude_target_regex Regular expressions matching target names whose
#'   recorded runtime should be set to zero without removing them from the DAG.
#' @return Tibble with one runtime weight per static target node.
#' @keywords internal
target_runtime_weights <- function(
  network,
  meta,
  branch_parallel = TRUE,
  exclude_target_regex = character()
) {
  weights <- network$vertices[, c("name", "type", "seconds", "branches"), drop = FALSE]
  weights$runtime_seconds <- as.numeric(weights$seconds)
  weights$recorded_runtime_seconds <- weights$runtime_seconds
  weights$runtime_source <- ifelse(is.na(weights$runtime_seconds), "missing", "recorded_target")

  if (isTRUE(branch_parallel)) {
    pattern_idx <- which(weights$type == "pattern")
    for (idx in pattern_idx) {
      meta_row <- meta[meta$name == weights$name[[idx]], , drop = FALSE]
      children <- if (nrow(meta_row)) unlist(meta_row$children[[1]], use.names = FALSE) else character()
      child_seconds <- meta$seconds[match(children, meta$name)]
      child_seconds <- child_seconds[!is.na(child_seconds)]

      if (length(child_seconds)) {
        weights$runtime_seconds[[idx]] <- max(child_seconds)
        weights$runtime_source[[idx]] <- "slowest_dynamic_branch"
      }
    }
  }

  weights$runtime_seconds_before_exclusion <- weights$runtime_seconds
  weights$excluded_from_benchmark <- benchmark_targets_match_regex(
    target_names = weights$name,
    regex = exclude_target_regex
  )
  weights$runtime_seconds[weights$excluded_from_benchmark] <- 0
  weights$runtime_source[weights$excluded_from_benchmark] <- "excluded"

  weights
}


benchmark_aggregation_tibble <- function(aggregations) {
  dataset_tibble_from_yaml <- read_dataset_config_tibble(config_file = "cfg_datasets.yaml")
  reaction_tibble <- build_reaction_tibble(dataset_tibble_from_yaml = dataset_tibble_from_yaml)
  aggregation_tibble_all_from_yaml <- read_aggregation_config_tibble(config_file = "cfg_aggregations.yaml")
  aggregation_tibble <- build_aggregation_tibble(
    aggregation_tibble_all_from_yaml = aggregation_tibble_all_from_yaml,
    reaction_tibble = reaction_tibble
  )

  missing_aggregations <- setdiff(aggregations, aggregation_tibble$aggregation)
  if (length(missing_aggregations)) {
    stop(
      "Aggregation(s) not found in cfg_aggregations.yaml: ",
      paste(missing_aggregations, collapse = ", "),
      call. = FALSE
    )
  }

  aggregation_tibble |>
    dplyr::mutate(
      benchmark_cellranger_count_dirs = purrr::map(
        .data$aggregation_reaction_IDs,
        \(reaction_ids) reaction_tibble$reaction_cellranger_count_dir[match(reaction_ids, reaction_tibble$reaction_ID)]
      )
    ) |>
    dplyr::slice(match(aggregations, .data$aggregation))
}


benchmark_aggregation_reaction_counts <- function(aggregations) {
  aggregation_tibble <- benchmark_aggregation_tibble(aggregations)
  reaction_counts <- lengths(aggregation_tibble$aggregation_reaction_IDs)
  stats::setNames(reaction_counts, aggregation_tibble$aggregation)[aggregations]
}


benchmark_cellranger_input_nuclei <- function(cellranger_count_dir) {
  cellranger_h5_file <- file.path(cellranger_count_dir, "outs", "filtered_feature_bc_matrix.h5")
  cellranger_h5_file_con <- hdf5r::H5File$new(cellranger_h5_file, mode = "r")
  on.exit(cellranger_h5_file_con$close_all())
  length(cellranger_h5_file_con[["matrix/barcodes"]][])
}


benchmark_aggregation_cellranger_input_nuclei <- function(aggregations) {
  aggregation_tibble <- benchmark_aggregation_tibble(aggregations)
  cellranger_input_nuclei <- purrr::map_int(
    aggregation_tibble$benchmark_cellranger_count_dirs,
    \(cellranger_count_dirs) {
      purrr::map_int(cellranger_count_dirs, benchmark_cellranger_input_nuclei) |>
        sum()
    }
  )

  stats::setNames(cellranger_input_nuclei, aggregation_tibble$aggregation)[aggregations]
}


#' Estimate wall time to a target from the recorded runtime critical path
#'
#' @param target_name Fully resolved static target name.
#' @param network Output from `targets::tar_network(targets_only = TRUE)`.
#' @param meta Output from `targets::tar_meta()`.
#' @param branch_parallel If `TRUE`, dynamic branches are assumed to run in
#'   parallel and pattern targets use their slowest branch runtime.
#' @param exclude_target_regex Regular expressions matching target names whose
#'   recorded runtime should be set to zero without removing them from the DAG.
#' @param skip_endpoint_runtime If `TRUE`, set the endpoint target runtime to
#'   zero and allow it to have missing recorded runtime. Upstream ancestors must
#'   still have recorded runtimes when `strict = TRUE`.
#' @param strict If `TRUE`, fail when any ancestor lacks recorded runtime.
#' @return List with a one-row summary tibble, the critical path, and all
#'   ancestor weights.
#' @keywords internal
estimate_target_walltime <- function(
  target_name,
  envir = parent.frame(),
  network = targets::tar_network(
    targets_only = TRUE,
    outdated = FALSE,
    reporter = "silent",
    callr_function = NULL,
    envir = envir
  ),
  meta = targets::tar_meta(),
  branch_parallel = TRUE,
  exclude_target_regex = character(),
  skip_endpoint_runtime = FALSE,
  strict = TRUE
) {
  if (!target_name %in% network$vertices$name) {
    stop("Target is not in the targets graph: ", target_name, call. = FALSE)
  }

  graph <- igraph::graph_from_data_frame(
    network$edges,
    vertices = network$vertices["name"],
    directed = TRUE
  )
  ancestor_names <- igraph::as_ids(igraph::subcomponent(graph, v = target_name, mode = "in"))
  ancestor_graph <- igraph::induced_subgraph(graph, vids = ancestor_names)
  node_order <- igraph::as_ids(igraph::topo_sort(ancestor_graph, mode = "out"))

  weights <- target_runtime_weights(
    network = network,
    meta = meta,
    branch_parallel = branch_parallel,
    exclude_target_regex = exclude_target_regex
  )
  weights <- weights[match(ancestor_names, weights$name), , drop = FALSE]
  endpoint_idx <- match(target_name, weights$name)
  weights$runtime_seconds_before_endpoint_skip <- weights$runtime_seconds
  weights$endpoint_runtime_skipped <- FALSE
  if (isTRUE(skip_endpoint_runtime)) {
    weights$runtime_seconds[[endpoint_idx]] <- 0
    weights$runtime_source[[endpoint_idx]] <- "skipped_endpoint"
    weights$endpoint_runtime_skipped[[endpoint_idx]] <- TRUE
  }
  missing_runtime <- weights$name[is.na(weights$runtime_seconds)]
  if (length(missing_runtime) && isTRUE(strict)) {
    stop(
      "Missing recorded runtime for ancestor target(s): ",
      paste(missing_runtime, collapse = ", "),
      call. = FALSE
    )
  }

  runtime <- weights$runtime_seconds
  runtime[is.na(runtime)] <- 0
  names(runtime) <- weights$name

  path_seconds <- stats::setNames(rep(NA_real_, length(node_order)), node_order)
  previous_node <- stats::setNames(rep(NA_character_, length(node_order)), node_order)

  for (node in node_order) {
    parents <- igraph::as_ids(igraph::neighbors(ancestor_graph, node, mode = "in"))
    if (!length(parents)) {
      path_seconds[[node]] <- runtime[[node]]
      next
    }

    parent_seconds <- path_seconds[parents]
    best_parent <- parents[[which.max(parent_seconds)]]
    path_seconds[[node]] <- runtime[[node]] + path_seconds[[best_parent]]
    previous_node[[node]] <- best_parent
  }

  path_names <- target_name
  while (!is.na(previous_node[[path_names[[1]]]])) {
    path_names <- c(previous_node[[path_names[[1]]]], path_names)
  }

  path <- weights[match(path_names, weights$name), , drop = FALSE]
  path$critical_path_seconds <- unname(path_seconds[path$name])
  excluded_runtime_seconds <- sum(
    weights$runtime_seconds_before_exclusion[weights$excluded_from_benchmark],
    na.rm = TRUE
  )
  skipped_endpoint_seconds <- weights$runtime_seconds_before_endpoint_skip[[endpoint_idx]]

  summary <- tibble::tibble(
    target = target_name,
    critical_path_seconds = unname(path_seconds[[target_name]]),
    critical_path_hours = unname(path_seconds[[target_name]]) / 3600,
    serial_sum_seconds = sum(runtime),
    serial_sum_hours = sum(runtime) / 3600,
    endpoint_seconds = runtime[[target_name]],
    endpoint_seconds_before_skip = skipped_endpoint_seconds,
    endpoint_runtime_skipped = isTRUE(skip_endpoint_runtime),
    ancestor_targets = length(ancestor_names),
    critical_path_targets = nrow(path),
    parallelized_pattern_targets = sum(weights$runtime_source == "slowest_dynamic_branch"),
    excluded_runtime_targets = sum(weights$excluded_from_benchmark),
    excluded_runtime_seconds = excluded_runtime_seconds,
    excluded_runtime_hours = excluded_runtime_seconds / 3600,
    missing_runtime_targets = length(missing_runtime)
  )

  list(summary = summary, critical_path = path, ancestors = weights)
}


#' Estimate wall time to aggregation-level multimodal Seurat export targets
#'
#' @param aggregations Aggregation names without the `multimodal_Seurat_object.`
#'   prefix.
#' @param store Targets store path.
#' @param branch_parallel If `TRUE`, dynamic branches are assumed to run in
#'   parallel and pattern targets use their slowest branch runtime.
#' @param aggregation_metadata Optional data frame with `aggregation`,
#'   `reaction_count`, and `cellranger_input_nuclei` columns. If `NULL`, derive
#'   these annotations from project configuration and raw Cell Ranger output.
#' @inheritParams estimate_target_walltime
#' @return Tibble with one row per aggregation. Critical-path details are stored
#'   in the `benchmark_details` attribute.
#' @keywords internal
estimate_multimodal_seurat_walltime <- function(
  aggregations,
  store = targets::tar_config_get("store"),
  branch_parallel = TRUE,
  exclude_target_regex = character(),
  skip_endpoint_runtime = FALSE,
  aggregation_metadata = NULL,
  envir = parent.frame()
) {
  network <- targets::tar_network(
    targets_only = TRUE,
    outdated = FALSE,
    reporter = "silent",
    callr_function = NULL,
    envir = envir,
    store = store
  )
  meta <- targets::tar_meta(store = store)

  estimates <- purrr::map(
    aggregations,
    \(aggregation) {
      estimate_target_walltime(
        target_name = paste0("multimodal_Seurat_object.", aggregation),
        network = network,
        meta = meta,
        branch_parallel = branch_parallel,
        exclude_target_regex = exclude_target_regex,
        skip_endpoint_runtime = skip_endpoint_runtime,
        envir = envir
      )
    }
  )
  names(estimates) <- aggregations

  summary <- purrr::map2_dfr(
    estimates,
    names(estimates),
    \(estimate, aggregation) {
      dplyr::mutate(estimate$summary, aggregation = aggregation, .before = "target")
    }
  )
  if (is.null(aggregation_metadata)) {
    reaction_counts <- benchmark_aggregation_reaction_counts(summary$aggregation)
    cellranger_input_nuclei <- benchmark_aggregation_cellranger_input_nuclei(summary$aggregation)
  } else {
    required_columns <- c("aggregation", "reaction_count", "cellranger_input_nuclei")
    missing_columns <- setdiff(required_columns, names(aggregation_metadata))
    if (length(missing_columns)) {
      stop(
        "aggregation_metadata is missing required column(s): ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    if (anyDuplicated(aggregation_metadata$aggregation)) {
      stop("aggregation_metadata contains duplicate aggregation names.", call. = FALSE)
    }
    aggregation_metadata <- aggregation_metadata[
      match(summary$aggregation, aggregation_metadata$aggregation),
      required_columns,
      drop = FALSE
    ]
    if (anyNA(aggregation_metadata$aggregation)) {
      stop("aggregation_metadata does not cover every requested aggregation.", call. = FALSE)
    }
    if (
      !is.numeric(aggregation_metadata$reaction_count) ||
        !is.numeric(aggregation_metadata$cellranger_input_nuclei) ||
        anyNA(aggregation_metadata[c("reaction_count", "cellranger_input_nuclei")]) ||
        any(aggregation_metadata$reaction_count <= 0) ||
        any(aggregation_metadata$cellranger_input_nuclei <= 0)
    ) {
      stop(
        "aggregation_metadata reaction and nuclei counts must be positive numeric values.",
        call. = FALSE
      )
    }
    reaction_counts <- stats::setNames(
      aggregation_metadata$reaction_count,
      aggregation_metadata$aggregation
    )
    cellranger_input_nuclei <- stats::setNames(
      aggregation_metadata$cellranger_input_nuclei,
      aggregation_metadata$aggregation
    )
  }
  summary <- dplyr::mutate(
    summary,
    reaction_count = unname(reaction_counts[.data$aggregation]),
    cellranger_input_nuclei = unname(cellranger_input_nuclei[.data$aggregation]),
    critical_path_minutes = .data$critical_path_seconds / 60,
    serial_sum_minutes = .data$serial_sum_seconds / 60,
    excluded_runtime_minutes = .data$excluded_runtime_seconds / 60,
    .after = "target"
  )
  attr(summary, "benchmark_details") <- estimates
  attr(summary, "benchmark_aggregations") <- aggregations
  attr(summary, "benchmark_branch_parallel") <- branch_parallel
  attr(summary, "benchmark_exclude_target_regex") <- exclude_target_regex
  attr(summary, "benchmark_skip_endpoint_runtime") <- skip_endpoint_runtime
  attr(summary, "benchmark_aggregation_metadata") <- summary[, c(
    "aggregation",
    "reaction_count",
    "cellranger_input_nuclei"
  )]
  attr(summary, "benchmark_target_label_suffixes") <- benchmark_plot_label_suffixes(summary)
  summary
}


#' Estimate and cache multimodal Seurat wall-time benchmarks
#'
#' @param aggregations Aggregation names without the `multimodal_Seurat_object.`
#'   prefix.
#' @param cache_file RDS file used to cache the benchmark result object. If
#'   `NULL`, cache below `store/benchmark/`.
#' @param force If `TRUE`, recompute the estimates even when `cache_file` exists.
#' @inheritParams estimate_multimodal_seurat_walltime
#' @return Benchmark summary tibble with `benchmark_details` attribute.
#' @keywords internal
cache_multimodal_seurat_walltime <- function(
  aggregations,
  cache_file = NULL,
  force = FALSE,
  store = targets::tar_config_get("store"),
  branch_parallel = TRUE,
  exclude_target_regex = character(),
  skip_endpoint_runtime = FALSE,
  aggregation_metadata = NULL,
  envir = parent.frame()
) {
  if (is.null(cache_file)) {
    cache_file <- file.path(store, "benchmark", "multimodal_seurat_walltime.rds")
  }

  if (file.exists(cache_file) && !isTRUE(force)) {
    summary <- readRDS(cache_file)
    cache_is_current <- "cellranger_input_nuclei" %in% colnames(summary) &&
      identical(attr(summary, "benchmark_aggregations"), aggregations) &&
      identical(attr(summary, "benchmark_branch_parallel"), branch_parallel) &&
      identical(attr(summary, "benchmark_exclude_target_regex"), exclude_target_regex) &&
      identical(attr(summary, "benchmark_skip_endpoint_runtime"), skip_endpoint_runtime)

    if (isTRUE(cache_is_current)) {
      return(summary)
    }
  }

  summary <- estimate_multimodal_seurat_walltime(
    aggregations = aggregations,
    store = store,
    branch_parallel = branch_parallel,
    exclude_target_regex = exclude_target_regex,
    skip_endpoint_runtime = skip_endpoint_runtime,
    aggregation_metadata = aggregation_metadata,
    envir = envir
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(summary, cache_file)
  summary
}


benchmark_plot_label_suffixes <- function(benchmark_results) {
  cached_suffixes <- attr(benchmark_results, "benchmark_target_label_suffixes")
  if (length(cached_suffixes)) {
    return(cached_suffixes)
  }

  reaction_IDs <- tryCatch(
    {
      aggregation_tibble <- benchmark_aggregation_tibble(benchmark_results$aggregation)
      unlist(aggregation_tibble$aggregation_reaction_IDs, use.names = FALSE)
    },
    error = \(condition) character()
  )

  unique(c(benchmark_results$aggregation, reaction_IDs))
}


benchmark_strip_target_label_suffixes <- function(target_names, suffixes) {
  suffixes <- suffixes[order(nchar(suffixes), decreasing = TRUE)]
  purrr::map_chr(
    target_names,
    \(target_name) {
      for (suffix in suffixes) {
        target_name <- sub(paste0("\\.", suffix, "$"), "", target_name)
      }
      target_name
    }
  )
}


benchmark_walltime_composition_tibble <- function(benchmark_results, top_n_targets = 10) {
  benchmark_details <- attr(benchmark_results, "benchmark_details")
  if (is.null(benchmark_details)) {
    stop(
      "Benchmark result is missing the 'benchmark_details' attribute. ",
      "Recompute it with estimate_multimodal_seurat_walltime() or cache_multimodal_seurat_walltime().",
      call. = FALSE
    )
  }

  missing_details <- setdiff(benchmark_results$aggregation, names(benchmark_details))
  if (length(missing_details)) {
    stop(
      "Benchmark details are missing aggregation(s): ",
      paste(missing_details, collapse = ", "),
      call. = FALSE
    )
  }

  suffixes <- benchmark_plot_label_suffixes(benchmark_results)
  path_tibble <- purrr::map_dfr(
    benchmark_results$aggregation,
    \(aggregation) {
      benchmark_details[[aggregation]]$critical_path |>
        tibble::as_tibble() |>
        dplyr::mutate(aggregation = aggregation, .before = 1)
    }
  ) |>
    dplyr::filter(.data$runtime_seconds > 0) |>
    dplyr::mutate(
      target_step = benchmark_strip_target_label_suffixes(.data$name, suffixes),
      runtime_minutes = .data$runtime_seconds / 60
    )

  top_targets <- path_tibble |>
    dplyr::summarise(runtime_minutes = sum(.data$runtime_minutes), .by = "target_step") |>
    dplyr::slice_max(.data$runtime_minutes, n = top_n_targets, with_ties = FALSE) |>
    dplyr::arrange(dplyr::desc(.data$runtime_minutes)) |>
    dplyr::pull(.data$target_step)
  target_levels <- c("Other", top_targets)

  path_tibble |>
    dplyr::mutate(
      target_step = dplyr::if_else(.data$target_step %in% top_targets, .data$target_step, "Other"),
      target_step = factor(.data$target_step, levels = target_levels)
    ) |>
    dplyr::summarise(runtime_minutes = sum(.data$runtime_minutes), .by = c("aggregation", "target_step")) |>
    dplyr::left_join(
      benchmark_results |>
        dplyr::select("aggregation", "cellranger_input_nuclei", "reaction_count", "critical_path_minutes"),
      by = "aggregation"
    )
}


#' Plot multimodal Seurat wall time by CellRanger input nuclei
#'
#' @param benchmark_results Output from `estimate_multimodal_seurat_walltime()`
#'   or `cache_multimodal_seurat_walltime()`.
#' @param time_col Numeric column to plot on the y axis.
#' @param top_n_targets Number of critical-path target steps to show as
#'   separate bar segments. Remaining steps are collapsed to `Other`.
#' @param target_step_labels Optional named character vector mapping target-step
#'   names to reader-facing labels.
#' @param color_critical_path_steps If `TRUE`, color and show a legend for
#'   critical-path bar segments. If `FALSE`, use a single neutral fill and hide
#'   the fill legend.
#' @return ggplot object.
#' @keywords internal
plot_multimodal_seurat_walltime <- function(
  benchmark_results,
  time_col = "critical_path_minutes",
  top_n_targets = 10,
  target_step_labels = NULL,
  color_critical_path_steps = TRUE
) {
  if (!identical(time_col, "critical_path_minutes")) {
    stop("The composition plot only supports time_col = 'critical_path_minutes'.", call. = FALSE)
  }

  required_cols <- c("aggregation", "reaction_count", "cellranger_input_nuclei", time_col)
  missing_cols <- setdiff(required_cols, names(benchmark_results))
  if (length(missing_cols)) {
    stop(
      "Benchmark result is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  composition_tibble <- benchmark_walltime_composition_tibble(
    benchmark_results = benchmark_results,
    top_n_targets = top_n_targets
  )
  plot_tibble <- benchmark_results |>
    dplyr::arrange(.data$cellranger_input_nuclei) |>
    dplyr::mutate(
      walltime_hours = .data[[time_col]] / 60
    )
  composition_tibble <- composition_tibble |>
    dplyr::mutate(runtime_hours = .data$runtime_minutes / 60)
  unique_x <- sort(unique(plot_tibble$cellranger_input_nuclei))
  if (any(unique_x <= 0)) {
    stop("cellranger_input_nuclei must be positive for a log10 x axis.", call. = FALSE)
  }
  y_plot_limit <- 10
  linear_anchor_x <- unique_x[[1]]
  linear_anchor_y <- plot_tibble$walltime_hours[[1]]
  if (linear_anchor_y <= 0) {
    stop("The smallest benchmark wall time must be positive for the linear reference.", call. = FALSE)
  }
  x_axis_anchor <- 10000
  linear_reference_tibble <- tibble::tibble(
    cellranger_input_nuclei = 10^seq(
      log10(x_axis_anchor),
      log10(max(unique_x)),
      length.out = 200
    )
  ) |>
    dplyr::mutate(
      walltime_hours = linear_anchor_y * .data$cellranger_input_nuclei / linear_anchor_x
    )
  linear_cap_x <- linear_anchor_x * y_plot_limit / linear_anchor_y
  if (max(linear_reference_tibble$walltime_hours) > y_plot_limit) {
    linear_reference_tibble <- linear_reference_tibble |>
      dplyr::filter(.data$walltime_hours <= y_plot_limit) |>
      dplyr::bind_rows(tibble::tibble(
        cellranger_input_nuclei = linear_cap_x,
        walltime_hours = y_plot_limit
      )) |>
      dplyr::arrange(.data$cellranger_input_nuclei)
    linear_arrow_tibble <- tibble::tibble(
      x = linear_anchor_x * y_plot_limit * 0.88 / linear_anchor_y,
      xend = linear_cap_x,
      y = y_plot_limit * 0.88,
      yend = y_plot_limit
    )
    linear_label_tibble <- tibble::tibble(
      cellranger_input_nuclei = linear_cap_x,
      walltime_hours = y_plot_limit,
      label = "Linear reference"
    )
  } else {
    linear_arrow_tibble <- tibble::tibble(x = numeric(), xend = numeric(), y = numeric(), yend = numeric())
    linear_label_tibble <- dplyr::slice_tail(linear_reference_tibble, n = 1) |>
      dplyr::mutate(label = "Linear reference")
  }
  bar_width_log10 <- if (length(unique_x) > 1) {
    min(diff(log10(unique_x))) * 0.55
  } else {
    0.2
  }
  composition_plot_tibble <- composition_tibble |>
    dplyr::arrange(.data$cellranger_input_nuclei, .data$target_step) |>
    dplyr::mutate(
      y_max = cumsum(.data$runtime_hours),
      y_min = .data$y_max - .data$runtime_hours,
      x_min = 10^(log10(.data$cellranger_input_nuclei) - bar_width_log10 / 2),
      x_max = 10^(log10(.data$cellranger_input_nuclei) + bar_width_log10 / 2),
      target_fill = if (isTRUE(color_critical_path_steps)) {
        as.character(.data$target_step)
      } else {
        "Critical path"
      },
      .by = "aggregation"
    )
  target_levels <- if (isTRUE(color_critical_path_steps)) {
    levels(composition_tibble$target_step)
  } else {
    "Critical path"
  }
  target_levels_without_other <- setdiff(target_levels, "Other")
  fill_values <- if (isTRUE(color_critical_path_steps)) {
    c(
      Other = "grey70",
      stats::setNames(scales::hue_pal()(length(target_levels_without_other)), target_levels_without_other)
    )
  } else {
    c("Critical path" = "grey70")
  }
  fill_labels <- function(labels) {
    if (!is.null(target_step_labels)) {
      mapped_labels <- target_step_labels[labels]
      labels[!is.na(mapped_labels)] <- unname(mapped_labels[!is.na(mapped_labels)])
    }
    stringr::str_wrap(labels, width = 34)
  }

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = composition_plot_tibble,
      ggplot2::aes(
        xmin = .data$x_min,
        xmax = .data$x_max,
        ymin = .data$y_min,
        ymax = .data$y_max,
        fill = .data$target_fill
      ),
      color = if (isTRUE(color_critical_path_steps)) "white" else NA,
      linewidth = if (isTRUE(color_critical_path_steps)) 0.2 else 0
    ) +
    ggplot2::geom_segment(
      data = plot_tibble,
      ggplot2::aes(
        x = .data$cellranger_input_nuclei,
        xend = .data$cellranger_input_nuclei,
        y = 0,
        yend = .data$walltime_hours
      ),
      color = "grey15",
      linetype = "dashed",
      linewidth = 0.35,
      alpha = 0.65
    ) +
    ggplot2::geom_line(
      data = linear_reference_tibble,
      ggplot2::aes(x = .data$cellranger_input_nuclei, y = .data$walltime_hours),
      linewidth = 0.65,
      linetype = "dotdash",
      color = "#8B1E3F"
    ) +
    ggplot2::geom_segment(
      data = linear_arrow_tibble,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed"),
      linewidth = 0.65,
      linetype = "dotdash",
      color = "#8B1E3F"
    ) +
    ggplot2::geom_text(
      data = linear_label_tibble,
      ggplot2::aes(x = .data$cellranger_input_nuclei, y = .data$walltime_hours, label = .data$label),
      hjust = 0,
      vjust = 1.15,
      size = 3,
      color = "#8B1E3F"
    ) +
    ggplot2::scale_x_log10(
      breaks = c(x_axis_anchor, unique_x),
      labels = \(breaks) {
        scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale())(breaks)
      },
      limits = c(x_axis_anchor, NA),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::annotation_logticks(sides = "b") +
    ggplot2::scale_fill_manual(
      values = fill_values[target_levels],
      labels = fill_labels,
      guide = if (isTRUE(color_critical_path_steps)) "legend" else "none"
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, y_plot_limit, by = 1),
      labels = scales::label_number(accuracy = 1),
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(
      x = "CellRanger input nuclei",
      y = "Wall time (hours)",
      fill = "Critical-path step"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_plot_limit), clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey88", linewidth = 0.35),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.88),
        color = "grey80",
        linewidth = 0.25
      ),
      axis.title = ggplot2::element_text(color = "grey15"),
      axis.text = ggplot2::element_text(color = "grey30")
    )
}

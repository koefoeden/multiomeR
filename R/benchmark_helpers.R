#' Build wall-time weights from recorded targets runtimes
#'
#' @param network Output from `targets::tar_network(targets_only = TRUE)`.
#' @param meta Output from `targets::tar_meta()`.
#' @param branch_parallel If `TRUE`, replace pattern-target runtime sums with
#'   the slowest recorded dynamic branch runtime.
#' @return Tibble with one runtime weight per static target node.
#' @keywords internal
target_runtime_weights <- function(network, meta, branch_parallel = TRUE) {
  weights <- network$vertices[, c("name", "type", "seconds", "branches"), drop = FALSE]
  weights$runtime_seconds <- as.numeric(weights$seconds)
  weights$runtime_source <- ifelse(is.na(weights$runtime_seconds), "missing", "recorded_target")

  if (!isTRUE(branch_parallel)) {
    return(weights)
  }

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

  weights
}


benchmark_aggregation_reaction_counts <- function(aggregations) {
  dataset_tibble_from_yaml <- read_config_tibble(config_file = "cfg_datasets.yaml", key_col = "dataset")
  reaction_tibble <- build_reaction_tibble(dataset_tibble_from_yaml = dataset_tibble_from_yaml)
  aggregation_tibble_all_from_yaml <- read_config_tibble(config_file = "cfg_aggregations.yaml", key_col = "aggregation")
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

  reaction_counts <- lengths(aggregation_tibble$aggregation_reaction_IDs)
  stats::setNames(reaction_counts, aggregation_tibble$aggregation)[aggregations]
}


#' Estimate wall time to a target from the recorded runtime critical path
#'
#' @param target_name Fully resolved static target name.
#' @param network Output from `targets::tar_network(targets_only = TRUE)`.
#' @param meta Output from `targets::tar_meta()`.
#' @param branch_parallel If `TRUE`, dynamic branches are assumed to run in
#'   parallel and pattern targets use their slowest branch runtime.
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
    branch_parallel = branch_parallel
  )
  weights <- weights[match(ancestor_names, weights$name), , drop = FALSE]
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

  summary <- tibble::tibble(
    target = target_name,
    critical_path_seconds = unname(path_seconds[[target_name]]),
    critical_path_hours = unname(path_seconds[[target_name]]) / 3600,
    serial_sum_seconds = sum(runtime),
    serial_sum_hours = sum(runtime) / 3600,
    endpoint_seconds = runtime[[target_name]],
    ancestor_targets = length(ancestor_names),
    critical_path_targets = nrow(path),
    parallelized_pattern_targets = sum(weights$runtime_source == "slowest_dynamic_branch"),
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
#' @return Tibble with one row per aggregation. Critical-path details are stored
#'   in the `benchmark_details` attribute.
#' @keywords internal
estimate_multimodal_seurat_walltime <- function(
  aggregations,
  store = targets::tar_config_get("store"),
  branch_parallel = TRUE,
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
  reaction_counts <- benchmark_aggregation_reaction_counts(summary$aggregation)
  summary <- dplyr::mutate(
    summary,
    reaction_count = unname(reaction_counts[.data$aggregation]),
    critical_path_minutes = .data$critical_path_seconds / 60,
    serial_sum_minutes = .data$serial_sum_seconds / 60,
    .after = "target"
  )
  attr(summary, "benchmark_details") <- estimates
  summary
}


#' Estimate and cache multimodal Seurat wall-time benchmarks
#'
#' @param aggregations Aggregation names without the `multimodal_Seurat_object.`
#'   prefix.
#' @param cache_file RDS file used to cache the benchmark result object.
#' @param force If `TRUE`, recompute the estimates even when `cache_file` exists.
#' @inheritParams estimate_multimodal_seurat_walltime
#' @return Benchmark summary tibble with `benchmark_details` attribute.
#' @keywords internal
cache_multimodal_seurat_walltime <- function(
  aggregations,
  cache_file = file.path("outputs", "benchmark", "multimodal_seurat_walltime.rds"),
  force = FALSE,
  store = targets::tar_config_get("store"),
  branch_parallel = TRUE,
  envir = parent.frame()
) {
  if (file.exists(cache_file) && !isTRUE(force)) {
    return(readRDS(cache_file))
  }

  summary <- estimate_multimodal_seurat_walltime(
    aggregations = aggregations,
    store = store,
    branch_parallel = branch_parallel,
    envir = envir
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(summary, cache_file)
  summary
}


#' Plot multimodal Seurat wall time by reaction count
#'
#' @param benchmark_results Output from `estimate_multimodal_seurat_walltime()`
#'   or `cache_multimodal_seurat_walltime()`.
#' @param time_col Numeric column to plot on the y axis.
#' @return ggplot object.
#' @keywords internal
plot_multimodal_seurat_walltime <- function(
  benchmark_results,
  time_col = "critical_path_minutes"
) {
  required_cols <- c("aggregation", "reaction_count", time_col)
  missing_cols <- setdiff(required_cols, names(benchmark_results))
  if (length(missing_cols)) {
    stop(
      "Benchmark result is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  plot_tibble <- benchmark_results |>
    dplyr::arrange(.data$reaction_count)

  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(x = .data$reaction_count, y = .data[[time_col]])
  ) +
    ggplot2::geom_line(linewidth = 0.7, color = "#3D5A5B") +
    ggplot2::geom_point(size = 2.8, color = "#1F2D2E") +
    ggplot2::scale_x_continuous(breaks = plot_tibble$reaction_count) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(accuracy = 0.1),
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(
      x = "Reactions",
      y = "Wall time (minutes)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey88", linewidth = 0.35),
      axis.title = ggplot2::element_text(color = "grey15"),
      axis.text = ggplot2::element_text(color = "grey30")
    )
}

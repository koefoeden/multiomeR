#' Get cache paths for a rendered targets graph
#'
#' @param pipeline_name Pipeline or module name used in the cache key.
#' @param target_graph_match Regex pattern passed to `tidyselect::matches()`.
#' @param cache_dir Directory for cached graph files.
#' @return Named list of graph cache paths.
#' @keywords internal
get_targets_graph_cache_paths <- function(pipeline_name, target_graph_match, cache_dir = file.path("website", "cache", "targets_graphs")) {
  cache_key <- digest::digest(
    list(
      pipeline_name = pipeline_name,
      target_graph_match = target_graph_match
    ),
    algo = "xxhash64"
  )
  file_stem <- paste0(pipeline_name, "_", cache_key)

  list(
    cache_dir = cache_dir,
    file_stem = file_stem,
    glimpse_object = file.path(cache_dir, paste0(file_stem, "_glimpse_object.qs2")),
    mermaid = file.path(cache_dir, paste0(file_stem, "_mermaid.mmd")),
    nodes = file.path(cache_dir, paste0(file_stem, "_targets_nodes.csv")),
    edges = file.path(cache_dir, paste0(file_stem, "_targets_graph.csv"))
  )
}

targets_graph_target_family <- function(name) {
  sub("[.].*$", "", name)
}

targets_graph_degree_table <- function(nodes, edges) {
  node_names <- nodes$name
  in_degree <- tabulate(match(edges$to, node_names), nbins = length(node_names))
  out_degree <- tabulate(match(edges$from, node_names), nbins = length(node_names))

  data.frame(
    name = node_names,
    target_family = targets_graph_target_family(node_names),
    in_degree = in_degree,
    out_degree = out_degree,
    total_degree = in_degree + out_degree,
    is_source_only = out_degree > 0L & in_degree == 0L,
    is_sink_only = in_degree > 0L & out_degree == 0L,
    is_disconnected = in_degree == 0L & out_degree == 0L,
    stringsAsFactors = FALSE
  )
}

targets_graph_node_table <- function(glimpse_graph) {
  nodes <- glimpse_graph$x$nodes
  edges <- glimpse_graph$x$edges
  degree_table <- targets_graph_degree_table(nodes, edges)
  metrics <- degree_table[
    match(nodes$name, degree_table$name),
    setdiff(names(degree_table), "name"),
    drop = FALSE
  ]
  row.names(metrics) <- NULL
  cbind(nodes, metrics)
}

targets_graph_prefixed_node_metrics <- function(degree_table, node_names, prefix) {
  metrics <- degree_table[
    match(node_names, degree_table$name),
    setdiff(names(degree_table), "name"),
    drop = FALSE
  ]
  names(metrics) <- paste0(prefix, names(metrics))
  row.names(metrics) <- NULL
  metrics
}

targets_graph_edge_table <- function(glimpse_graph) {
  edges <- glimpse_graph$x$edges
  degree_table <- targets_graph_degree_table(glimpse_graph$x$nodes, edges)

  from_metrics <- targets_graph_prefixed_node_metrics(degree_table, edges$from, "from_")
  to_metrics <- targets_graph_prefixed_node_metrics(degree_table, edges$to, "to_")
  edges <- cbind(edges, from_metrics, to_metrics)
  edges$crosses_target_family <- edges$from_target_family != edges$to_target_family
  edges
}

targets_graph_mermaid_shape_open <- function(type) {
  unname(c(
    object = "{{",
    "function" = ">",
    stem = "([",
    pattern = "["
  )[type])
}

targets_graph_mermaid_shape_close <- function(type) {
  unname(c(
    object = "}}",
    "function" = "]",
    stem = "])",
    pattern = "]"
  )[type])
}

targets_graph_mermaid_vertex_text <- function(nodes) {
  sprintf(
    "%s%s\"%s\"%s:::%s",
    nodes$name,
    targets_graph_mermaid_shape_open(nodes$type),
    nodes$label,
    targets_graph_mermaid_shape_close(nodes$type),
    nodes$status
  )
}

targets_graph_mermaid_lines <- function(glimpse_graph) {
  nodes <- glimpse_graph$x$nodes
  if (nrow(nodes) < 1L) {
    return("")
  }

  text <- targets_graph_mermaid_vertex_text(nodes)
  names(text) <- nodes$name

  edges <- glimpse_graph$x$edges
  disconnected <- setdiff(nodes$name, c(edges$from, edges$to))
  disconnected <- unname(text[disconnected])
  disconnected <- if (length(disconnected)) paste0("    ", disconnected) else character(0L)
  edges$from <- unname(text[edges$from])
  edges$to <- unname(text[edges$to])

  c(
    "%%{init: {'flowchart': {'defaultRenderer': 'elk'}}}%%",
    "flowchart TB",
    "  style Graph fill:#FFFFFF00,stroke:#000000;",
    "  subgraph Graph",
    "    direction TB",
    sprintf("    %s --> %s", edges$from, edges$to),
    disconnected,
    "  end"
  )
}

targets_graph_mermaid_cache_is_current <- function(lines) {
  identical(
    lines[seq_len(min(length(lines), 2L))],
    c("%%{init: {'flowchart': {'defaultRenderer': 'elk'}}}%%", "flowchart TB")
  )
}

#' Build cached targets graph assets
#'
#' Creates the same graph assets previously generated directly inside the
#' Quarto target-graph snippet: a `tar_glimpse()` object, a Mermaid graph, and
#' node and edge tables enriched with graph metrics. This helper is intentionally
#' separate from the docs so future output-gallery or graph pages can reuse it.
#'
#' @param pipeline_name Pipeline or module name used in the cache key.
#' @param target_graph_match Regex pattern passed to `tidyselect::matches()`.
#' @param cache_dir Directory for cached graph files.
#' @param force If `TRUE`, recreate cached assets even when files already exist.
#' @param allow_match Regex pattern for the `tar_glimpse()` allow argument.
#' @param formatter Optional function applied to the `tar_glimpse()` object.
#' @return Named list with graph assets and cache paths.
#' @keywords internal
build_targets_graph_assets <- function(
  pipeline_name,
  target_graph_match,
  cache_dir = file.path("website", "cache", "targets_graphs"),
  force = FALSE,
  allow_match = target_graph_match,
  formatter = NULL
) {
  paths <- get_targets_graph_cache_paths(
    pipeline_name = pipeline_name,
    target_graph_match = target_graph_match,
    cache_dir = cache_dir
  )
  dir.create(paths$cache_dir, recursive = TRUE, showWarnings = FALSE)

  glimpse_graph <- if (file.exists(paths$glimpse_object) && !force) {
    qs2::qs_read(paths$glimpse_object)
  } else {
    graph <- rlang::inject(
      targets::tar_glimpse(
        names = tidyselect::matches(!!target_graph_match),
        allow = tidyselect::matches(!!allow_match),
        shortcut = TRUE,
        physics = TRUE
      )
    )
    qs2::qs_save(graph, paths$glimpse_object)
    graph
  }

  mermaid_lines <- if (file.exists(paths$mermaid) && !force) {
    lines <- readLines(paths$mermaid, warn = FALSE)
    if (targets_graph_mermaid_cache_is_current(lines)) {
      lines
    } else {
      lines <- targets_graph_mermaid_lines(glimpse_graph)
      writeLines(lines, paths$mermaid)
      lines
    }
  } else {
    lines <- targets_graph_mermaid_lines(glimpse_graph)
    writeLines(lines, paths$mermaid)
    lines
  }

  nodes <- targets_graph_node_table(glimpse_graph)
  edges <- targets_graph_edge_table(glimpse_graph)
  utils::write.csv(nodes, paths$nodes, row.names = FALSE)
  if (!is.null(edges) && all(c("from", "to") %in% names(edges))) {
    utils::write.csv(edges, paths$edges, row.names = FALSE)
  }

  list(
    glimpse_graph = glimpse_graph,
    mermaid_lines = mermaid_lines,
    nodes = nodes,
    edges = edges,
    formatted_graph = if (is.null(formatter)) NULL else formatter(glimpse_graph),
    paths = paths
  )
}

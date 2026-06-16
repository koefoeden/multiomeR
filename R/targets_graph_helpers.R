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
  disconnected <- if (length(disconnected)) paste0("  ", disconnected) else character(0L)
  edges$from <- unname(text[edges$from])
  edges$to <- unname(text[edges$to])

  c(
    "%%{init: {'flowchart': {'defaultRenderer': 'elk'}}}%%",
    "flowchart TB",
    sprintf("  %s --> %s", edges$from, edges$to),
    disconnected
  )
}

targets_graph_mermaid_cache_is_current <- function(lines) {
  identical(
    lines[seq_len(min(length(lines), 2L))],
    c("%%{init: {'flowchart': {'defaultRenderer': 'elk'}}}%%", "flowchart TB")
  )
}

targets_graph_part_of_graph_tag <- function(graph_id) {
  if (!grepl("^[A-Za-z0-9_]+$", graph_id)) {
    stop("graph_id must contain only letters, numbers, and underscores: ", graph_id, call. = FALSE)
  }
  paste0("[part_of_graph:", graph_id, "]")
}

targets_graph_part_of_graph_ids <- function(manifest = targets::tar_manifest(callr_function = NULL)) {
  descriptions <- as.character(manifest$description)
  descriptions[is.na(descriptions)] <- ""
  tags <- regmatches(
    descriptions,
    gregexpr("\\[part_of_graph:[A-Za-z0-9_]+\\]", descriptions)
  )
  tags <- unlist(tags, use.names = FALSE)
  sort(unique(sub("\\]$", "", sub("^\\[part_of_graph:", "", tags))))
}

targets_graph_part_of_graph_names <- function(graph_id, manifest = targets::tar_manifest(callr_function = NULL)) {
  tag <- targets_graph_part_of_graph_tag(graph_id)
  descriptions <- as.character(manifest$description)
  descriptions[is.na(descriptions)] <- ""
  manifest$name[grepl(tag, descriptions, fixed = TRUE)]
}

targets_graph_default_label_suffixes <- function() {
  suffixes <- character()
  if (exists("reaction_tibble", inherits = TRUE)) {
    reaction_tibble_obj <- get("reaction_tibble", inherits = TRUE)
    if ("reaction_ID" %in% names(reaction_tibble_obj)) {
      suffixes <- c(suffixes, reaction_tibble_obj$reaction_ID)
    }
  }
  if (exists("dataset_tibble", inherits = TRUE)) {
    dataset_tibble_obj <- get("dataset_tibble", inherits = TRUE)
    if ("dataset" %in% names(dataset_tibble_obj)) {
      suffixes <- c(suffixes, dataset_tibble_obj$dataset)
    }
  }
  if (exists("aggregation_tibble", inherits = TRUE)) {
    aggregation_tibble_obj <- get("aggregation_tibble", inherits = TRUE)
    if ("aggregation" %in% names(aggregation_tibble_obj)) {
      suffixes <- c(suffixes, aggregation_tibble_obj$aggregation)
    }
  }
  if (exists("known_aggregation_modules", inherits = TRUE)) {
    suffixes <- c(suffixes, get("known_aggregation_modules", inherits = TRUE))
  }
  suffixes <- unique(as.character(suffixes))
  suffixes[!is.na(suffixes) & nzchar(suffixes)]
}

targets_graph_regex_escape <- function(x) {
  gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", x, perl = TRUE)
}

targets_graph_clean_node_labels <- function(glimpse_graph, suffixes = targets_graph_default_label_suffixes()) {
  suffixes <- unique(as.character(suffixes))
  suffixes <- suffixes[!is.na(suffixes) & nzchar(suffixes)]
  if (!length(suffixes)) {
    return(glimpse_graph)
  }

  suffixes <- suffixes[order(nchar(suffixes), decreasing = TRUE)]
  pattern <- paste0("\\.(", paste(targets_graph_regex_escape(suffixes), collapse = "|"), ")(?=\\.|$)")
  glimpse_graph$x$nodes$label <- gsub(pattern, "", glimpse_graph$x$nodes$label, perl = TRUE)
  glimpse_graph
}

targets_graph_adjacency_get <- function(adjacency, node) {
  value <- adjacency[[node]]
  if (is.null(value)) {
    character()
  } else {
    unname(value)
  }
}

targets_graph_nearest_retained <- function(start, adjacency, removed) {
  retained <- character()
  seen <- start
  queue <- targets_graph_adjacency_get(adjacency, start)

  while (length(queue)) {
    current <- queue[[1]]
    queue <- queue[-1]
    if (current %in% seen) {
      next
    }
    seen <- c(seen, current)

    if (current %in% removed) {
      queue <- c(queue, targets_graph_adjacency_get(adjacency, current))
    } else {
      retained <- c(retained, current)
    }
  }

  unique(retained)
}

targets_graph_bypass_removed_edges <- function(edges, removed) {
  if (!nrow(edges)) {
    return(edges)
  }

  kept_edges <- edges[
    !edges$from %in% removed & !edges$to %in% removed,
    ,
    drop = FALSE
  ]

  forward <- split(as.character(edges$to), as.character(edges$from))
  reverse <- split(as.character(edges$from), as.character(edges$to))
  bypass_pairs <- data.frame(from = character(), to = character())

  for (node in removed) {
    parents <- targets_graph_nearest_retained(node, reverse, removed)
    children <- targets_graph_nearest_retained(node, forward, removed)
    if (length(parents) && length(children)) {
      bypass_pairs <- rbind(
        bypass_pairs,
        expand.grid(from = parents, to = children, stringsAsFactors = FALSE)
      )
    }
  }

  bypass_pairs <- unique(bypass_pairs[bypass_pairs$from != bypass_pairs$to, , drop = FALSE])
  kept_pairs <- paste(kept_edges$from, kept_edges$to, sep = "\r")
  bypass_pairs <- bypass_pairs[!paste(bypass_pairs$from, bypass_pairs$to, sep = "\r") %in% kept_pairs, , drop = FALSE]
  if (!nrow(bypass_pairs)) {
    return(kept_edges)
  }

  bypass_edges <- bypass_pairs
  for (column in setdiff(names(edges), c("from", "to"))) {
    bypass_edges[[column]] <- edges[[column]][[1]]
  }
  bypass_edges <- bypass_edges[, names(edges), drop = FALSE]
  unique(rbind(kept_edges, bypass_edges))
}

targets_graph_prune_to_part_of_graph <- function(glimpse_graph, graph_id) {
  tag <- targets_graph_part_of_graph_tag(graph_id)
  nodes <- glimpse_graph$x$nodes
  descriptions <- as.character(nodes$description)
  descriptions[is.na(descriptions)] <- ""
  keep_node <- grepl(tag, descriptions, fixed = TRUE)
  if (!any(keep_node)) {
    stop("No nodes in the graph have tag ", tag, call. = FALSE)
  }

  removed <- nodes$name[!keep_node]
  glimpse_graph$x$nodes <- nodes[keep_node, , drop = FALSE]
  glimpse_graph$x$edges <- targets_graph_bypass_removed_edges(glimpse_graph$x$edges, removed)
  attr(glimpse_graph, "targets_graph_prune_summary") <- list(
    graph_id = graph_id,
    kept_nodes = sum(keep_node),
    removed_nodes = length(removed),
    edges = nrow(glimpse_graph$x$edges)
  )
  glimpse_graph
}

targets_graph_part_of_graph_network <- function(
  graph_id,
  manifest = targets::tar_manifest(callr_function = NULL),
  shortcut = FALSE,
  physics = FALSE
) {
  target_names <- targets_graph_part_of_graph_names(graph_id, manifest = manifest)
  if (!length(target_names)) {
    stop("No targets in the manifest have tag ", targets_graph_part_of_graph_tag(graph_id), call. = FALSE)
  }

  rlang::inject(
    targets::tar_glimpse(
      names = tidyselect::any_of(!!target_names),
      shortcut = !!shortcut,
      physics = !!physics
    )
  )
}

targets_graph_write_part_of_graph_mermaid <- function(
  graph_id,
  output_file,
  manifest = targets::tar_manifest(callr_function = NULL),
  label_suffixes = targets_graph_default_label_suffixes(),
  shortcut = FALSE,
  physics = FALSE
) {
  graph <- targets_graph_part_of_graph_network(
    graph_id = graph_id,
    manifest = manifest,
    shortcut = shortcut,
    physics = physics
  ) |>
    targets_graph_prune_to_part_of_graph(graph_id = graph_id) |>
    targets_graph_clean_node_labels(suffixes = label_suffixes)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(targets_graph_mermaid_lines(graph), output_file)
  invisible(graph)
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

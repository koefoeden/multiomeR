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

targets_graph_mermaid_label_text <- function(label) {
  label <- gsub("<", "&lt;", label, fixed = TRUE)
  gsub(">", "&gt;", label, fixed = TRUE)
}

targets_graph_mermaid_vertex_text <- function(nodes) {
  sprintf(
    "%s%s\"%s\"%s:::%s",
    nodes$name,
    targets_graph_mermaid_shape_open(nodes$type),
    targets_graph_mermaid_label_text(nodes$label),
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
  suffixes <- data.frame(suffix = character(), placeholder = character())
  if (exists("GEM_well_tibble", inherits = TRUE)) {
    GEM_well_tibble_obj <- get("GEM_well_tibble", inherits = TRUE)
    if ("GEM_well_ID" %in% names(GEM_well_tibble_obj)) {
      suffixes <- rbind(
        suffixes,
        data.frame(
          suffix = GEM_well_tibble_obj$GEM_well_ID,
          placeholder = "<GEM_well_ID>"
        )
      )
    }
  }
  if (exists("dataset_tibble", inherits = TRUE)) {
    dataset_tibble_obj <- get("dataset_tibble", inherits = TRUE)
    if ("dataset" %in% names(dataset_tibble_obj)) {
      suffixes <- rbind(
        suffixes,
        data.frame(
          suffix = dataset_tibble_obj$dataset,
          placeholder = "<dataset_name>"
        )
      )
    }
  }
  if (exists("aggregation_tibble", inherits = TRUE)) {
    aggregation_tibble_obj <- get("aggregation_tibble", inherits = TRUE)
    if ("aggregation" %in% names(aggregation_tibble_obj)) {
      suffixes <- rbind(
        suffixes,
        data.frame(
          suffix = aggregation_tibble_obj$aggregation,
          placeholder = "<aggregation_name>"
        )
      )
    }
  }
  if (exists("known_aggregation_modules", inherits = TRUE)) {
    suffixes <- rbind(
      suffixes,
      data.frame(
        suffix = get("known_aggregation_modules", inherits = TRUE),
        placeholder = "<module_name>"
      )
    )
  }
  suffixes$suffix <- as.character(suffixes$suffix)
  suffixes$placeholder <- as.character(suffixes$placeholder)
  suffixes <- suffixes[!is.na(suffixes$suffix) & nzchar(suffixes$suffix), , drop = FALSE]
  suffixes[!duplicated(suffixes$suffix), , drop = FALSE]
}

targets_graph_regex_escape <- function(x) {
  gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", x, perl = TRUE)
}

targets_graph_normalize_label_suffixes <- function(suffixes) {
  if (is.data.frame(suffixes)) {
    suffixes <- suffixes[, c("suffix", "placeholder"), drop = FALSE]
  } else {
    suffixes <- data.frame(suffix = suffixes, placeholder = "", stringsAsFactors = FALSE)
  }
  suffixes$suffix <- as.character(suffixes$suffix)
  suffixes$placeholder <- as.character(suffixes$placeholder)
  suffixes <- suffixes[!is.na(suffixes$suffix) & nzchar(suffixes$suffix), , drop = FALSE]
  suffixes <- suffixes[!duplicated(suffixes$suffix), , drop = FALSE]
  suffixes[order(nchar(suffixes$suffix), decreasing = TRUE), , drop = FALSE]
}

targets_graph_clean_node_labels <- function(glimpse_graph, suffixes = targets_graph_default_label_suffixes()) {
  suffixes <- targets_graph_normalize_label_suffixes(suffixes)
  if (!nrow(suffixes)) {
    return(glimpse_graph)
  }

  labels <- glimpse_graph$x$nodes$label
  for (i in seq_len(nrow(suffixes))) {
    labels <- gsub(
      pattern = paste0("\\.", targets_graph_regex_escape(suffixes$suffix[[i]]), "(?=\\.|$)"),
      replacement = paste0(".", suffixes$placeholder[[i]]),
      x = labels,
      perl = TRUE
    )
  }
  glimpse_graph$x$nodes$label <- labels
  glimpse_graph
}

targets_graph_collapse_duplicate_labels <- function(glimpse_graph) {
  nodes <- glimpse_graph$x$nodes
  edges <- glimpse_graph$x$edges
  duplicate_label <- duplicated(nodes$label)
  if (!any(duplicate_label)) {
    return(glimpse_graph)
  }

  representative <- nodes$name[match(nodes$label, nodes$label)]
  names(representative) <- nodes$name
  glimpse_graph$x$nodes <- nodes[!duplicate_label, , drop = FALSE]

  edges$from <- unname(representative[edges$from])
  edges$to <- unname(representative[edges$to])
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  glimpse_graph$x$edges <- edges[!duplicated(edges[, c("from", "to"), drop = FALSE]), , drop = FALSE]
  attr(glimpse_graph, "targets_graph_label_collapse_summary") <- list(
    kept_nodes = nrow(glimpse_graph$x$nodes),
    collapsed_nodes = sum(duplicate_label),
    edges = nrow(glimpse_graph$x$edges)
  )
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
    targets_graph_clean_node_labels(suffixes = label_suffixes) |>
    targets_graph_collapse_duplicate_labels()

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(targets_graph_mermaid_lines(graph), output_file)
  invisible(graph)
}

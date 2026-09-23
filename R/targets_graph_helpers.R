targets_graph_mermaid_lines <- function(glimpse_graph) {
  nodes <- glimpse_graph$x$nodes
  if (nrow(nodes) < 1L) {
    return("")
  }

  shape_open <- c(object = "{{", "function" = ">", stem = "([", pattern = "[")
  shape_close <- c(object = "}}", "function" = "]", stem = "])", pattern = "]")
  label <- gsub(">", "&gt;", gsub("<", "&lt;", nodes$label, fixed = TRUE), fixed = TRUE)
  text <- sprintf(
    "%s%s\"%s\"%s:::%s",
    nodes$name, shape_open[nodes$type], label, shape_close[nodes$type], nodes$status
  )
  names(text) <- nodes$name

  edges <- glimpse_graph$x$edges
  disconnected <- unname(text[setdiff(nodes$name, c(edges$from, edges$to))])
  c(
    "%%{init: {'flowchart': {'defaultRenderer': 'elk'}}}%%",
    "flowchart TB",
    sprintf("  %s --> %s", unname(text[edges$from]), unname(text[edges$to])),
    if (length(disconnected)) paste0("  ", disconnected)
  )
}

targets_graph_part_of_graph_tag <- function(graph_id) {
  paste0("[part_of_graph:", graph_id, "]")
}

targets_graph_part_of_graph_ids <- function(manifest = targets::tar_manifest(callr_function = NULL)) {
  tags <- unlist(regmatches(
    manifest$description,
    gregexpr("\\[part_of_graph:[A-Za-z0-9_]+\\]", manifest$description)
  ))
  sort(unique(sub("^\\[part_of_graph:(.*)\\]$", "\\1", tags)))
}

# Uses the mapping values defined by _targets.R in the calling session.
targets_graph_default_label_suffixes <- function() {
  manifest <- read_config_parameter_manifest("cfg_pipeline_parameters.tsv", scope = "aggregation")
  modules <- manifest$allowed_values[[match("modules", manifest$param_name)]]
  data.frame(
    suffix = c(GEM_well_tibble$GEM_well_ID, aggregation_tibble$aggregation, modules),
    placeholder = rep(
      c("<GEM_well_ID>", "<aggregation_name>", "<module_name>"),
      c(nrow(GEM_well_tibble), nrow(aggregation_tibble), length(modules))
    )
  )
}

targets_graph_clean_node_labels <- function(glimpse_graph, suffixes) {
  suffixes <- suffixes[!is.na(suffixes$suffix) & nzchar(suffixes$suffix), , drop = FALSE]
  suffixes <- suffixes[!duplicated(suffixes$suffix), , drop = FALSE]
  suffixes <- suffixes[order(nchar(suffixes$suffix), decreasing = TRUE), , drop = FALSE]
  labels <- glimpse_graph$x$nodes$label
  for (i in seq_len(nrow(suffixes))) {
    escaped_suffix <- gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", suffixes$suffix[[i]], perl = TRUE)
    labels <- gsub(
      paste0("\\.", escaped_suffix, "(?=\\.|$)"),
      paste0(".", suffixes$placeholder[[i]]),
      labels,
      perl = TRUE
    )
  }
  glimpse_graph$x$nodes$label <- labels
  glimpse_graph
}

targets_graph_collapse_duplicate_labels <- function(glimpse_graph) {
  nodes <- glimpse_graph$x$nodes
  representative <- stats::setNames(nodes$name[match(nodes$label, nodes$label)], nodes$name)
  edges <- glimpse_graph$x$edges
  edges$from <- unname(representative[edges$from])
  edges$to <- unname(representative[edges$to])
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  glimpse_graph$x$nodes <- nodes[!duplicated(nodes$label), , drop = FALSE]
  glimpse_graph$x$edges <- edges[!duplicated(edges[, c("from", "to")]), , drop = FALSE]
  glimpse_graph
}

# Breadth-first search through removed nodes to the nearest retained nodes.
targets_graph_nearest_retained <- function(start, adjacency, removed) {
  retained <- character()
  seen <- start
  queue <- unname(adjacency[[start]])
  while (length(queue)) {
    current <- queue[[1]]
    queue <- queue[-1]
    if (current %in% seen) {
      next
    }
    seen <- c(seen, current)
    if (current %in% removed) {
      queue <- c(queue, unname(adjacency[[current]]))
    } else {
      retained <- c(retained, current)
    }
  }
  unique(retained)
}

targets_graph_bypass_removed_edges <- function(edges, removed) {
  kept_edges <- edges[!edges$from %in% removed & !edges$to %in% removed, , drop = FALSE]
  forward <- split(as.character(edges$to), as.character(edges$from))
  reverse <- split(as.character(edges$from), as.character(edges$to))
  bypass_pairs <- unique(do.call(rbind, c(
    list(data.frame(from = character(), to = character())),
    lapply(removed, function(node) {
      expand.grid(
        from = targets_graph_nearest_retained(node, reverse, removed),
        to = targets_graph_nearest_retained(node, forward, removed),
        stringsAsFactors = FALSE
      )
    })
  )))
  bypass_pairs <- bypass_pairs[
    bypass_pairs$from != bypass_pairs$to &
      !paste(bypass_pairs$from, bypass_pairs$to, sep = "\r") %in% paste(kept_edges$from, kept_edges$to, sep = "\r"),
    ,
    drop = FALSE
  ]
  if (!nrow(bypass_pairs)) {
    return(kept_edges)
  }
  for (column in setdiff(names(edges), c("from", "to"))) {
    bypass_pairs[[column]] <- edges[[column]][[1]]
  }
  unique(rbind(kept_edges, bypass_pairs[, names(edges), drop = FALSE]))
}

targets_graph_prune_to_part_of_graph <- function(glimpse_graph, graph_id) {
  tag <- targets_graph_part_of_graph_tag(graph_id)
  nodes <- glimpse_graph$x$nodes
  keep_node <- grepl(tag, nodes$description, fixed = TRUE)
  if (!any(keep_node)) {
    stop("No nodes in the graph have tag ", tag, call. = FALSE)
  }
  glimpse_graph$x$nodes <- nodes[keep_node, , drop = FALSE]
  glimpse_graph$x$edges <- targets_graph_bypass_removed_edges(glimpse_graph$x$edges, nodes$name[!keep_node])
  glimpse_graph
}

targets_graph_write_part_of_graph_mermaid <- function(
  graph_id,
  output_file,
  manifest = targets::tar_manifest(callr_function = NULL),
  label_suffixes = targets_graph_default_label_suffixes()
) {
  target_names <- manifest$name[grepl(targets_graph_part_of_graph_tag(graph_id), manifest$description, fixed = TRUE)]
  graph <- rlang::inject(targets::tar_glimpse(names = tidyselect::any_of(!!target_names))) |>
    targets_graph_prune_to_part_of_graph(graph_id = graph_id) |>
    targets_graph_clean_node_labels(suffixes = label_suffixes) |>
    targets_graph_collapse_duplicate_labels()

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(targets_graph_mermaid_lines(graph), output_file)
  invisible(graph)
}

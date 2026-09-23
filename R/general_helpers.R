#' Collapse duplicate names
#'
#' Group equally named entries, concatenate their values one level, and
#' recursively collapse nested named lists. Groups follow `split()` name order.
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


label_plot_variable <- function(variable) {
  labels <- c(nCount_RNA = "RNA UMI count", nFeature_RNA = "Detected RNA genes",
    RNA_mito_percent = "Mitochondrial RNA (%)", GEM_well_ID = "GEM well", donor_id = "Donor",
    log10_nCount_RNA = "Log10 RNA UMI count", nCount_ATAC = "ATAC peak count",
    log10_nCount_ATAC = "Log10 ATAC peak count", ATAC.weight = "ATAC WNN weight",
    PCA_harmony_SNN_cluster = "GEX cluster", PCA_harmony_SNN_cluster_named = "Named GEX cluster",
    PCA_harmony_SNN_cluster_cell_type = "GEX-assigned cell type",
    LSI_harmony_SNN_cluster = "ATAC cluster", LSI_harmony_SNN_cluster_named = "Named ATAC cluster",
    LSI_harmony_SNN_cluster_cell_type = "ATAC-cluster cell type (GEX-derived)",
    WNN_harmony_SNN_cluster = "WNN cluster", WNN_harmony_SNN_cluster_named = "Named WNN cluster",
    WNN_harmony_SNN_cluster_cell_type = "WNN-cluster cell type (GEX-derived)")
  dplyr::coalesce(unname(labels[variable]), gsub("[_.]", " ", variable))
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
#' @param legend_continuous Continuous legend mode passed to `BPCells::plot_embedding()`,
#'   for example `quantile`.
#' @param quantile_range Numeric colour-clipping quantiles.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_UMAP_from_metadata <- function(
  metadata_tibble,
  variable,
  value_source = c("metadata", "feature"),
  feature_matrix = NULL,
  umap_cols = c("LSI_UMAP_1", "LSI_UMAP_2"),
  legend_continuous = "quantile",
  quantile_range = c(0.01, 0.99)
) {
  value_source <- match.arg(value_source)
  values <- if (value_source == "metadata") {
    metadata_tibble[[variable]]
  } else {
    as.numeric(as.matrix(feature_matrix[variable, metadata_tibble$barcode_w_prefix, drop = FALSE]))
  }
  plot_col_data <- stats::setNames(data.frame(values), variable)
  plot_embedding_matrix <- as.matrix(metadata_tibble[umap_cols])
  if (is.numeric(values)) {
    keep_rows <- is.finite(values)
    if (!any(keep_rows)) {
      return(structure(list(), class = c("empty_plot_list", "list")))
    }
    plot_col_data <- plot_col_data[keep_rows, , drop = FALSE]
    plot_embedding_matrix <- plot_embedding_matrix[keep_rows, , drop = FALSE]
  }

  # Point size and alpha shrink log-linearly from 5,000 to 500,000 cells.
  log_n_fraction <- (min(max(log10(nrow(metadata_tibble)), log10(5000)), log10(5e5)) - log10(5000)) /
    (log10(5e5) - log10(5000))
  plot <- BPCells::plot_embedding(
    source = plot_col_data,
    embedding = plot_embedding_matrix,
    features = variable,
    size = 2 - (2 - 1) * log_n_fraction,
    rasterize = TRUE,
    raster_pixels = 1024,
    randomize_order = TRUE,
    quantile_range = quantile_range,
    labels_discrete = FALSE,
    legend_continuous = legend_continuous,
    return_plot_list = TRUE,
    apply_styling = TRUE
  )
  if (inherits(plot, "ggplot") && length(plot$layers) > 0) {
    plot$layers[[1]]$aes_params$alpha <- 1 - (1 - 0.5) * log_n_fraction
  }
  plot + ggplot2::labs(
    title = paste(if (value_source == "feature") variable else label_plot_variable(variable),
      "on", gsub("LSI", "ATAC", gsub("_", " ", sub("_UMAP.*", "", umap_cols[[1]]))), "UMAP"),
    subtitle = stringr::str_wrap("Look for coherent local patterns or sample-specific separation; distances between islands do not measure biological difference.", width = 100),
    caption = stringr::str_wrap(paste(
      "Each point is a cell from the supplied metadata; only finite values are shown for numeric features.",
      if (is.numeric(values)) paste0(
        "Colour limits use quantiles ", paste(quantile_range, collapse = "-"), "; values outside are clipped for display.",
        if (legend_continuous == "quantile") " Legend endpoints identify quantiles, not absolute values.") else
        "Colours identify categories; categorical labels do not establish independent biological validation.",
      if (value_source == "feature") "Feature values come directly from the supplied matrix without normalization in this plotting helper."), width = 110))
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
  modalities <- c(GEX = "PCA_harmony_SNN", WNN = "WNN_harmony_SNN", ATAC = "LSI_harmony_SNN")
  embeddings <- list(GEX = c("GEX_UMAP_1", "GEX_UMAP_2"),
    WNN = c("WNN_UMAP_1", "WNN_UMAP_2"), ATAC = c("LSI_UMAP_1", "LSI_UMAP_2"))
  panels <- purrr::imap(modalities, function(cluster_prefix, modality) {
    cluster_col <- paste0(cluster_prefix, "_cluster_", cluster_col_suffix)
    row_title <- paste0(modality, "-derived ",
      if (cluster_col_suffix == "named") "clusters" else "cell-type labels")
    strip <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = row_title, hjust = 0,
        size = 3.5, fontface = "bold", colour = "grey20") +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(.01, .01))) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "grey92", colour = "grey75"),
        plot.margin = ggplot2::margin(0, 4, 0, 4))
    plots <- purrr::imap(embeddings, function(umap_cols, embedding_name) {
      centers <- metadata_tibble |>
        dplyr::transmute(label = as.character(.data[[cluster_col]]),
          x = .data[[umap_cols[[1]]]], y = .data[[umap_cols[[2]]]]) |>
        dplyr::filter(!is.na(label), is.finite(x), is.finite(y)) |>
        dplyr::summarise(x = median(x), y = median(y), .by = label)
      plot_UMAP_from_metadata(metadata_tibble, variable = cluster_col, umap_cols = umap_cols) +
        ggrepel::geom_text_repel(data = centers, ggplot2::aes(x, y, label = label),
          inherit.aes = FALSE, size = 2, colour = "grey20", seed = 1,
          box.padding = 0.15, point.padding = 0.05, segment.size = 0.2,
          segment.alpha = 0.5, max.overlaps = Inf) +
        ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(), legend.position = "none", axis.line = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour = "grey85", fill = NA),
          plot.title = ggplot2::element_text(size = 11, hjust = 0.5)) +
        ggplot2::labs(title = paste(embedding_name, "UMAP"), subtitle = NULL, caption = NULL,
          x = NULL, y = NULL)
    })
    c(list(strip), unname(plots))
  }) |> purrr::flatten()
  patchwork::wrap_plots(panels, design = "AAA\nBCD\nEEE\nFGH\nIII\nJKL",
    heights = do.call(grid::unit.c, rep(list(grid::unit(7, "mm"), grid::unit(1, "null")), 3))) +
    patchwork::plot_annotation(
      title = paste(if (cluster_col_suffix == "named") "Cluster labels" else "Cell-type labels",
        "across GEX, ATAC and WNN embeddings"),
      subtitle = "Row strips identify the source of the labels; columns identify the embedding. Compare each row for splits or mixing.",
      caption = stringr::str_wrap(paste(
        "All nine panels use the same supplied cells. Small labels are repelled from each group's median coordinates; colours identify groups within each row.",
        "GEX-derived cell-type names are not independent evidence of agreement. UMAP island distances and orientations are not directly comparable across embeddings."
      ), 150))
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

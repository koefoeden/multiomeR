plot_nuclei_per_donor_id <- function(
  metadata_tibble,
  fill_by = "PCA_harmony_SNN_cluster"
) {
  metadata_tibble %>%
    ggplot2::ggplot(ggplot2::aes(x = donor_id, fill = .data[[fill_by]])) +
    ggplot2::geom_bar() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust = 1))
}

#' Plot markers volcano simple
#'
#' Draw faceted marker volcano plots with top up/down genes labelled per cluster.
#'
#' @param markers_tibble Marker-results tibble with at least `cluster`, `gene`,
#'   `avg_log2FC`, `p_val`, and `p_val_adj` columns.
#' @return A faceted ggplot with log-fold change on x, `-log10(p_val)` on y,
#'   adjusted-significance color, and up to 20 labels per direction and cluster.
#' @keywords internal

plot_markers_volcano_simple <- function(markers_tibble) {
  markers_tibble_formatted <- markers_tibble %>%
    dplyr::mutate(
      p_val = dplyr::case_when(p_val == 0 ~ .Machine$double.xmin, .default = p_val),
      p_val_adj = dplyr::case_when(
        p_val_adj == 0 ~ .Machine$double.xmin,
        .default = p_val_adj
      )
    ) %>%
    dplyr::arrange(p_val, dplyr::desc(abs(avg_log2FC)))

  # 20 top labels in each direction per cluster
  top_labels <- markers_tibble_formatted %>%
    dplyr::mutate(direction = dplyr::case_when(avg_log2FC > 0 ~ "up", TRUE ~ "down")) %>%
    dplyr::group_by(cluster, direction) %>%
    dplyr::slice_head(n = 20)

  plot <- markers_tibble_formatted %>%
    ggplot2::ggplot(ggplot2::aes(
      x = avg_log2FC,
      y = -log10(p_val),
      color = p_val_adj < 0.05,
      label = gene
    )) +
    ggrastr::geom_point_rast(alpha = 0.5, size = 0.5) +
    ggplot2::theme(legend.position = "none") +
    ggrepel::geom_text_repel(
      data = top_labels,
      ggplot2::aes(label = gene),
      max.overlaps = 30,
      size = 2
    ) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::geom_vline(xintercept = 0, lty = 2) +
    ggplot2::facet_wrap(~cluster, scales = "free") +
    ggplot2::scale_x_continuous(limits = symmetric_limits)

  return(plot)
}

#' Plot categorical bars plot
#'
#' Plot categorical metadata composition within each cluster.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param metadata_cols Character vector of metadata columns to test or plot.
#' @param optional_metadata_cols Metadata columns to include only when present, so shared plotting code can span datasets with different annotations.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_categorical_bars_plot <- function(
  metadata_tibble,
  metadata_cols,
  optional_metadata_cols = NULL,
  cluster_col
) {
  if (!cluster_col %in% colnames(metadata_tibble)) {
    stop("Cluster column not found in metadata: ", cluster_col)
  }

  metadata_cols <- metadata_cols %||% character()
  missing_metadata_cols <- setdiff(metadata_cols, colnames(metadata_tibble))
  if (length(missing_metadata_cols) > 0) {
    stop("Required metadata column(s) not found: ", paste(missing_metadata_cols, collapse = ", "))
  }
  optional_metadata_cols <- intersect(optional_metadata_cols %||% character(), colnames(metadata_tibble))
  plot_metadata_cols <- unique(c(metadata_cols, optional_metadata_cols))

  if (length(plot_metadata_cols) == 0) {
    return(list())
  }

  metadata <- metadata_tibble %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(plot_metadata_cols), as.character)) %>%
    dplyr::select(dplyr::all_of(c(plot_metadata_cols, cluster_col)))

  plot_metadata_cols %>%
    purrr::set_names() %>%
    purrr::map(
      ~ metadata %>%
        dplyr::summarise(n_nuclei = dplyr::n(), .by = dplyr::all_of(c(.x, cluster_col))) %>%
        ggplot2::ggplot(ggplot2::aes(y = .data[[cluster_col]], x = n_nuclei, fill = .data[[.x]])) +
        ggplot2::geom_col() +
        ggplot2::theme(legend.position = "bottom")
    )
}

#' Plot cluster confusion matrix
#'
#' Plot row-normalized overlap between two cluster or annotation columns.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param source_col Metadata column defining heatmap rows.
#' @param target_col Metadata column defining heatmap columns.
#' @param source_label Axis label for source rows.
#' @param target_label Axis label for target columns.
#' @param title Optional plot title.
#' @return A ggplot tile heatmap where fill is each target count divided by the
#'   source-row total and tile labels show raw counts.
#' @keywords internal

plot_cluster_confusion_matrix <- function(
  metadata_tibble,
  source_col,
  target_col,
  source_label,
  target_label,
  title = NULL
) {
  missing_cols <- setdiff(c(source_col, target_col), names(metadata_tibble))
  if (length(missing_cols) > 0) {
    stop("Missing confusion-matrix metadata column(s): ", paste(missing_cols, collapse = ", "))
  }

  counts_tibble <- metadata_tibble |>
    dplyr::transmute(
      source = as.character(.data[[source_col]]),
      target = as.character(.data[[target_col]])
    ) |>
    dplyr::filter(!is.na(source), !is.na(target), source != "", target != "") |>
    dplyr::count(source, target, name = "n")

  source_levels <- gtools::mixedsort(unique(counts_tibble$source))
  target_levels <- gtools::mixedsort(unique(counts_tibble$target))

  counts_tibble |>
    tidyr::complete(source = source_levels, target = target_levels, fill = list(n = 0L)) |>
    dplyr::mutate(source_total = sum(n), .by = source) |>
    dplyr::mutate(
      source = factor(source, levels = rev(source_levels)),
      target = factor(target, levels = target_levels),
      source_fraction = dplyr::if_else(source_total > 0, n / source_total, 0),
      n_label = dplyr::case_when(
        n == 0 ~ "",
        n >= 1000 ~ paste0(scales::number(n / 1000, accuracy = 0.1), "k"),
        .default = scales::number(n, accuracy = 1)
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = target, y = source, fill = source_fraction)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = n_label), size = 2.5) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#2166AC",
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(title = title, x = target_label, y = source_label, fill = "row fraction") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

#' Plot ATAC vs RNA weight boxplots
#'
#' Compare WNN ATAC weights and RNA/ATAC depth metrics across clusters.
#'
#' @param metadata Cell metadata containing WNN `ATAC.weight`, RNA/ATAC depth
#'   columns, and the cluster label column.
#' @param cluster_label_col Single column name used for cluster label col; the column must exist in the relevant metadata tibble.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_ATAC_vs_RNA_weight_boxplots <- function(
  metadata,
  cluster_label_col = "WNN_harmony_SNN_cluster"
) {
  clusters_sorted_by_median_ATAC_weights <- metadata %>%
    tibble::as_tibble() %>%
    dplyr::group_by(.data[[cluster_label_col]]) %>%
    dplyr::summarise(median_ATAC_weight = stats::median(ATAC.weight)) %>%
    dplyr::arrange(dplyr::desc(median_ATAC_weight)) %>%
    dplyr::pull(1) %>%
    as.vector()

  boxplot <- metadata %>%
    tibble::as_tibble() %>%
    dplyr::select(dplyr::all_of(c(
      cluster_label_col,
      "ATAC.weight",
      "log10_nCount_RNA",
      "log10_nCount_ATAC"
    ))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("ATAC.weight", "log10_nCount_RNA", "log10_nCount_ATAC"))
    ) %>%
    dplyr::mutate(
      name = factor(
        name,
        levels = c("log10_nCount_RNA", "ATAC.weight", "log10_nCount_ATAC")
      ),
      sorted_cluster = factor(
        .data[[cluster_label_col]],
        levels = clusters_sorted_by_median_ATAC_weights
      )
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = value, y = sorted_cluster)) +
    ggplot2::geom_boxplot(ggplot2::aes(color = .data[[cluster_label_col]])) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_wrap(~name, scales = "free_x") +
    ggplot2::labs(x = "log10(counts) / WNN weight / log10(counts)", y = "Cluster")

  return(boxplot)
}

#' Plot marker expression dot BPCells
#'
#' Plot BPCells marker-expression dot plots from GEX metadata and marker sets.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param features Character vector of feature names to extract from the matrix row names; missing features are handled by the called helper.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @param scale_factor Scale factor used when normalizing counts, coverage, or marker scores.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_marker_expression_dot_BPCells <- function(feature_matrix, metadata_tibble, features, group_col, scale_factor = 10000) {
  features <- intersect(features, rownames(feature_matrix))
  if (length(features) == 0) {
    stop("No requested marker features were found in the feature matrix.")
  }

  metadata <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", group_col)), dplyr::any_of("nCount_RNA")) |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% colnames(feature_matrix), !is.na(.data[[group_col]])) |>
    dplyr::arrange(match(.data$barcode_w_prefix, colnames(feature_matrix)))

  counts_matrix <- feature_matrix[, metadata$barcode_w_prefix, drop = FALSE]
  cell_counts <- if ("nCount_RNA" %in% colnames(metadata)) {
    metadata$nCount_RNA
  } else {
    BPCells::colSums(counts_matrix)
  }
  log_norm_matrix <- counts_matrix |>
    (\(matrix_in) if (inherits(matrix_in, "matrix")) Matrix::Matrix(matrix_in, sparse = TRUE) else matrix_in)() |>
    BPCells::multiply_cols(ifelse(cell_counts > 0, scale_factor / cell_counts, 0)) |>
    BPCells::log1p_slow()

  BPCells::plot_dot(
    source = log_norm_matrix,
    features = features,
    groups = metadata[[group_col]],
    group_order = levels(as.factor(metadata[[group_col]])),
    gene_mapping = NULL
  ) +
    ggplot2::labs(y = group_col)
}

#' Plot marker gene activity dot BPCells
#'
#' Plot gene-activity marker dot plots for ATAC-derived activity matrices.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param features Character vector of feature names to extract from the matrix row names; missing features are handled by the called helper.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_marker_gene_activity_dot_BPCells <- function(feature_matrix, metadata_tibble, features, group_col) {
  features <- intersect(features, rownames(feature_matrix))
  if (length(features) == 0) {
    stop("No requested marker features were found in the feature matrix.")
  }

  metadata <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", group_col))) |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% colnames(feature_matrix), !is.na(.data[[group_col]])) |>
    dplyr::arrange(match(.data$barcode_w_prefix, colnames(feature_matrix)))

  BPCells::plot_dot(
    source = BPCells::log1p_slow(feature_matrix[, metadata$barcode_w_prefix, drop = FALSE]),
    features = features,
    groups = metadata[[group_col]],
    group_order = levels(as.factor(metadata[[group_col]])),
    gene_mapping = NULL
  ) +
    ggplot2::labs(y = group_col)
}

#' Plot module scores dot for metadata
#'
#' Validate YAML/TSV configuration tables before constructing the active targets graph.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param marker_genes_list Named list of marker genes or modules; names become labels in dot plots and metadata-score summaries.
#' @param cluster_by Single column name used for cluster by; the column must exist in the relevant metadata tibble.
#' @param positive_threshold Value above which a feature/module score is counted as positive in dot-plot summaries.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_module_scores_dot_for_metadata <- function(metadata_tibble, marker_genes_list, cluster_by, positive_threshold = 0) {
  module_names <- intersect(names(marker_genes_list), colnames(metadata_tibble))
  if (length(module_names) == 0) {
    stop("No marker module score columns were found in metadata.")
  }

  plot_tibble <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c(cluster_by, module_names))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(module_names), names_to = "module", values_to = "score") |>
    dplyr::filter(!is.na(.data[[cluster_by]]), is.finite(.data$score)) |>
    dplyr::summarise(
      mean_score = mean(.data$score),
      pct_positive = 100 * mean(.data$score > positive_threshold),
      .by = c(dplyr::all_of(cluster_by), "module")
    ) |>
    dplyr::mutate(
      module = factor(.data$module, levels = module_names),
      cluster = factor(.data[[cluster_by]], levels = levels(as.factor(metadata_tibble[[cluster_by]])))
    )

  plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$module, y = .data$cluster, color = .data$mean_score, size = .data$pct_positive)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_area(limits = c(0, 100), max_size = 7) +
    ggplot2::scale_color_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0) +
    ggplot2::labs(x = "Marker module", y = cluster_by, color = "Mean score", size = "% score > 0") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot feature scores heatmap from matrix
#'
#' Summarize feature-matrix scores by metadata group in a heatmap.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param features Character vector of feature names to extract from the matrix row names; missing features are handled by the called helper.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_feature_scores_heatmap_from_matrix <- function(feature_matrix, metadata_tibble, features, group_col) {
  requested_features <- features
  if (length(requested_features) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void()
    )
  }
  features <- intersect(requested_features, rownames(feature_matrix))
  if (length(features) == 0) {
    stop(
      "None of the requested score features were found in the feature matrix: ",
      paste(requested_features, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", group_col))) |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% colnames(feature_matrix), !is.na(.data[[group_col]])) |>
    dplyr::arrange(match(.data$barcode_w_prefix, colnames(feature_matrix)))

  score_matrix <- as.matrix(feature_matrix[features, metadata$barcode_w_prefix, drop = FALSE])
  plot_tibble <- t(score_matrix) |>
    tibble::as_tibble(.name_repair = "minimal") |>
    dplyr::mutate(group = metadata[[group_col]]) |>
    tidyr::pivot_longer(cols = dplyr::all_of(features), names_to = "feature", values_to = "score") |>
    dplyr::summarise(
      mean_score = mean(.data$score),
      .by = c("group", "feature")
    ) |>
    dplyr::mutate(
      feature = factor(.data$feature, levels = features),
      group = factor(.data$group, levels = levels(as.factor(metadata[[group_col]])))
    )

  plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$feature, y = .data$group, fill = .data$mean_score)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0) +
    ggplot2::labs(x = "Feature", y = group_col, fill = "Mean score") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Plot similarity matrix from GRanges list
#'
#' Plot pairwise overlap similarity between named GRanges collections.
#'
#' @param GRanges_list_in GRanges object containing GRanges list in coordinates and metadata.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_similarity_matrix_from_GRanges_list <- function(GRanges_list_in) {
  GRanges_list <- GRanges_list_in |>
    as.list() |>
    purrr::map(\(gr) GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE))
  set_names <- as.character(seq_along(GRanges_list))

  all_peaks <- GenomicRanges::GRangesList(GRanges_list) |>
    unlist(use.names = FALSE)
  cluster_idx <- rep(seq_along(GRanges_list), lengths(GRanges_list))

  bins <- GenomicRanges::disjoin(all_peaks, ignore.strand = TRUE)
  hits <- GenomicRanges::findOverlaps(bins, all_peaks, ignore.strand = TRUE)
  query_hits <- S4Vectors::queryHits(hits)
  subject_hits <- S4Vectors::subjectHits(hits)

  membership_matrix <- Matrix::sparseMatrix(
    i = query_hits,
    j = cluster_idx[subject_hits],
    x = rep(1.0, length(query_hits)),
    dims = c(length(bins), length(GRanges_list))
  )

  weighted_membership_matrix <- membership_matrix
  weighted_membership_matrix@x <- as.numeric(GenomicRanges::width(bins)[weighted_membership_matrix@i + 1L])

  intersection_matrix <- Matrix::crossprod(membership_matrix, weighted_membership_matrix) |>
    as.matrix()
  set_widths <- Matrix::colSums(weighted_membership_matrix) |>
    as.numeric()

  union_matrix <- outer(set_widths, set_widths, "+") - intersection_matrix
  min_width_matrix <- outer(set_widths, set_widths, pmin)

  jaccard_matrix <- intersection_matrix / union_matrix
  fraction_matrix <- intersection_matrix / min_width_matrix

  result_matrix <- jaccard_matrix
  result_matrix[lower.tri(result_matrix)] <- fraction_matrix[lower.tri(fraction_matrix)]
  dimnames(result_matrix) <- list(set_names, set_names)

  result_tibble <- result_matrix %>%
    tibble::as_tibble() %>%
    tibble::rownames_to_column(var = "x") %>%
    tidyr::pivot_longer(cols = 2:dplyr::last_col(), names_to = "y") %>%
    dplyr::mutate(
      x = as.numeric(x),
      y = as.numeric(y)
    )

  tile_plot <- result_tibble %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::labs(x = "Cluster", y = "Cluster") +
    ggplot2::ggtitle("Jaccard index (upper) vs overlap fraction (lower)") +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 3) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, length(set_names), 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(1, length(set_names), 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1))

  return(tile_plot)
}

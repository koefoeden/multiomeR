#' Plot GWAS chromVAR peak weights summary
#'
#' Plot how many peaks receive nonzero GWAS chromVAR weights per trait.
#'
#' @param peak_weight_records List of GWAS peak-weight records, each containing
#'   `GWAS_ID` and a named `peak_weights_vec`.
#' @param overlap_threshold Number of nonzero peak overlaps used as the dashed
#'   threshold in the overlap-fraction facet.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_chromVAR_peak_weights_summary <- function(peak_weight_records, overlap_threshold = 100) {
  total_peaks <- length(peak_weight_records[[1]]$peak_weights_vec)
  overlap_threshold_fraction <- overlap_threshold / total_peaks

  plot_tibble <- peak_weight_records |>
    get_GWAS_chromVAR_peak_weight_summary_tibble() %>%
    tidyr::pivot_longer(!dplyr::matches("GWAS_ID|_sheet"))

  vline_tibble <- tibble::tibble(
    name = "frac_overlapped_peaks",
    xintercept = overlap_threshold_fraction
  )

  plot_tibble %>%
    ggplot2::ggplot(ggplot2::aes(y = GWAS_ID, x = value)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      data = vline_tibble,
      ggplot2::aes(xintercept = xintercept),
      inherit.aes = FALSE,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~name, scales = "free_x") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = stringr::str_glue("Peak x GWAS posterior-probability overlap summary, from {total_peaks} total peaks"),
      subtitle = stringr::str_wrap("Look for traits with limited peak overlap before interpreting enrichment; low representation can limit downstream sensitivity.", width = 100),
      caption = stringr::str_glue(
        "Dashed line in the frac_overlapped_peaks facet marks {overlap_threshold} overlapped peaks ({scales::percent(overlap_threshold_fraction)} of all peaks)."
      )
    )
}

#' Plot GWAS credible-set similarity
#'
#' @param similarity_matrix Output of `get_GWAS_credible_set_similarity_matrix()`.
#' @param GWAS_inputs_tibble GWAS metadata with `GWAS_ID` and `Category`.
#' @return A ggplot of the lower triangle, ordered by configured category.
#' @keywords internal

plot_GWAS_credible_set_similarity <- function(similarity_matrix, GWAS_inputs_tibble) {
  ordered_metadata <- dplyr::arrange(GWAS_inputs_tibble, .data$Category, .data$GWAS_ID)
  GWAS_IDs <- ordered_metadata$GWAS_ID
  GWAS_labels <- stringr::str_replace_all(GWAS_IDs, "_", " ")
  category_breaks <- get_plot_group_breaks(ordered_metadata$Category)
  n_GWAS <- length(GWAS_IDs)
  # Category separators stop at the diagonal of the lower triangle.
  separator_tibble <- tibble::tibble(
    x = c(rep(0.5, length(category_breaks)), category_breaks),
    xend = c(category_breaks, category_breaks),
    y = c(category_breaks, category_breaks),
    yend = c(category_breaks, rep(n_GWAS + 0.5, length(category_breaks)))
  )
  plot_tibble <- tidyr::expand_grid(row = seq_len(n_GWAS), column = seq_len(n_GWAS)) |>
    dplyr::filter(.data$row > .data$column) |>
    dplyr::mutate(similarity = similarity_matrix[cbind(GWAS_IDs[.data$row], GWAS_IDs[.data$column])])

  ggplot2::ggplot(plot_tibble, ggplot2::aes(x = column, y = row, fill = similarity)) +
    ggplot2::geom_tile(color = "grey85", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = dplyr::if_else(similarity >= 0.005, sprintf("%.2f", similarity), "")),
      size = 2
    ) +
    ggplot2::geom_segment(
      data = separator_tibble,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      color = "grey25",
      linewidth = 0.4
    ) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_GWAS), labels = GWAS_labels, expand = c(0, 0)) +
    ggplot2::scale_y_reverse(breaks = seq_len(n_GWAS), labels = GWAS_labels, expand = c(0, 0)) +
    ggplot2::scale_fill_gradient(
      name = "Shared PIP mass",
      low = "white",
      high = "#08519C",
      limits = c(0, NA),
      transform = "sqrt"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "GWAS credible-set similarity",
      subtitle = "Similar pairs share peak weights, so their enrichment results are not independent.",
      caption = paste(
        "Shared mass of the per-GWAS normalized posterior probabilities of identical variants (one minus the total variation distance);",
        "not a colocalization test. Lines separate configured GWAS categories; labels show similarities of at least 0.005."
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "right"
    )
}

#' Get SCAVENGE TRS UMAP plots
#'
#' Build UMAP overlay plots from cell-level SCAVENGE TRS scores.
#'
#' @param TRS_tibble Cell-level TRS tibble with `barcode_w_prefix`, `GWAS_ID`,
#'   and `score` columns.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param umap_cols Two metadata columns used as UMAP x/y coordinates.
#' @param label_col Optional metadata column used to label group centroids on
#'   the UMAP.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

get_SCAVENGE_TRS_UMAP_plots <- function(TRS_tibble, metadata_tibble, umap_cols, label_col = NULL) {
  if (nrow(TRS_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  umap_cols <- unlist(as.list(umap_cols), use.names = FALSE)
  GWAS_ID <- unique(TRS_tibble$GWAS_ID)
  score_col <- stringr::str_c("score_", GWAS_ID[[1]])
  SCAVENGE_metadata_tibble <- metadata_tibble |>
    dplyr::left_join(dplyr::select(TRS_tibble, barcode_w_prefix, score), by = "barcode_w_prefix") |>
    dplyr::rename(!!score_col := score)

  plot <- plot_UMAP_from_metadata(
    metadata_tibble = SCAVENGE_metadata_tibble,
    variable = score_col,
    umap_cols = umap_cols,
    legend_continuous = "value",
    quantile_range = c(0, 1)
  )

  if (inherits(plot, "empty_plot_list")) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot <- plot +
    ggplot2::scale_color_viridis_c(name = "TRS score") +
    ggplot2::labs(title = paste("SCAVENGE trait relevance:", GWAS_ID[[1]]),
      subtitle = stringr::str_wrap("Look for localized high trait-relevance scores; these are propagated scores, not probabilities of causal involvement.", width = 100),
      caption = stringr::str_wrap(paste(plot$labels$caption,
        "Colours show SCAVENGE TRS on the supplied graph embedding. Cell-type labels, when present, come from GEX-derived annotations."), width = 110))

  if (!is.null(label_col)) {
    label_tibble <- SCAVENGE_metadata_tibble |>
      dplyr::filter(
        !is.na(.data[[label_col]]),
        is.finite(.data[[umap_cols[[1]]]]),
        is.finite(.data[[umap_cols[[2]]]])
      ) |>
      dplyr::summarise(
        label = dplyr::first(as.character(.data[[label_col]])),
        UMAP_1 = stats::median(.data[[umap_cols[[1]]]]),
        UMAP_2 = stats::median(.data[[umap_cols[[2]]]]),
        .by = dplyr::all_of(label_col)
      ) |>
      dplyr::select(label, UMAP_1, UMAP_2)

    if (nrow(label_tibble) > 0) {
      plot <- plot +
        ggrepel::geom_label_repel(
          data = label_tibble,
          ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
          inherit.aes = FALSE,
          size = 2.4,
          label.size = 0.15,
          label.padding = grid::unit(0.08, "lines"),
          min.segment.length = 0,
          max.overlaps = Inf,
          seed = 1
        )
    }
  }

  plot
}

plot_SCAVENGE_summary_score_intervals <- function(summary_tibble) {
  if (nrow(summary_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot_tibble <- summary_tibble |>
    dplyr::mutate(
      facet_id = stringr::str_c(grouping_col, GWAS_ID, sep = "___"),
      group_reorder = tidytext::reorder_within(cluster, median_score, facet_id)
    )

  split(plot_tibble, plot_tibble$grouping_col) |>
    purrr::map(\(group_tibble) {
      ggplot2::ggplot(group_tibble, ggplot2::aes(x = group_reorder, fill = cluster)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = min_score, ymax = max_score), color = "grey35") +
        ggplot2::geom_crossbar(ggplot2::aes(y = median_score, ymin = q25_score, ymax = q75_score), width = 0.6) +
        ggplot2::facet_wrap(~GWAS_ID, scales = "free") +
        tidytext::scale_x_reordered(labels = function(x) gsub("___.*$", "", x)) +
        ggplot2::labs(title = "SCAVENGE score distributions by group",
          subtitle = "Compare median scores and spread; wide ranges can indicate a subset of high-scoring cells.",
          caption = stringr::str_wrap("Centre: group median; box: 25th-75th percentiles; line: full minimum-maximum range. These are cell-score distributions, not confidence intervals. Groups are ordered separately within each trait; facet scales may differ.", width = 110),
          x = NULL, y = "SCAVENGE TRS", fill = "Group") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
    })
}

add_GWAS_heatmap_categories <- function(heatmap_data, GWAS_tibble) {
  heatmap_data |>
    dplyr::inner_join(GWAS_tibble |> dplyr::select(GWAS_ID, Category), by = "GWAS_ID")
}

make_named_heatmap_palette <- function(values, palette = "Set3") {
  values <- sort(unique(stats::na.omit(values)))
  if (length(values) == 0) {
    return(character())
  }

  max_brewer_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  colors <- if (length(values) <= max_brewer_n) {
    RColorBrewer::brewer.pal(max(3, length(values)), palette)[seq_along(values)]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(max_brewer_n, palette))(length(values))
  }
  rlang::set_names(colors, values)
}

gwas_heatmap_metadata_theme <- function(show_y = FALSE, show_x = FALSE) {
  ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = if (show_x) ggplot2::element_text(angle = 45, hjust = 1) else ggplot2::element_blank(),
      axis.text.y = if (show_y) ggplot2::element_text(size = 8, color = "grey20") else ggplot2::element_blank(),
      axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      legend.spacing.y = grid::unit(1, "mm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(t = 4, r = 5, b = 0, l = 2),
      strip.clip = "off",
      strip.placement = "outside",
      strip.text.x = ggplot2::element_text(angle = 45, hjust = 0, size = 8, margin = ggplot2::margin(b = 2)),
      strip.text.x.bottom = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2))
    )
}

#' Order GWAS rows and clusters for score heatmaps
#'
#' Rows follow GWAS category and ID. Clusters are grouped into the first
#' matching compartment pattern (hyphens and underscores match each other;
#' unmatched clusters form `Other`) and hierarchically ordered by their scores
#' within each compartment.
#' @keywords internal
get_ordered_GWAS_score_plot_data <- function(data_per_GWAS_and_cluster_df, compartments_patterns, score_col = "median_score") {
  row_metadata <- data_per_GWAS_and_cluster_df |>
    dplyr::distinct(Category, GWAS_ID) |>
    dplyr::distinct(GWAS_ID, .keep_all = TRUE)
  score_mat <- data_per_GWAS_and_cluster_df |>
    dplyr::select(GWAS_ID, cluster, dplyr::all_of(score_col)) |>
    tidyr::pivot_wider(names_from = cluster, values_from = dplyr::all_of(score_col), values_fill = 0) |>
    tibble::column_to_rownames("GWAS_ID") |>
    as.matrix()

  clusters <- colnames(score_mat)
  compartment <- if (is.null(compartments_patterns)) {
    "Cluster"
  } else {
    patterns <- stringr::str_replace_all(unlist(compartments_patterns), "-", "_")
    matched <- vapply(stringr::str_replace_all(clusters, "-", "_"), \(cluster) {
      hit <- which(stringr::str_detect(cluster, patterns))
      if (length(hit)) names(compartments_patterns)[[hit[[1]]]] else "Other"
    }, character(1))
    factor(unname(matched), levels = c(names(compartments_patterns), "Other"))
  }
  cluster_metadata <- tibble::tibble(cluster = clusters, compartment = compartment)
  cluster_levels <- cluster_metadata |>
    dplyr::arrange(compartment, cluster) |>
    dplyr::group_split(compartment, .keep = FALSE) |>
    purrr::map(\(group_df) {
      group_clusters <- group_df$cluster
      if (length(group_clusters) > 2) {
        group_clusters[stats::hclust(stats::dist(t(score_mat[, group_clusters, drop = FALSE])))$order]
      } else {
        group_clusters
      }
    }) |>
    unlist(use.names = FALSE)
  row_levels <- rev(dplyr::arrange(row_metadata, Category, GWAS_ID)$GWAS_ID)
  cluster_support <- data_per_GWAS_and_cluster_df |>
    dplyr::select(cluster, dplyr::any_of(c("n_cells", "n_counts", "n_features", "counts_per_feature"))) |>
    dplyr::summarise(
      dplyr::across(dplyr::everything(), \(x) max(x, na.rm = TRUE)),
      .by = cluster
    )

  list(
    metadata = row_metadata |> dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = row_levels)) |> dplyr::arrange(GWAS_ID),
    scores = data_per_GWAS_and_cluster_df |> dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = row_levels), cluster = factor(cluster, levels = cluster_levels)),
    clusters = cluster_metadata |>
      dplyr::left_join(cluster_support, by = "cluster") |>
      dplyr::mutate(cluster = factor(cluster, levels = cluster_levels)),
    row_levels = row_levels,
    cluster_levels = cluster_levels
  )
}

get_plot_group_breaks <- function(ordered_values) {
  group_lengths <- rle(as.character(ordered_values))$lengths
  if (length(group_lengths) <= 1) {
    return(numeric())
  }
  cumsum(group_lengths)[seq_len(length(group_lengths) - 1)] + 0.5
}

#' Plot a GWAS feature matrix
#'
#' Draw a GWAS-by-feature heatmap or significance-filtered dotplot.
#'
#' @param score_plot_data Ordered score tibble containing GWAS IDs, feature IDs,
#'   and fill values.
#' @param feature_metadata Metadata for plotted features, including compartment
#'   ordering used for vertical separators.
#' @param feature_col Feature/cluster column plotted on the x axis.
#' @param fill_col Numeric column mapped to tile fill.
#' @param fill_label Legend label for the score color scale.
#' @param fill_scale Color-scale type. Use `sequential` for nonnegative scores
#'   and `diverging` for signed deviations.
#' @param fill_limits Optional numeric fill-scale limits.
#' @param support_label_col Optional text column drawn on top of heatmap tiles.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_feature_heatmap <- function(
  score_plot_data,
  feature_metadata,
  feature_col,
  fill_col,
  fill_label,
  fill_scale = c("diverging", "sequential"),
  fill_limits = NULL,
  support_label_col = NULL
) {
  fill_scale <- match.arg(fill_scale)
  row_categories <- score_plot_data |>
    dplyr::distinct(GWAS_ID, Category) |>
    dplyr::mutate(GWAS_ID = as.character(GWAS_ID))
  row_breaks <- get_plot_group_breaks(row_categories$Category[match(levels(score_plot_data$GWAS_ID), row_categories$GWAS_ID)])
  feature_breaks <- feature_metadata |>
    dplyr::arrange(.data[[feature_col]]) |>
    dplyr::pull(compartment) |>
    get_plot_group_breaks()

  feature_plot <- score_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[feature_col]], y = GWAS_ID)) +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0))

  feature_plot <- feature_plot +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data[[fill_col]]),
      color = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey30", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = feature_breaks, color = "grey30", linewidth = 0.35)

  if (!is.null(support_label_col) && support_label_col %in% colnames(score_plot_data)) {
    feature_plot <- feature_plot +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[[support_label_col]]),
        size = 3.4,
        color = "white",
        fontface = "bold",
        na.rm = TRUE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[[support_label_col]]),
        size = 2.8,
        color = "grey10",
        fontface = "bold",
        na.rm = TRUE
      )
  }

  score_guide <- ggplot2::guide_colorbar(
    title.position = "top",
    barwidth = grid::unit(32, "mm"),
    barheight = grid::unit(3, "mm")
  )
  if (fill_scale == "sequential") {
    feature_plot <- feature_plot +
      ggplot2::scale_fill_gradient(
        low = "#F7FBFF",
        high = "#08519C",
        limits = fill_limits,
        name = fill_label,
        guide = score_guide
      )
  } else {
    feature_plot <- feature_plot +
      ggplot2::scale_fill_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        limits = fill_limits,
        name = fill_label,
        guide = score_guide
      )
  }

  compact_x_axis <- dplyr::n_distinct(score_plot_data[[feature_col]]) <= 8
  feature_plot +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(
      xlim = c(0.5, nlevels(score_plot_data[[feature_col]]) + 0.5),
      ylim = c(0.5, nlevels(score_plot_data$GWAS_ID) + 0.5),
      clip = "off"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = if (compact_x_axis) 0 else 45,
        hjust = if (compact_x_axis) 0.5 else 1,
        vjust = if (compact_x_axis) 0.5 else 1,
        color = "grey20"
      ),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.justification = "right",
      legend.box.just = "right",
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 2),
      strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 2))
    )
}

format_GWAS_bar_number <- function(values) {
  vapply(values, \(value) {
    if (!is.finite(value)) {
      return(NA_character_)
    }
    divisor <- if (value >= 1e6) 1e6 else if (value >= 1e3) 1e3 else 1
    suffix <- if (divisor == 1e6) "M" else if (divisor == 1e3) "K" else ""
    paste0(format(signif(value / divisor, 3), trim = TRUE, scientific = FALSE), suffix)
  }, character(1))
}

plot_GWAS_feature_support_tracks <- function(feature_metadata, feature_col = "cluster") {
  support_plot_data <- feature_metadata |>
    dplyr::arrange(.data[[feature_col]]) |>
    dplyr::mutate(
      label = format_GWAS_bar_number(n_cells)
    )

  support_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[feature_col]], y = n_cells)) +
    ggplot2::geom_col(fill = "#6B7280", color = "white", linewidth = 0.3, width = 1, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      angle = 90,
      hjust = -0.15,
      color = "grey20",
      size = 2.8,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      position = "right",
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = ggplot2::expansion(mult = c(0, 0.8))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = "Nuclei") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(angle = 90),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 2, r = 5, b = 0, l = 2)
    )
}

#' Plot GWAS metadata tracks
#'
#' Plot GWAS category, finemapping method, sample-size, loci, and ancestry tracks.
#'
#' @param ordered_metadata GWAS metadata tibble already ordered/factored by
#'   `GWAS_ID`, with category, method, loci, sample-size, and ancestry columns.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_GWAS_metadata_tracks <- function(ordered_metadata) {
  method_colors <- c("SuSie" = "#238B45", "SuSiE-inf" = "#41B6C4", "PICS" = "#F16913")
  unknown_methods <- setdiff(
    sort(unique(stats::na.omit(ordered_metadata$finemappingMethod))),
    names(method_colors)
  )
  method_colors <- c(
    method_colors,
    make_named_heatmap_palette(unknown_methods, palette = "Dark2")
  )
  row_levels <- levels(ordered_metadata$GWAS_ID)
  row_categories <- ordered_metadata |>
    dplyr::distinct(GWAS_ID, Category) |>
    dplyr::mutate(GWAS_ID = as.character(GWAS_ID))
  row_breaks <- get_plot_group_breaks(row_categories$Category[match(row_levels, row_categories$GWAS_ID)])
  bar_plot_data <- ordered_metadata |>
    dplyr::transmute(GWAS_ID, Loci = n_credible_set_loci, Samples = sample_size) |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "track", values_to = "value") |>
    dplyr::mutate(label = format_GWAS_bar_number(value))
  # Open Targets and local metadata both carry all six ancestry fractions.
  ancestry_columns <- grep("^ancestry_(EUR|EAS|AFR|AMR|SAS|OTH)$", colnames(ordered_metadata), value = TRUE)
  ancestry_reported <- rowSums(!is.na(ordered_metadata[ancestry_columns])) > 0
  ancestry_plot_data <- ordered_metadata[ancestry_reported, , drop = FALSE] |>
    dplyr::select(GWAS_ID, dplyr::all_of(ancestry_columns)) |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "ancestry_group", values_to = "fraction") |>
    dplyr::filter(!is.na(fraction)) |>
    dplyr::mutate(ancestry_group = stringr::str_remove(ancestry_group, "^ancestry_")) |>
    dplyr::bind_rows(
      tibble::tibble(
        GWAS_ID = ordered_metadata$GWAS_ID[!ancestry_reported],
        ancestry_group = "Unreported",
        fraction = 1
      )
    )

  category_plot <- ordered_metadata |>
    ggplot2::ggplot(ggplot2::aes(x = 1, y = GWAS_ID, fill = Category)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(
      labels = \(x) stringr::str_replace_all(x, "_", " "),
      drop = FALSE,
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(
      values = make_named_heatmap_palette(ordered_metadata$Category, "Set3"),
      labels = \(x) stringr::str_replace_all(x, "_", " "),
      name = "Category",
      guide = ggplot2::guide_legend(ncol = min(4L, dplyr::n_distinct(ordered_metadata$Category)), byrow = TRUE, title.position = "top", order = 1)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme(show_y = TRUE)

  method_plot <- ordered_metadata |>
    ggplot2::ggplot(ggplot2::aes(x = 1, y = GWAS_ID, fill = finemappingMethod)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = method_colors,
      name = "Method",
      guide = ggplot2::guide_legend(ncol = 1, title.position = "top", order = 2)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme()

  bar_plot <- bar_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = GWAS_ID)) +
    ggplot2::geom_col(fill = "grey45", width = 0.8, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust = -0.1,
      color = "grey20",
      size = 3,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::facet_grid(. ~ track, scales = "free_x", switch = "x") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.6))) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme()

  ancestry_plot <- ancestry_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = fraction, y = GWAS_ID, fill = ancestry_group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_x_continuous(breaks = c(0, 1), labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = c(
        "EUR" = "#4DAF4A",
        "EAS" = "#377EB8",
        "AFR" = "#984EA3",
        "AMR" = "#FF7F00",
        "SAS" = "#E41A1C",
        "OTH" = "#999999",
        "Unreported" = "#D9D9D9"
      ),
      name = "Ancestry",
      guide = ggplot2::guide_legend(ncol = 2, byrow = TRUE, title.position = "top", order = 3)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme(show_x = TRUE)

  tracks <- lapply(list(category_plot, method_plot, bar_plot, ancestry_plot), \(plot) {
    plot + ggplot2::coord_cartesian(ylim = c(0.5, length(row_levels) + 0.5), clip = "off")
  })
  patchwork::wrap_plots(tracks, nrow = 1, widths = c(0.48, 0.48, 1.25, 1.15), guides = "collect") &
    ggplot2::theme(legend.justification = "left", legend.box.just = "left")
}

#' Plot GWAS scores by cluster
#'
#' Combine GWAS metadata tracks with an ordered cluster heatmap or dotplot.
#'
#' @param data_per_GWAS_and_cluster_df Score summary tibble with one row per
#'   GWAS/cluster combination.
#' @param GWAS_metadata_tracks_plot Four-track patchwork returned by
#'   `plot_GWAS_metadata_tracks()`, aligned to the same GWAS ordering.
#' @param compartments_patterns Optional named regex patterns used to group
#'   clusters into compartments.
#' @param fill_col Numeric column mapped to heatmap fill.
#' @param fill_label Legend label for the heatmap fill.
#' @param fill_scale Color-scale type passed to `plot_GWAS_feature_heatmap()`.
#' @param fill_limits Optional numeric fill-scale limits.
#' @param support_label_col Optional text column drawn on top of heatmap tiles.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_by_cluster_heatmap <- function(
  data_per_GWAS_and_cluster_df,
  GWAS_metadata_tracks_plot,
  compartments_patterns = NULL,
  fill_col = "median_score",
  fill_label = "Score",
  fill_scale = c("diverging", "sequential"),
  fill_limits = NULL,
  support_label_col = NULL
) {
  if (nrow(data_per_GWAS_and_cluster_df) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }
  fill_scale <- match.arg(fill_scale)

  ordered_data <- get_ordered_GWAS_score_plot_data(
    data_per_GWAS_and_cluster_df,
    compartments_patterns,
    score_col = fill_col
  )

  score_plot <- plot_GWAS_feature_heatmap(
    score_plot_data = ordered_data$scores,
    feature_metadata = ordered_data$clusters,
    feature_col = "cluster",
    fill_col = fill_col,
    fill_label = fill_label,
    fill_scale = fill_scale,
    fill_limits = fill_limits,
    support_label_col = support_label_col
  )

  # Keep all body panels at the same layout level so axes and legends cannot
  # independently resize the metadata rows or the score matrix.
  plots <- lapply(seq_len(4L), \(i) GWAS_metadata_tracks_plot[[i]])
  plots[[5]] <- score_plot
  design <- "ABCDE"
  heights <- grid::unit(1, "null")
  if ("n_cells" %in% colnames(ordered_data$clusters) && any(!is.na(ordered_data$clusters$n_cells))) {
    plots[[6]] <- plot_GWAS_feature_support_tracks(ordered_data$clusters)
    design <- "####F\nABCDE"
    heights <- grid::unit(c(28, 1), c("mm", "null"))
  }

  patchwork::wrap_plots(
    plots, design = design, widths = c(0.35, 0.35, 2, 1.2, 9),
    heights = heights, guides = "collect"
  ) & ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "top",
    legend.margin = ggplot2::margin(4, 8, 4, 0)
  )
}

#' Save a GWAS-by-group heatmap with a height that grows with the GWAS count
#' @keywords internal
save_GWAS_heatmap <- function(plot, GWAS_IDs) {
  save_plots_structured(plot, filetype = "png", width = 17, height = max(5.5, 0.3 * dplyr::n_distinct(GWAS_IDs) + 3.5))
}

#' Plot one WNN SCAVENGE median-TRS heatmap
#' @param heatmap_data Cached WNN group scores, optionally with `support_label`
#'   tile labels and the `support_note` that explains them.
#' @param GWAS_metadata_tracks_plot Shared GWAS metadata tracks.
#' @param compartments_patterns Optional cell-type compartment patterns.
#' @param grouping WNN cell types or named clusters.
#' @param scaled Min-max scale medians within each GWAS for display only.
#' @return One annotated heatmap, or an empty plot list.
plot_WNN_TRS_heatmap <- function(heatmap_data, GWAS_metadata_tracks_plot,
  compartments_patterns = NULL, grouping = c("cell_types", "clusters"), scaled = FALSE) {
  grouping <- match.arg(grouping)
  group_column <- if (grouping == "cell_types") "WNN_harmony_SNN_cluster_cell_type" else "WNN_harmony_SNN_cluster_named"
  data <- heatmap_data |>
    dplyr::filter(grouping_col == group_column) |>
    dplyr::select(-grouping_col)
  if (nrow(data) == 0L) return(structure(list(), class = c("empty_plot_list", "list")))
  if (scaled) data <- data |>
    dplyr::mutate(median_score = min_max_scale_vec(median_score), .by = GWAS_ID)
  plot_GWAS_by_cluster_heatmap(data,
    GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
    compartments_patterns = compartments_patterns,
    fill_label = if (scaled) "Relative median TRS" else "Median TRS",
    fill_scale = "sequential", fill_limits = if (scaled) c(0, 1) else NULL,
    support_label_col = "support_label"
  ) + patchwork::plot_annotation(
    title = paste(if (scaled) "Relative" else "Median", "SCAVENGE trait scores by WNN",
      if (grouping == "cell_types") "cell type" else "cluster"),
    subtitle = if (scaled) "Compare groups within each trait; colours do not measure absolute differences between traits."
      else "Look for groups with high median trait scores.",
    caption = stringr::str_wrap(paste(c(
      if (scaled) "Colour: median cell-level TRS min-max scaled to 0-1 across the displayed WNN groups within each GWAS."
      else "Colour: median cell-level TRS within each WNN group, without display scaling.",
      "All available scores are shown.",
      unique(data[["support_note"]])
    ), collapse = " "), 150)
  )
}

#' Plot ordinary or absolute-effect chromVAR deviations from cached summaries
#' @param deviation_tibble Cached scores containing deviation, relative_deviation and support_label.
#' @param GWAS_metadata_tracks_plot Metadata tracks for the same GWAS collection.
#' @param compartments_patterns Optional cell-type compartment patterns.
#' @param standardize Use the precomputed within-GWAS standardized deviation column.
#' @param beta_weighted Label the input as PIP times absolute-effect weighted.
#' @return Annotated heatmap, or an explicit empty-result plot.
plot_GWAS_chromVAR_deviation_heatmap <- function(deviation_tibble,
  GWAS_metadata_tracks_plot, compartments_patterns = NULL,
  standardize = FALSE, beta_weighted = FALSE) {
  if (nrow(deviation_tibble) == 0L) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No eligible cell-type GWAS scores"))
  }
  quantity <- if (standardize) "Standardized deviation" else "Deviation"
  plot_GWAS_by_cluster_heatmap(
    deviation_tibble,
    GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
    compartments_patterns = compartments_patterns,
    fill_col = if (standardize) "relative_deviation" else "deviation",
    fill_label = paste(if (beta_weighted) "PIP × |β|" else "Ordinary", quantity),
    support_label_col = "support_label"
  ) + patchwork::plot_annotation(
    title = paste(if (beta_weighted) "Effect-magnitude-weighted" else "Ordinary",
      "GWAS-linked accessibility by cell type", if (standardize) "— standardized" else "— unscaled"),
    subtitle = paste(
      if (standardize) "Compare cell types within each trait; red indicates above-average deviation."
      else "Positive values indicate accessibility above the weighted background expectation.",
      "Stars show chromVAR z-score support, not adjusted-p significance or independent donor evidence.", sep = "\n"),
    caption = stringr::str_wrap(paste(
      "ATAC counts are summed by WNN cell type.",
      if (beta_weighted) paste(
        "PIP × |β| weighting uses variant effects when complete, otherwise the highest-PIP available effect per locus; only eligible GWAS are included.",
        "Genetic effect direction is discarded. Weights preserve retained PIP mass before peak weights are capped at one."),
      if (standardize) "Colour: deviation centred and divided by its SD across cell types within each GWAS."
      else "Colour: weighted observed-minus-expected accessibility divided by the depth-adjusted expectation, without centring or SD scaling across cell types.",
      "Stars: * Z >= 1.645; ** Z >= 2.326 (unadjusted upper-tail normal P <= 0.05 and <= 0.01) against the betterChromVAR background for the corresponding weighting. Nuclei counts describe input support; this is a descriptive pooled comparison."
    ), width = 150)
  )
}

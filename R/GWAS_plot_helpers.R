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
      caption = stringr::str_glue(
        "Dashed line in the frac_overlapped_peaks facet marks {overlap_threshold} overlapped peaks ({scales::percent(overlap_threshold_fraction)} of all peaks)."
      )
    )
}

scale_GWAS_heatmap_scores <- function(heatmap_data, group_cols = character(), score_col = "median_score") {
  if (nrow(heatmap_data) == 0) {
    return(heatmap_data)
  }

  heatmap_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("GWAS_ID", group_cols)))) |>
    dplyr::mutate(!!score_col := min_max_scale_vec(.data[[score_col]])) |>
    dplyr::ungroup()
}

#' Add cluster-level permutation significance to SCAVENGE dotplot data
#'
#' Classify visible point sizes from within-grouping BH-adjusted empirical
#' P-values for the cluster median TRS statistic.
#'
#' @param dotplot_data SCAVENGE group-summary data containing
#'   `permutation_p_adj_BH`.
#' @return `dotplot_data` with a `significance` factor. Values above 0.05 are
#'   `NA` and therefore omitted from the dotplot.
#' @keywords internal

add_SCAVENGE_dotplot_significance <- function(dotplot_data) {
  dotplot_data |>
    dplyr::mutate(
      significance = factor(
        dplyr::case_when(
          permutation_p_adj_BH < 0.01 ~ "P < 0.01",
          permutation_p_adj_BH <= 0.05 ~ "P <= 0.05"
        ),
        levels = c("P <= 0.05", "P < 0.01")
      )
    )
}

#' Get SCAVENGE TRS UMAP plots
#'
#' Build UMAP overlay plots from cell-level SCAVENGE TRS scores.
#'
#' @param TRS_tibble Cell-level TRS tibble with `barcode_w_prefix`, `GWAS_ID`,
#'   `score`, and significance columns.
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
    ggplot2::labs(subtitle = NULL)

  if (!is.null(label_col)) {
    if (!is.character(label_col) || length(label_col) != 1 || is.na(label_col) || !nzchar(label_col)) {
      stop("`label_col` must be NULL or a non-empty length-1 character vector.", call. = FALSE)
    }
    if (!label_col %in% colnames(SCAVENGE_metadata_tibble)) {
      stop("Required label column not found: ", label_col, call. = FALSE)
    }

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

plot_SCAVENGE_summary_sig_proportion <- function(summary_tibble) {
  if (nrow(summary_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot_tibble <- summary_tibble |>
    dplyr::mutate(
      facet_id = stringr::str_c(grouping_col, GWAS_ID, sep = "___"),
      group_reorder = tidytext::reorder_within(cluster, prop_sig, facet_id)
    )

  split(plot_tibble, plot_tibble$grouping_col) |>
    purrr::map(\(group_tibble) {
      ggplot2::ggplot(group_tibble, ggplot2::aes(x = group_reorder, y = prop_sig, fill = cluster)) +
        ggplot2::facet_wrap(~GWAS_ID, scales = "free") +
        tidytext::scale_x_reordered(labels = function(x) gsub("___.*$", "", x)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(x = "", y = "Proportion of enriched cells") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
    })
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
        ggplot2::labs(x = "", y = "SCAVENGE TRS") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
    })
}

add_GWAS_heatmap_categories <- function(heatmap_data, GWAS_tibble) {
  heatmap_data |>
    dplyr::inner_join(GWAS_tibble |> dplyr::select(GWAS_ID, Category, dplyr::any_of("variant_weighting_mode")), by = "GWAS_ID")
}

#' Assign compartment
#'
#' Assign categorical compartments from regex patterns applied to type labels.
#'
#' @param in_tibble Input tibble containing `type_col`.
#' @param patterns Named character vector of regex patterns; names become
#'   compartment labels.
#' @param type_col Column whose values are matched after replacing `-` with `_`.
#' @param new_col Name of the output compartment column.
#' @return `in_tibble` with `new_col` added; unmatched rows receive `NA`.
#' @keywords internal

assign_compartment <- function(in_tibble, patterns, type_col, new_col = "compartment") {
  match_col <- paste0(type_col, "__compartment_match")
  patterns <- stringr::str_replace_all(patterns, "-", "_")
  conds <- purrr::map2(
    patterns,
    names(patterns),
    ~ rlang::expr(stringr::str_detect(.data[[match_col]], !!.x) ~ !!.y)
  )

  in_tibble |>
    dplyr::mutate(!!match_col := stringr::str_replace_all(as.character(.data[[type_col]]), "-", "_")) |>
    dplyr::mutate(
      !!new_col := dplyr::case_when(
        !!!conds,
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::select(-dplyr::all_of(match_col))
}

#' Get compartment metadata
#'
#' Build one-row-per-value metadata with ordered compartment labels.
#'
#' @param values Character/factor values to de-duplicate into metadata rows.
#' @param compartments_patterns Optional named regex vector used by
#'   `assign_compartment()`.
#' @param type_col Name of the value column in the returned metadata.
#' @param default_compartment Compartment label used when no pattern mapping is
#'   supplied.
#' @return Metadata tibble with `type_col` and factor `compartment`.
#' @keywords internal

get_compartment_metadata <- function(values, compartments_patterns, type_col, default_compartment) {
  metadata <- tibble::tibble(!!type_col := as.character(values)) |>
    dplyr::distinct()

  if (is.null(compartments_patterns)) {
    return(metadata |> dplyr::mutate(compartment = default_compartment))
  }

  metadata |>
    assign_compartment(compartments_patterns, type_col = type_col) |>
    dplyr::mutate(
      compartment = dplyr::coalesce(compartment, "Other"),
      compartment = factor(compartment, levels = c(names(compartments_patterns), "Other"))
    )
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

get_heatmap_legend_ncol <- function(values, max_row_chars = 70) {
  values <- sort(unique(stats::na.omit(as.character(values))))
  if (length(values) == 0) {
    return(1L)
  }

  label_widths <- nchar(values) + 6
  for (ncol in seq.int(length(values), 1L)) {
    nrow <- ceiling(length(values) / ncol)
    row_widths <- vapply(seq_len(nrow), \(row_idx) sum(label_widths[seq(row_idx, length(values), by = nrow)]), numeric(1))
    if (max(row_widths) <= max_row_chars) {
      return(ncol)
    }
  }
  1L
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

#' Order GWAS score plot ids
#'
#' Order GWAS rows and feature columns for clustered score heatmaps.
#'
#' @param score_data GWAS-by-feature score tibble with GWAS ID, cluster/feature,
#'   and score columns.
#' @param row_metadata GWAS metadata tibble used to annotate and order GWAS rows.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param score_col Numeric score column used for clustering/order calculation.
#' @param compartments_patterns Optional named regex patterns used to group
#'   feature columns into compartments before plotting.
#' @return List with ordered GWAS IDs, ordered feature IDs, and feature metadata.
#' @keywords internal

order_GWAS_score_plot_ids <- function(score_data, row_metadata, cluster_col = "cluster", score_col = "median_score", compartments_patterns = NULL) {
  score_mat <- score_data |>
    dplyr::select(GWAS_ID, dplyr::all_of(c(cluster_col, score_col))) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(cluster_col), values_from = dplyr::all_of(score_col), values_fill = 0) |>
    tibble::column_to_rownames("GWAS_ID") |>
    as.matrix()

  cluster_metadata <- get_compartment_metadata(
    colnames(score_mat),
    compartments_patterns,
    type_col = cluster_col,
    default_compartment = "Cluster"
  )

  cluster_order <- cluster_metadata |>
    dplyr::arrange(compartment, .data[[cluster_col]]) |>
    dplyr::group_split(compartment, .keep = FALSE) |>
    purrr::map(\(group_df) {
      clusters <- group_df[[cluster_col]]
      if (length(clusters) > 2) {
        clusters[stats::hclust(stats::dist(t(score_mat[, clusters, drop = FALSE])))$order]
      } else {
        clusters
      }
    }) |>
    unlist(use.names = FALSE)

  row_order <- row_metadata |>
    dplyr::arrange(Category, GWAS_ID) |>
    dplyr::pull(GWAS_ID)

  list(row_order = row_order, cluster_order = cluster_order, cluster_metadata = cluster_metadata)
}

get_ordered_GWAS_score_plot_data <- function(data_per_GWAS_and_cluster_df, compartments_patterns, score_col = "median_score") {
  row_metadata <- data_per_GWAS_and_cluster_df |>
    dplyr::distinct(Category, GWAS_ID) |>
    dplyr::distinct(GWAS_ID, .keep_all = TRUE)
  axes <- order_GWAS_score_plot_ids(
    data_per_GWAS_and_cluster_df,
    row_metadata,
    score_col = score_col,
    compartments_patterns = compartments_patterns
  )
  row_levels <- rev(axes$row_order)
  cluster_levels <- axes$cluster_order
  cluster_support <- data_per_GWAS_and_cluster_df |>
    dplyr::select(cluster, dplyr::any_of(c("n_cells", "n_counts", "n_features", "counts_per_feature"))) |>
    dplyr::summarise(
      dplyr::across(dplyr::everything(), \(x) max(x, na.rm = TRUE)),
      .by = cluster
    )

  list(
    metadata = row_metadata |> dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = row_levels)) |> dplyr::arrange(GWAS_ID),
    scores = data_per_GWAS_and_cluster_df |> dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = row_levels), cluster = factor(cluster, levels = cluster_levels)),
    clusters = axes$cluster_metadata |>
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
#' @param fill_col Numeric column mapped to tile fill or point color.
#' @param fill_label Legend label for the score color scale.
#' @param fill_midpoint Midpoint for the diverging fill scale.
#' @param fill_scale Color-scale type. Use `sequential` for nonnegative scores
#'   and `diverging` for signed deviations.
#' @param fill_limits Optional numeric fill-scale limits.
#' @param title Optional plot title.
#' @param support_label_col Optional text column drawn on top of heatmap tiles.
#' @param point_size_col Optional factor column mapped to point size. When set,
#'   draw a red sequential dotplot instead of heatmap tiles.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_feature_heatmap <- function(
  score_plot_data,
  feature_metadata,
  feature_col,
  fill_col,
  fill_label,
  fill_midpoint = 0,
  fill_scale = c("diverging", "sequential"),
  fill_limits = NULL,
  title = NULL,
  support_label_col = NULL,
  point_size_col = NULL
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

  if (is.null(point_size_col)) {
    feature_plot <- feature_plot +
      ggplot2::geom_tile(
        ggplot2::aes(fill = .data[[fill_col]]),
        color = "white",
        linewidth = 0.3
      )
  } else {
    feature_plot <- feature_plot +
      ggplot2::geom_point(
        data = score_plot_data |> dplyr::filter(!is.na(.data[[point_size_col]])),
        ggplot2::aes(color = .data[[fill_col]], size = .data[[point_size_col]]),
        alpha = 0.95,
        na.rm = TRUE
      ) +
      ggplot2::scale_size_manual(
        values = c("P <= 0.05" = 2.4, "P < 0.01" = 4.8),
        drop = FALSE,
        name = "BH-adjusted P-value",
        guide = ggplot2::guide_legend(title.position = "top")
      )
  }

  feature_plot <- feature_plot +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey30", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = feature_breaks, color = "grey30", linewidth = 0.35)

  if (is.null(point_size_col) && !is.null(support_label_col) && support_label_col %in% colnames(score_plot_data)) {
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
  if (!is.null(point_size_col)) {
    feature_plot <- feature_plot +
      ggplot2::scale_color_gradient(
        low = "#FEE5D9",
        high = "#A50F15",
        limits = fill_limits,
        name = fill_label,
        guide = score_guide
      )
  } else if (fill_scale == "sequential") {
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
        midpoint = fill_midpoint,
        limits = fill_limits,
        name = fill_label,
        guide = score_guide
      )
  }

  compact_x_axis <- dplyr::n_distinct(score_plot_data[[feature_col]]) <= 8
  feature_plot +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::coord_cartesian(clip = "off") +
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
  if (!"n_cells" %in% colnames(feature_metadata) || all(is.na(feature_metadata$n_cells))) {
    return(patchwork::plot_spacer())
  }

  support_plot_data <- feature_metadata |>
    dplyr::arrange(.data[[feature_col]]) |>
    dplyr::mutate(
      panel_max = max(n_cells, na.rm = TRUE),
      label = format_GWAS_bar_number(n_cells),
      label_inside = n_cells / panel_max >= 0.28,
      label_y = dplyr::if_else(label_inside, n_cells - 0.025 * panel_max, n_cells + 0.025 * panel_max),
      label_hjust = dplyr::if_else(label_inside, 1, 0),
      label_color = dplyr::if_else(label_inside, "white", "grey20")
    )

  support_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[feature_col]], y = n_cells)) +
    ggplot2::geom_col(fill = "#6B7280", color = "white", linewidth = 0.3, width = 1, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(y = label_y, label = label, hjust = label_hjust, color = label_color),
      angle = 90,
      size = 2.8,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      position = "right",
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::scale_color_identity() +
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
  row_levels <- if (is.factor(ordered_metadata$GWAS_ID)) levels(ordered_metadata$GWAS_ID) else unique(as.character(ordered_metadata$GWAS_ID))
  row_categories <- ordered_metadata |>
    dplyr::distinct(GWAS_ID, Category) |>
    dplyr::mutate(GWAS_ID = as.character(GWAS_ID))
  row_breaks <- get_plot_group_breaks(row_categories$Category[match(row_levels, row_categories$GWAS_ID)])
  bar_plot_data <- ordered_metadata |>
    dplyr::transmute(GWAS_ID, Loci = n_credible_set_loci, Samples = sample_size) |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "track", values_to = "value") |>
    dplyr::mutate(
      panel_max = max(value, na.rm = TRUE),
      label = format_GWAS_bar_number(value),
      label_inside = value / panel_max >= 0.28,
      label_x = dplyr::if_else(label_inside, value - 0.025 * panel_max, value + 0.025 * panel_max),
      label_hjust = dplyr::if_else(label_inside, 1, 0),
      label_color = dplyr::if_else(label_inside, "white", "grey20"),
      .by = track
    )
  ancestry_columns <- grep(
    "^ancestry_(EUR|EAS|AFR|AMR|SAS|OTH)$",
    colnames(ordered_metadata),
    value = TRUE
  )
  ancestry_reported <- if (length(ancestry_columns) == 0L) {
    rep(FALSE, nrow(ordered_metadata))
  } else {
    rowSums(!is.na(ordered_metadata[ancestry_columns])) > 0
  }
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
      guide = ggplot2::guide_legend(ncol = 1, title.position = "top")
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
      guide = ggplot2::guide_legend(ncol = 1, title.position = "top")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme()

  bar_plot <- bar_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = GWAS_ID)) +
    ggplot2::geom_col(fill = "grey45", width = 0.8, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(x = label_x, label = label, hjust = label_hjust, color = label_color),
      size = 3,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::facet_grid(. ~ track, scales = "free_x", switch = "x") +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
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
      guide = ggplot2::guide_legend(ncol = min(2L, get_heatmap_legend_ncol(ancestry_plot_data$ancestry_group, max_row_chars = 35)), byrow = TRUE, title.position = "top")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme(show_x = TRUE)

  patchwork::wrap_plots(category_plot, method_plot, bar_plot, ancestry_plot, nrow = 1, widths = c(0.48, 0.48, 1.25, 1.15), guides = "collect") &
    ggplot2::theme(legend.justification = "left", legend.box.just = "left")
}

#' Plot GWAS scores by cluster
#'
#' Combine GWAS metadata tracks with an ordered cluster heatmap or dotplot.
#'
#' @param data_per_GWAS_and_cluster_df Score summary tibble with one row per
#'   GWAS/cluster combination.
#' @param GWAS_metadata_tracks_plot Patchwork/ggplot metadata track aligned to
#'   the same GWAS ordering.
#' @param compartments_patterns Optional named regex patterns used to group
#'   clusters into compartments.
#' @param scaled Logical; when `TRUE`, use 0.5 as the diverging color midpoint
#'   for min-max scaled scores.
#' @param fill_col Numeric column mapped to heatmap fill.
#' @param fill_label Legend label for the heatmap fill.
#' @param fill_scale Color-scale type passed to `plot_GWAS_feature_heatmap()`.
#' @param fill_limits Optional numeric fill-scale limits.
#' @param support_label_col Optional text column drawn on top of heatmap tiles.
#' @param point_size_col Optional factor column mapped to dot size.
#' @param show_feature_support Logical; when `TRUE`, draw a nuclei-count support
#'   annotation if `n_cells` is available.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_by_cluster_heatmap <- function(
  data_per_GWAS_and_cluster_df,
  GWAS_metadata_tracks_plot,
  compartments_patterns = NULL,
  scaled = FALSE,
  fill_col = "median_score",
  fill_label = "Score",
  fill_scale = c("diverging", "sequential"),
  fill_limits = NULL,
  support_label_col = NULL,
  point_size_col = NULL,
  show_feature_support = TRUE
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
    fill_midpoint = if (scaled) 0.5 else 0,
    fill_scale = fill_scale,
    fill_limits = fill_limits,
    support_label_col = support_label_col,
    point_size_col = point_size_col
  )

  if (isTRUE(show_feature_support) && "n_cells" %in% colnames(ordered_data$clusters) && any(!is.na(ordered_data$clusters$n_cells))) {
    left_panel <- patchwork::plot_spacer() / GWAS_metadata_tracks_plot + patchwork::plot_layout(heights = c(0.14, 1))
    right_panel <- plot_GWAS_feature_support_tracks(ordered_data$clusters) / score_plot + patchwork::plot_layout(heights = c(0.14, 1))
    return(patchwork::wrap_plots(left_panel, right_panel, nrow = 1, widths = c(4.8, 9)))
  }

  patchwork::wrap_plots(GWAS_metadata_tracks_plot, score_plot, nrow = 1, widths = c(4.8, 9))
}

#' Plot grouped GWAS scores by cluster
#'
#' Split GWAS cluster heatmaps or dotplots by a grouping column and return a
#' named plot list.
#'
#' @param data_per_GWAS_and_cluster_df Score summary tibble containing `split_col`
#'   in addition to GWAS and cluster score fields.
#' @param split_col Column used to split the data into one heatmap per value.
#' @param GWAS_metadata_tracks_plot Patchwork/ggplot metadata track aligned to
#'   the same GWAS ordering.
#' @param name_suffix Optional suffix appended to names of returned plot-list
#'   elements.
#' @param compartments_patterns Optional named regex patterns used to group
#'   clusters into compartments.
#' @param scaled Logical; when `TRUE`, use scaled-score color midpoint behavior.
#' @param fill_label Legend label for the heatmap fill.
#' @param fill_scale Color-scale type passed to `plot_GWAS_by_cluster_heatmap()`.
#' @param fill_limits Optional numeric fill-scale limits.
#' @param support_label_col Optional text column drawn on top of heatmap tiles.
#' @param point_size_col Optional factor column mapped to dot size.
#' @param caption Optional caption added below each grouped heatmap.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_grouped_GWAS_by_cluster_heatmaps <- function(
  data_per_GWAS_and_cluster_df,
  split_col,
  GWAS_metadata_tracks_plot,
  name_suffix = NULL,
  compartments_patterns = NULL,
  scaled = FALSE,
  fill_label = "Score",
  fill_scale = c("diverging", "sequential"),
  fill_limits = NULL,
  support_label_col = NULL,
  point_size_col = NULL,
  caption = NULL
) {
  if (nrow(data_per_GWAS_and_cluster_df) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }
  fill_scale <- match.arg(fill_scale)

  grouped_plot_data <- dplyr::group_by(data_per_GWAS_and_cluster_df, .data[[split_col]])
  plot_names <- dplyr::group_keys(grouped_plot_data)[[split_col]]
  if (!is.null(name_suffix)) {
    plot_names <- stringr::str_c(plot_names, "_", name_suffix)
  }

  grouped_plot_data |>
    dplyr::group_split() |>
    purrr::set_names(plot_names) |>
    purrr::map(\(group_data) {
      plot <- plot_GWAS_by_cluster_heatmap(
        group_data |> dplyr::select(-dplyr::all_of(split_col)),
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = compartments_patterns,
        scaled = scaled,
        fill_label = fill_label,
        fill_scale = fill_scale,
        fill_limits = fill_limits,
        support_label_col = support_label_col,
        point_size_col = point_size_col
      )
      if (!is.null(caption)) {
        plot <- plot + patchwork::plot_annotation(caption = caption)
      }
      plot
    })
}

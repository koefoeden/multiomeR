#' Append per-action QC counts; overlapping exclusions go to the first listed action.
summarize_QC_cell_retention <- function(before_metadata, after_metadata, GEM_well_IDs,
                                      stage, discarded_barcodes, previous_stages = NULL) {
  if (length(discarded_barcodes) == 0L) discarded_barcodes <- list("No exclusions" = character())
  if (is.null(names(discarded_barcodes)) || anyDuplicated(names(discarded_barcodes))) {
    stop("Retention discard actions must have unique names.")
  }
  remaining <- before_metadata |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", "GEM_well_ID")))
  if (anyDuplicated(remaining$barcode_w_prefix) ||
    !all(after_metadata$barcode_w_prefix %in% remaining$barcode_w_prefix)) {
    stop("Retention stages require unique input barcodes and a retained subset: ", stage)
  }
  action_counts <- vector("list", length(discarded_barcodes))
  for (action_index in seq_along(discarded_barcodes)) {
    action <- names(discarded_barcodes)[action_index]
    barcodes <- discarded_barcodes[[action_index]]
    counts <- tibble::tibble(GEM_well_ID = GEM_well_IDs) |>
      dplyr::left_join(dplyr::count(remaining, .data$GEM_well_ID, name = "input_cells"),
        by = "GEM_well_ID") |>
      dplyr::left_join(remaining |>
        dplyr::filter(.data$barcode_w_prefix %in% barcodes) |>
        dplyr::count(.data$GEM_well_ID, name = "excluded_cells"), by = "GEM_well_ID") |>
      dplyr::mutate(
        input_cells = tidyr::replace_na(.data$input_cells, 0L),
        excluded_cells = tidyr::replace_na(.data$excluded_cells, 0L),
        retained_cells = .data$input_cells - .data$excluded_cells,
        retained_fraction = .data$retained_cells / dplyr::na_if(.data$input_cells, 0L),
        stage = stage, discard_action = action,
        action_order = action_index
      )
    remaining <- dplyr::filter(remaining, !.data$barcode_w_prefix %in% barcodes)
    action_counts[[action_index]] <- counts
  }
  current_stage <- dplyr::bind_rows(action_counts)
  if (!setequal(remaining$barcode_w_prefix, after_metadata$barcode_w_prefix)) {
    stop("Discard actions do not reproduce the retained barcodes: ", stage)
  }
  if (!is.null(previous_stages)) {
    previous <- previous_stages |>
      dplyr::filter(.data$stage == tail(unique(previous_stages$stage), 1)) |>
      dplyr::filter(.data$action_order == max(.data$action_order))
    current_input <- dplyr::filter(current_stage, .data$action_order == 1L)
    if (!identical(current_input$input_cells,
      previous$retained_cells[match(current_input$GEM_well_ID, previous$GEM_well_ID)])) {
      stop("Retention stages do not connect per GEM well: ", stage)
    }
  }
  dplyr::bind_rows(previous_stages, current_stage)
}

#' Split the two losses between clustered metadata and its scDblFinder-filtered subset.
get_doublet_QC_discarded_barcodes <- function(clustered_metadata, retained_metadata,
                                           scDblFinder_results_df, class_col,
                                           remove_called_doublets) {
  barcodes <- if ("barcode_w_prefix" %in% names(scDblFinder_results_df)) {
    scDblFinder_results_df$barcode_w_prefix
  } else rownames(scDblFinder_results_df)
  list(
    "Called doublets" = if (remove_called_doublets) {
      barcodes[which(scDblFinder_results_df[[class_col]] == "doublet")]
    } else character(),
    "High-doublet clusters" = setdiff(clustered_metadata$barcode_w_prefix,
      retained_metadata$barcode_w_prefix)
  )
}

#' Plot count-proportional continuing and distinct discard branches at each checkpoint.
plot_QC_cell_retention <- function(retention_tibble) {
  actions <- retention_tibble |>
    dplyr::group_by(.data$stage, .data$action_order, .data$discard_action) |>
    dplyr::summarize(input_cells = sum(.data$input_cells),
      retained_cells = sum(.data$retained_cells),
      excluded_cells = sum(.data$excluded_cells), .groups = "drop") |>
    dplyr::mutate(stage_index = match(.data$stage, unique(retention_tibble$stage))) |>
    dplyr::arrange(.data$stage_index, .data$action_order)
  stages <- actions |>
    dplyr::group_by(.data$stage_index, .data$stage) |>
    dplyr::summarize(input_cells = dplyr::first(.data$input_cells),
      retained_cells = dplyr::last(.data$retained_cells), .groups = "drop")
  count_scale <- max(c(stages$input_cells, 1))
  branches <- actions |>
    dplyr::filter(.data$excluded_cells > 0) |>
    dplyr::group_by(.data$stage_index) |>
    dplyr::mutate(offset = 0.12 * count_scale * rev(seq_len(dplyr::n())),
      label_y = (.data$input_cells + .data$retained_cells) / 2 + .data$offset) |>
    dplyr::ungroup()
  progress <- seq(0, 1, length.out = 40)
  bend <- progress^2 * (3 - 2 * progress)
  continuing <- purrr::map_dfr(seq_len(nrow(stages)), function(index) {
    tibble::tibble(x = c(index - 1 + progress, rev(index - 1 + progress)),
      y = c(rep(0, length(progress)), rep(stages$retained_cells[index], length(progress))),
      flow = paste(index, "retained"), outcome = "Continuing")
  })
  discarded <- purrr::map_dfr(seq_len(nrow(branches)), function(index) {
    branch <- branches[index, ]
    x <- branch$stage_index - 1 + 0.55 * progress
    tibble::tibble(x = c(x, rev(x)),
      y = c(branch$retained_cells + branch$offset * bend,
        rev(branch$input_cells + branch$offset * bend)),
      flow = paste(index, "excluded"), outcome = "Discarded")
  })
  counts <- tibble::tibble(x = 0:nrow(stages),
    cells = c(stages$input_cells[1], stages$retained_cells))
  ggplot2::ggplot(actions) +
    ggplot2::geom_polygon(data = dplyr::bind_rows(continuing, discarded),
      ggplot2::aes(.data$x, .data$y, group = .data$flow, fill = .data$outcome), alpha = 0.8) +
    ggplot2::geom_text(data = counts,
      ggplot2::aes(x = .data$x, y = -0.055 * count_scale,
        label = scales::comma(.data$cells)), size = 3.8) +
    ggplot2::geom_text(data = branches,
      ggplot2::aes(x = .data$stage_index - 0.42, y = .data$label_y,
        label = paste0(stringr::str_wrap(.data$discard_action, width = 22),
          "\n−", scales::comma(.data$excluded_cells))),
      size = 3.2, hjust = 0, colour = "#A45132") +
    ggplot2::scale_fill_manual(values = c(Continuing = "#348D98", Discarded = "#D58B65")) +
    ggplot2::scale_x_continuous(breaks = 0:nrow(stages),
      labels = c("Cell Ranger\ncalled", stringr::str_replace(stages$stage, "_", "_\n")),
      expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::scale_y_continuous(labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(title = "Nuclei retained through QC", x = NULL, y = "Nuclei",
      fill = NULL, caption = stringr::str_wrap(paste(
        "Counts summed across all configured GEM wells. Ribbon thickness represents nuclei.",
        "Overlapping exclusions are counted once, in listed action order (top to bottom at each stage).",
        "Numbers below show nuclei continuing; zero-loss actions remain in the table but have no discard branch."
      ), width = 85)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(), legend.position = "top")
}

save_QC_cell_retention_plot <- function(retention_tibble) {
  plot <- plot_QC_cell_retention(retention_tibble)
  branches_per_stage <- plot$data |>
    dplyr::filter(.data$excluded_cells > 0) |>
    dplyr::count(.data$stage)
  save_plots_structured(plot,
    width = max(8, 4 * dplyr::n_distinct(retention_tibble$stage) + 2),
    height = max(6, 3 + max(c(0, branches_per_stage$n)))
  )
}

plot_nuclei_per_donor_id <- function(
  metadata_tibble,
  fill_by = "PCA_harmony_SNN_cluster"
) {
  metadata_tibble %>%
    ggplot2::ggplot(ggplot2::aes(x = donor_id, fill = .data[[fill_by]])) +
    ggplot2::geom_bar() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust = 1))
}

#' Plot one per-dataset QC violin
#'
#' @param plot_tibble Per-feature plotting data with `GEM_well_ID`, `dataset`,
#'   `feature`, and `value` columns.
#' @param feature_thresholds Threshold intervals for the plotted feature.
#' @param show_dataset_legend Whether to show the legend for dataset fill colors.
#' @return A ggplot ready for saving or composition.
#' @keywords internal

plot_per_dataset_QC_violin <- function(
  plot_tibble,
  feature_thresholds,
  show_dataset_legend = FALSE
) {
  GEM_well_levels <- unique(as.character(plot_tibble$GEM_well_ID))
  plot_tibble <- plot_tibble |>
    dplyr::mutate(
      GEM_well_ID = factor(.data$GEM_well_ID, levels = GEM_well_levels),
      GEM_well_position = match(.data$GEM_well_ID, GEM_well_levels)
    )
  plot <- plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(
      x = GEM_well_position,
      y = value,
      fill = dataset,
      group = GEM_well_ID
    ))

  if (nrow(feature_thresholds) > 0) {
    if ("GEM_well_ID" %in% colnames(feature_thresholds)) {
      feature_thresholds <- feature_thresholds |>
        dplyr::filter(.data$GEM_well_ID %in% GEM_well_levels) |>
        dplyr::mutate(
          GEM_well_position = match(.data$GEM_well_ID, GEM_well_levels),
          xmin = .data$GEM_well_position - 0.45,
          xmax = .data$GEM_well_position + 0.45
        )
      plot <- plot +
        ggplot2::geom_rect(
          data = feature_thresholds,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = "grey60",
          alpha = 0.25,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_segment(
          data = feature_thresholds,
          ggplot2::aes(
            x = xmin,
            xend = xmax,
            y = threshold,
            yend = threshold
          ),
          color = "grey35",
          linetype = "dashed",
          inherit.aes = FALSE
        )
    } else {
      plot <- plot +
        ggplot2::geom_rect(
          data = feature_thresholds,
          ggplot2::aes(ymin = ymin, ymax = ymax),
          xmin = -Inf,
          xmax = Inf,
          fill = "grey60",
          alpha = 0.25,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_hline(
          data = feature_thresholds,
          ggplot2::aes(yintercept = threshold),
          color = "grey35",
          linetype = "dashed",
          inherit.aes = FALSE
        )
    }
  }

  plot +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(GEM_well_levels),
      labels = GEM_well_levels
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_blank(),
      legend.position = if (show_dataset_legend) "bottom" else "none"
    )
}

#' Plot per-dataset QC violins
#'
#' @param metadata_tibble Per-cell metadata containing `GEM_well_ID`,
#'   `dataset`, and requested QC feature columns.
#' @param QC_exclude_vector Character vector of QC exclusion expressions.
#' @param feature_names Character vector of QC features to plot.
#' @param GEM_well_QC_exclude_list Optional named list of GEM-well-specific QC
#'   exclusion expression vectors. Cannot be supplied together with
#'   `QC_exclude_vector`.
#' @param show_dataset_legend Whether to show the legend for dataset fill colors.
#' @return Named list with one violin plot per available QC feature.
#' @keywords internal

plot_per_dataset_QC_violins <- function(
  metadata_tibble,
  QC_exclude_vector = NULL,
  feature_names,
  GEM_well_QC_exclude_list = NULL,
  show_dataset_legend = FALSE
) {
  if (!is.null(QC_exclude_vector) && !is.null(GEM_well_QC_exclude_list)) {
    stop(
      "Supply either QC_exclude_vector or GEM_well_QC_exclude_list, not both.",
      call. = FALSE
    )
  }
  threshold_tibble <- if (is.null(GEM_well_QC_exclude_list)) {
    get_QC_exclude_threshold_tibble(
      QC_exclude_vector = QC_exclude_vector,
      feature_names = feature_names
    )
  } else {
    get_GEM_well_QC_exclude_threshold_tibble(
      GEM_well_QC_exclude_list = GEM_well_QC_exclude_list,
      feature_names = feature_names
    )
  }
  features <- intersect(feature_names, colnames(metadata_tibble))

  features |>
    purrr::set_names() |>
    purrr::map(\(feature) {
      plot_tibble <- metadata_tibble |>
        dplyr::transmute(
          GEM_well_ID,
          dataset,
          feature = .env$feature,
          value = unname(.data[[feature]])
        ) |>
        dplyr::filter(
          .by = GEM_well_ID,
          value >= stats::quantile(value, probs = 0.02, na.rm = TRUE),
          value <= stats::quantile(value, probs = 0.98, na.rm = TRUE)
        )

      plot_per_dataset_QC_violin(
        plot_tibble = plot_tibble,
        feature_thresholds = dplyr::filter(
          threshold_tibble,
          .data$feature == .env$feature
        ),
        show_dataset_legend = show_dataset_legend
      )
    })
}

#' Compare manifest QC metrics across wells, clusters or cell types
#'
#' Return one annotated plot per metric. Optional GEM-well exclusion regions
#' remain associated with their individual wells.
#'
#' @param metadata_tibble Nucleus-level metadata containing grouping, fill and
#'   QC metric columns.
#' @param QC_metric_manifest_tibble Informational QC metric manifest with metric
#'   labels, checkpoint availability, plotting quantiles, and global plotting status.
#' @param checkpoints Checkpoint names whose metrics should be plotted.
#' @param group_col Metadata column defining the violin groups.
#' @param fill_col Metadata column defining violin colors.
#' @param QC_exclude_per_GEM_well_list Named list mapping GEM-well IDs to
#'   character vectors of QC exclusion expressions.
#' @return Named list containing one annotated ggplot per available QC metric.
#' @keywords internal

plot_QC_metric_violins <- function(
  metadata_tibble,
  QC_metric_manifest_tibble,
  checkpoints,
  group_col,
  fill_col = group_col,
  QC_exclude_per_GEM_well_list = NULL
) {
  metric_manifest <- QC_metric_manifest_tibble |>
    dplyr::filter(
      .data$available_from_checkpoint %in% checkpoints,
      .data$do_plot,
      .data$metric_id %in% colnames(metadata_tibble)
    ) |>
    dplyr::filter(purrr::map_lgl(
      .data$metric_id,
      \(metric_id) {
        values <- metadata_tibble[[metric_id]]
        is.numeric(values) && any(is.finite(values))
      }
    ))
  if (nrow(metric_manifest) == 0L) {
    stop("No available QC metrics were found in the metadata.", call. = FALSE)
  }

  threshold_tibble <- if (!is.null(QC_exclude_per_GEM_well_list)) {
    stopifnot(group_col == "GEM_well_ID")
    get_GEM_well_QC_exclude_threshold_tibble(
      GEM_well_QC_exclude_list = QC_exclude_per_GEM_well_list,
      feature_names = metric_manifest$metric_id
    )
  }
  group_levels <- unique(as.character(metadata_tibble[[group_col]]))

  metric_manifest |>
    dplyr::select(
      metric_id,
      display_name,
      description,
      plot_min_q,
      plot_max_q
    ) |>
    purrr::pmap(\(metric_id, display_name, description, plot_min_q, plot_max_q) {
      plot_tibble <- metadata_tibble |>
        dplyr::transmute(
          group = factor(.data[[group_col]], levels = group_levels),
          group_position = match(.data[[group_col]], group_levels),
          fill = .data[[fill_col]],
          value = unname(.data[[metric_id]])
        ) |>
        dplyr::filter(is.finite(.data$value)) |>
        dplyr::mutate(
          .by = group,
          plot_min = if (is.na(plot_min_q)) {
            -Inf
          } else {
            stats::quantile(.data$value, probs = plot_min_q)
          },
          plot_max = if (is.na(plot_max_q)) {
            Inf
          } else {
            stats::quantile(.data$value, probs = plot_max_q)
          }
        ) |>
        dplyr::filter(
          .data$value >= .data$plot_min,
          .data$value <= .data$plot_max
        ) |>
        dplyr::select(-plot_min, -plot_max)
      feature_thresholds <- if (!is.null(threshold_tibble)) {
        threshold_tibble |>
          dplyr::filter(
            .data$feature == .env$metric_id,
            .data$GEM_well_ID %in% group_levels
          ) |>
          dplyr::mutate(
            group_position = match(.data$GEM_well_ID, group_levels),
            xmin = .data$group_position - 0.45,
            xmax = .data$group_position + 0.45
          )
      }

      plot <- plot_tibble |>
        ggplot2::ggplot(ggplot2::aes(
          x = group_position,
          y = value,
          fill = fill,
          group = group
        ))
      if (!is.null(feature_thresholds) && nrow(feature_thresholds) > 0L) {
        plot <- plot +
          ggplot2::geom_rect(
            data = feature_thresholds,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey60",
            alpha = 0.25,
            inherit.aes = FALSE
          ) +
          ggplot2::geom_segment(
            data = feature_thresholds,
            ggplot2::aes(
              x = xmin,
              xend = xmax,
              y = threshold,
              yend = threshold
            ),
            color = "grey35",
            linetype = "dashed",
            inherit.aes = FALSE
          )
      }

      plot +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::scale_x_continuous(
          breaks = seq_along(group_levels),
          labels = group_levels
        ) +
        ggplot2::labs(
          title = display_name,
          subtitle = stringr::str_wrap(
            stringr::str_c(metric_id, ": ", description),
            width = 120
          ),
          x = NULL,
          y = display_name,
          fill = if (fill_col == "dataset") "Dataset" else fill_col
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = if (fill_col == group_col) "none" else "bottom"
        )
    }) |>
    stats::setNames(metric_manifest$metric_id)
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

#' Plot within-cell-type cluster markers with explicit contrast interpretation.
plot_cluster_marker_volcano <- function(marker_result) {
  clusters <- marker_result$clusters
  plot <- plot_markers_volcano_simple(marker_result$markers)
  if (length(clusters) == 2L) {
    comparison <- paste0(
      "Cluster ", clusters[[1]], " - cluster ", clusters[[2]],
      ".\nPositive log2 fold change: higher expression in cluster ", clusters[[1]], "."
    )
    plot <- plot + ggplot2::facet_null()
  } else {
    comparison <- "Each cluster versus the pooled remaining clusters within this cell type."
  }
  plot + ggplot2::labs(
    title = marker_result$cell_type,
    subtitle = paste(comparison, "BH-adjusted across genes separately for each contrast.", sep = "\n"),
    x = "Log2 fold change (cluster / background)",
    y = "-log10(p-value)"
  )
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
#' @param source_cell_type_col,target_cell_type_col Cell-type annotation columns
#'   used to order both axes and outline matching blocks.
#' @return A ggplot tile heatmap where fill is each target count divided by the
#'   source-row total and tile labels show raw counts.
#' @keywords internal

plot_cluster_confusion_matrix <- function(
  metadata_tibble,
  source_col,
  target_col,
  source_label,
  target_label,
  title = NULL,
  source_cell_type_col = source_col,
  target_cell_type_col = target_col
) {
  missing_cols <- setdiff(c(source_col, target_col, source_cell_type_col, target_cell_type_col), names(metadata_tibble))
  if (length(missing_cols) > 0) {
    stop("Missing confusion-matrix metadata column(s): ", paste(missing_cols, collapse = ", "))
  }

  cluster_metadata <- metadata_tibble |>
    dplyr::transmute(
      source = as.character(.data[[source_col]]),
      target = as.character(.data[[target_col]]),
      source_cell_type = as.character(.data[[source_cell_type_col]]),
      target_cell_type = as.character(.data[[target_cell_type_col]])
    ) |>
    dplyr::filter(!is.na(source), !is.na(target), source != "", target != "")
  cell_type_order <- gtools::mixedsort(unique(c(cluster_metadata$source_cell_type, cluster_metadata$target_cell_type)))
  source_key <- cluster_metadata |>
    dplyr::distinct(source, source_cell_type) |>
    dplyr::arrange(match(source_cell_type, cell_type_order), match(source, gtools::mixedsort(source)))
  target_key <- cluster_metadata |>
    dplyr::distinct(target, target_cell_type) |>
    dplyr::arrange(match(target_cell_type, cell_type_order), match(target, gtools::mixedsort(target)))
  source_levels <- source_key$source
  target_levels <- target_key$target

  plot_data <- cluster_metadata |>
    dplyr::count(source, target, name = "n") |>
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
      ),
      count_text_size = pmin(2.5, 85 / max(length(source_levels), length(target_levels)), 7.5 / pmax(1, nchar(n_label)))
    )
  matching_blocks <- plot_data |>
    dplyr::mutate(
      source_cell_type = source_key$source_cell_type[match(source, source_key$source)],
      target_cell_type = target_key$target_cell_type[match(target, target_key$target)]
    ) |>
    dplyr::filter(!is.na(source_cell_type), source_cell_type != "", source_cell_type == target_cell_type) |>
    dplyr::summarise(
      xmin = min(as.integer(target)) - 0.5,
      xmax = max(as.integer(target)) + 0.5,
      ymin = min(as.integer(source)) - 0.5,
      ymax = max(as.integer(source)) + 0.5,
      .by = source_cell_type
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = target, y = source, fill = source_fraction)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::geom_rect(
      data = matching_blocks,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, color = "#A66F00", linewidth = 1, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_label, color = source_fraction > 0.55, size = count_text_size)
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    ggplot2::scale_size_identity() +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#2166AC",
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = title, x = target_label, y = source_label, fill = "Row fraction"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

#' Compare two modalities at cluster and cell-type resolution
#'
#' @param metadata_tibble Annotated cells retained at the checkpoint.
#' @param source_label,target_label Modalities: RNA, ATAC, or WNN.
#' @return A paired patchwork with a shared subtitle and legend.
plot_modality_confusion_matrices <- function(metadata_tibble, source_label, target_label) {
  prefixes <- c(RNA = "PCA", ATAC = "LSI", WNN = "WNN")
  source_cluster <- paste0(prefixes[source_label], "_harmony_SNN_cluster_named")
  target_cluster <- paste0(prefixes[target_label], "_harmony_SNN_cluster_named")
  source_cell_type <- paste0(prefixes[source_label], "_harmony_SNN_cluster_cell_type")
  target_cell_type <- paste0(prefixes[target_label], "_harmony_SNN_cluster_cell_type")
  cluster_plot <- plot_cluster_confusion_matrix(
    metadata_tibble, source_cluster, target_cluster,
    paste(source_label, "SNN cluster"), paste(target_label, "SNN cluster"),
    title = "SNN clusters",
    source_cell_type_col = source_cell_type, target_cell_type_col = target_cell_type
  )
  cell_type_plot <- plot_cluster_confusion_matrix(
    metadata_tibble, source_cell_type, target_cell_type,
    paste(source_label, "cell type"), paste(target_label, "cell type"),
    title = "Cell types"
  )
  patchwork::wrap_plots(
    cluster_plot, cell_type_plot, nrow = 1, guides = "collect",
    widths = c(nlevels(cluster_plot$data$target), nlevels(cell_type_plot$data$target))
  ) +
    patchwork::plot_annotation(
      title = paste(source_label, "vs", target_label),
      subtitle = "Colors: fraction of cells within each row. Labels: cell counts (k = 1,000). Outlines: matching cell types."
    ) &
    ggplot2::theme(legend.position = "bottom")
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

get_marker_cell_type_order <- function(groups, marker_set_names, group_cell_types = NULL) {
  observed <- levels(droplevels(as.factor(groups)))
  if (is.null(group_cell_types)) group_cell_types <- stats::setNames(observed, observed)
  matched <- unlist(lapply(marker_set_names, function(label)
    observed[which(group_cell_types[observed] == label)]), use.names = FALSE)
  c(matched, setdiff(observed, matched))
}

get_marker_group_cell_types <- function(metadata, group_col, cell_type_col = group_col) {
  mapping <- unique(metadata[, unique(c(group_col, cell_type_col))])
  mapping <- mapping[!is.na(mapping[[group_col]]), , drop = FALSE]
  if (anyDuplicated(mapping[[group_col]])) stop("Each plotted group must have one assigned cell type.")
  stats::setNames(as.character(mapping[[cell_type_col]]), as.character(mapping[[group_col]]))
}

# Shared presentation; callers supply matching rectangles before adding dots.
style_marker_dot_plot <- function(plot, matching_rows, group_label, score_type) {
  is_module <- score_type == "Adjusted UCell marker-set scores"
  plot +
    ggplot2::geom_rect(data = matching_rows,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = "grey95", color = "grey75", linewidth = 0.25, inherit.aes = FALSE) +
    ggplot2::labs(
      title = paste(score_type, "by", group_label), y = group_label,
      subtitle = if (is_module) paste0(
        "Red: evidence above matched background; blue: below it; white: near background. Read down columns for broad or weak evidence.\n",
        "Columns group similar cluster profiles; rows follow assigned cell types to retain the shaded staircase. Similarity does not establish biological identity.") else paste0(
        "Read across each marker set to check whether several genes support the same groups or one gene dominates.\n",
        "Look for coherent signal in matching rows, broad expression elsewhere, and markers with little detection."),
      caption = paste0(if (is_module)
        "Colour: cached mean UCell minus the matched-control 95th percentile, before doublet filtering. Dot area: percentage of cells with raw UCell > 0.\n" else
        "Colour: per-gene scaled mean signal across the plotted groups. Dot area: percentage of cells with detected signal.\n",
        "Shaded boxes mark assigned cell-type matches to marker sets; unassigned groups have no match. Labels derived from these markers are not independent validation.\n",
        if (is_module) "Column order: Euclidean distance and Ward clustering of unstandardized cluster-adjusted profiles; each cluster has equal weight. Negative scores are retained." else
        "Only available positive markers are shown. Shared genes can occur in multiple sets; compare individual genes before interpreting a set as specific.")) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title.position = "plot", plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 0))
}

plot_grouped_marker_dot_BPCells <- function(source, marker_genes_list, groups, group_col,
                                           group_cell_types = NULL, group_label = group_col,
                                           score_type = "Marker expression") {
  group_order <- get_marker_cell_type_order(groups, names(marker_genes_list), group_cell_types)
  if (is.null(group_cell_types)) group_cell_types <- stats::setNames(group_order, group_order)
  cell_type_order <- unique(unname(group_cell_types[group_order]))
  marker_set_order <- c(intersect(cell_type_order, names(marker_genes_list)), setdiff(names(marker_genes_list), cell_type_order))
  markers <- tibble::enframe(marker_genes_list[marker_set_order], name = "marker_set", value = "feature") |>
    tidyr::unnest_longer(feature) |>
    dplyr::filter(!stringr::str_ends(.data$feature, "-")) |>
    dplyr::mutate(feature = stringr::str_remove(.data$feature, "[+]$")) |>
    dplyr::filter(.data$feature %in% rownames(source)) |>
    dplyr::distinct(.data$marker_set, .data$feature) |>
    dplyr::mutate(
      marker_set = factor(.data$marker_set, levels = marker_set_order),
      marker_feature = factor(dplyr::row_number())
    )
  if (nrow(markers) == 0L) {
    stop("No positive marker features were found in the feature matrix.")
  }

  plot_data <- BPCells::plot_dot(
    source = source,
    features = unique(markers$feature),
    groups = groups,
    group_order = group_order,
    gene_mapping = NULL,
    return_data = TRUE
  ) |>
    dplyr::inner_join(markers, by = "feature", relationship = "many-to-many") |>
    dplyr::mutate(group = factor(.data$group, levels = rev(group_order)))
  matching_rows <- plot_data |>
    dplyr::filter(unname(group_cell_types[as.character(.data$group)]) == as.character(.data$marker_set)) |>
    dplyr::distinct(.data$marker_set, .data$group) |>
    dplyr::mutate(xmin = -Inf, xmax = Inf,
      ymin = as.integer(.data$group) - 0.45, ymax = as.integer(.data$group) + 0.45)

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$marker_feature, y = .data$group, color = .data$average, size = .data$percent))
  style_marker_dot_plot(plot, matching_rows, group_label, score_type) +
    ggplot2::geom_point() +
    ggplot2::scale_size_area() +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::scale_x_discrete(labels = stats::setNames(markers$feature, markers$marker_feature)) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(marker_set),
      switch = "x",
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(
        marker_set = \(labels) stringr::str_wrap(stringr::str_replace_all(labels, "_", " "), width = 14)
      )
    ) +
    ggplot2::labs(
      x = "Positive markers by cell-type set", size = "% Detected", color = "Mean Z-score"
    ) +
    ggplot2::theme(strip.placement = "outside")
}

#' Plot marker expression dot BPCells
#'
#' Plot BPCells marker-expression dot plots from GEX metadata and marker sets.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param marker_genes_list Named cell-type marker sets; trailing `-` entries are excluded, while trailing `+` and unsigned entries are positive markers.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @param scale_factor Scale factor used when normalizing counts, coverage, or marker scores.
#' @param cell_type_col Optional metadata column mapping plotted groups to assigned cell types for marker-set highlights.
#' @param group_label Human-readable grouping label for the title and y axis.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_marker_expression_dot_BPCells <- function(feature_matrix, metadata_tibble, marker_genes_list, group_col,
                                               scale_factor = 10000, cell_type_col = NULL, group_label = group_col) {
  metadata <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", group_col, cell_type_col)), dplyr::any_of("nCount_RNA")) |>
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

  plot_grouped_marker_dot_BPCells(
    source = log_norm_matrix,
    marker_genes_list = marker_genes_list,
    groups = metadata[[group_col]],
    group_col = group_col, group_label = group_label,
    group_cell_types = get_marker_group_cell_types(metadata, group_col,
      if (is.null(cell_type_col)) group_col else cell_type_col)
  )
}

#' Plot marker gene activity dot BPCells
#'
#' Plot gene-activity marker dot plots for ATAC-derived activity matrices.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param marker_genes_list Named cell-type marker sets; trailing `-` entries are excluded, while trailing `+` and unsigned entries are positive markers.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_marker_gene_activity_dot_BPCells <- function(feature_matrix, metadata_tibble, marker_genes_list, group_col) {
  metadata <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", group_col))) |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% colnames(feature_matrix), !is.na(.data[[group_col]])) |>
    dplyr::arrange(match(.data$barcode_w_prefix, colnames(feature_matrix)))

  plot_grouped_marker_dot_BPCells(
    source = BPCells::log1p_slow(feature_matrix[, metadata$barcode_w_prefix, drop = FALSE]),
    marker_genes_list = marker_genes_list,
    groups = metadata[[group_col]],
    group_col = group_col, score_type = "Marker gene activity"
  )
}

#' Plot cached annotation evidence, with one distance order shared by both views.
#' Cell-type rows average cluster-adjusted scores by cell count, not a new null.
plot_UCell_annotation_dot <- function(annotation, metadata_tibble, group_by = c("cluster", "cell_type"),
                                       cluster_column = "PCA_harmony_SNN_cluster") {
  group_by <- match.arg(group_by)
  evidence <- annotation$evidence
  profiles <- stats::xtabs(excess ~ cluster + label, evidence)
  ordering <- stats::hclust(stats::dist(t(profiles)), method = "ward.D2")
  module_names <- colnames(profiles)[ordering$order]
  cell_type_col <- paste0(cluster_column, "_cell_type")
  group_col <- if (group_by == "cluster") paste0(cluster_column, "_named") else cell_type_col
  group_label <- if (group_by == "cluster") "GEX cluster" else "GEX cell type"
  group_cell_types <- get_marker_group_cell_types(metadata_tibble, group_col, cell_type_col)
  group_order <- get_marker_cell_type_order(metadata_tibble[[group_col]], module_names, group_cell_types)
  cluster_groups <- unique(data.frame(cluster = as.character(metadata_tibble[[cluster_column]]),
    group = as.character(metadata_tibble[[group_col]])))
  stopifnot(!anyDuplicated(cluster_groups$cluster),
    setequal(cluster_groups$cluster, evidence$cluster),
    setequal(metadata_tibble$barcode_w_prefix, rownames(annotation$cell_scores)))
  positive <- dplyr::bind_rows(lapply(module_names, function(label) {
    percent <- tapply(annotation$cell_scores[metadata_tibble$barcode_w_prefix, label] > 0,
      as.character(metadata_tibble[[group_col]]), mean) * 100
    data.frame(group = names(percent), label, pct_positive = unname(percent))
  }))
  plot_tibble <- dplyr::left_join(evidence, cluster_groups, by = "cluster")
  plot_tibble <- if (group_by == "cell_type") plot_tibble |>
    dplyr::summarise(adjusted_score = stats::weighted.mean(.data$excess, .data$cells),
      .by = c("group", "label")) else plot_tibble |>
    dplyr::mutate(adjusted_score = .data$excess)
  plot_tibble <- plot_tibble |>
    dplyr::left_join(positive, by = c("group", "label")) |>
    dplyr::mutate(module = factor(.data$label, levels = module_names),
      cluster = factor(.data$group, levels = rev(group_order)))

  matching_rows <- plot_tibble |>
    dplyr::filter(unname(group_cell_types[as.character(.data$cluster)]) == as.character(.data$module)) |>
    dplyr::mutate(xmin = as.integer(.data$module) - 0.45, xmax = as.integer(.data$module) + 0.45,
      ymin = as.integer(.data$cluster) - 0.45, ymax = as.integer(.data$cluster) + 0.45)
  plot <- ggplot2::ggplot(plot_tibble, ggplot2::aes(x = .data$module, y = .data$cluster, color = .data$adjusted_score, size = .data$pct_positive))
  plot <- style_marker_dot_plot(plot, matching_rows, group_label, "Adjusted UCell marker-set scores")
  plot +
    ggplot2::geom_point() +
    ggplot2::scale_size_area(limits = c(0, 100), max_size = 7) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-1, 1) * max(abs(evidence$excess))) +
    ggplot2::labs(x = "Marker set", color = "Adjusted UCell score", size = "% raw score > 0",
      caption = paste0(plot$labels$caption, "\n", if (group_by == "cell_type")
        "Cell-type rows are cell-count-weighted means of the cached cluster-adjusted scores; no pooled background is re-estimated." else
        "Cluster colours use exactly the evidence used for annotation. Positive evidence still needs sufficient separation from competing labels."))
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
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
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

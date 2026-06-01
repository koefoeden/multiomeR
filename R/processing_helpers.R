get_excluded_BCs <- function(metadata_tibble, QC_exclude_vector) {
  if (is.null(QC_exclude_vector)) {
    return(list())
  }
  QC_exclude_vector %>%
    purrr::set_names() %>%
    purrr::map(.f = \(cut_off_expression) {
      # turn character vector of QC-cut-off into an evaluatable statement, used in filter step below.
      filter_expression <- rlang::parse_expr(cut_off_expression)
      metadata_tibble %>%
        dplyr::filter(!!filter_expression) %>%
        dplyr::pull(barcode_w_prefix)
    })
}

#' Get QC exclude threshold tibble
#'
#' Parse simple QC cutoff expressions into threshold intervals for plotting.
#'
#' @param QC_exclude_vector Character vector of dplyr filter expressions such as `nCount_RNA < 1000`; each expression defines one exclusion set.
#' @param feature_names Character vector of valid feature names; threshold expressions referencing other names are ignored.
#' @return A tibble with feature, threshold, and y-range columns for expressions
#'   of the form `<feature> <op> <number>`. Unsupported expressions are ignored.
#' @keywords internal

get_QC_exclude_threshold_tibble <- function(QC_exclude_vector, feature_names) {
  empty_thresholds <- tibble::tibble(
    feature = character(),
    threshold = double(),
    ymin = double(),
    ymax = double()
  )

  if (is.null(QC_exclude_vector)) {
    return(empty_thresholds)
  }

  purrr::map_dfr(QC_exclude_vector, \(cut_off_expression) {
    filter_expression <- rlang::parse_expr(cut_off_expression)
    if (!rlang::is_call(filter_expression)) {
      return(empty_thresholds)
    }

    operator <- rlang::call_name(filter_expression)
    expression_args <- rlang::call_args(filter_expression)
    if (
      !operator %in% c("<", "<=", ">", ">=") ||
        length(expression_args) != 2 ||
        !rlang::is_symbol(expression_args[[1]]) ||
        !is.numeric(expression_args[[2]])
    ) {
      return(empty_thresholds)
    }

    feature <- rlang::as_name(expression_args[[1]])
    threshold <- as.numeric(expression_args[[2]])

    if (!feature %in% feature_names || length(threshold) != 1 || !is.finite(threshold)) {
      return(empty_thresholds)
    }

    tibble::tibble(
      feature = feature,
      threshold = threshold,
      ymin = if (operator %in% c("<", "<=")) -Inf else threshold,
      ymax = if (operator %in% c("<", "<=")) threshold else Inf
    )
  })
}

#' Plot upset from excluded BCs list
#'
#' Plot overlaps among QC-excluded barcode sets as an UpSet-style chart.
#'
#' @param QC_excluded_BCs_list Named list of barcode vectors, one element per QC exclusion reason.
#' @param n_total Denominator used when reporting exclusion-set and intersection percentages.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_upset_from_excluded_BCs_list <- function(QC_excluded_BCs_list, n_total) {
  format_QC_percent <- function(percent) {
    dplyr::case_when(
      percent >= 10 ~ stringr::str_c(round(percent), "%"),
      percent >= 1 ~ stringr::str_c(round(percent, 1), "%"),
      .default = stringr::str_c(signif(percent, 2), "%")
    )
  }

  QC_excluded_BCs_list <- purrr::discard(QC_excluded_BCs_list, \(barcodes) length(barcodes) == 0)
  set_names <- names(QC_excluded_BCs_list)

  if (length(set_names) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No excluded barcodes") +
        ggplot2::theme_void() +
        ggplot2::labs(title = stringr::str_glue("Excluded barcodes by QC filter type (input barcodes: {n_total})"))
    )
  }
  if (length(set_names) > 30) {
    stop("At most 30 QC exclusion categories are supported in the compact UpSet plot.")
  }

  set_bits <- bitwShiftL(1L, seq_along(set_names) - 1L)
  names(set_bits) <- set_names
  membership_dt <- data.table::rbindlist(
    purrr::map2(
      .x = QC_excluded_BCs_list,
      .y = set_bits,
      ~ data.table::data.table(barcode = unique(.x), bit = .y)
    ),
    use.names = TRUE
  )

  barcode_masks <- membership_dt[, .(mask = sum(bit)), by = "barcode"]
  intersection_counts <- barcode_masks[, .(n = .N), by = "mask"][order(-n, mask)]
  intersection_counts[, intersection_idx := .I]
  rm(membership_dt, barcode_masks)

  set_size_tibble <- tibble::tibble(
    set = set_names,
    bit = unname(set_bits),
    set_idx = seq_along(set_names)
  ) |>
    dplyr::mutate(
      y = length(set_names) - .data$set_idx + 1,
      n = purrr::map_int(.data$bit, \(bit) sum(intersection_counts$n[bitwAnd(intersection_counts$mask, bit) > 0])),
      percent = 100 * .data$n / n_total
    )

  matrix_tibble <- tidyr::crossing(
    intersection_idx = intersection_counts$intersection_idx,
    set = set_names
  ) |>
    dplyr::left_join(dplyr::select(set_size_tibble, set, bit, y), by = "set") |>
    dplyr::left_join(dplyr::select(intersection_counts, intersection_idx, mask), by = "intersection_idx") |>
    dplyr::mutate(present = bitwAnd(.data$mask, .data$bit) > 0)
  matrix_segments <- matrix_tibble |>
    dplyr::filter(.data$present) |>
    dplyr::summarise(ymin = min(.data$y), ymax = max(.data$y), .by = "intersection_idx") |>
    dplyr::filter(.data$ymin != .data$ymax)

  intersection_tibble <- intersection_counts |>
    tibble::as_tibble() |>
    dplyr::mutate(
      percent = 100 * .data$n / n_total,
      percent_label = format_QC_percent(.data$percent)
    )

  bar_plot <- intersection_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$intersection_idx, y = .data$percent)) +
    ggplot2::geom_col(width = 0.65, fill = "gray23") +
    ggplot2::geom_text(ggplot2::aes(label = .data$percent_label), angle = 45, hjust = 0, vjust = -0.5, size = 2.2) +
    ggplot2::scale_x_continuous(limits = c(0.5, nrow(intersection_tibble) + 0.5), expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(labels = format_QC_percent, expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "Intersection size (% of input)") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  matrix_plot <- matrix_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$intersection_idx, y = .data$y)) +
    ggplot2::geom_point(color = "grey85", size = 2.2) +
    ggplot2::geom_segment(
      data = matrix_segments,
      ggplot2::aes(xend = .data$intersection_idx, y = .data$ymin, yend = .data$ymax),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(data = dplyr::filter(matrix_tibble, .data$present), color = "gray23", size = 2.2) +
    ggplot2::scale_x_continuous(limits = c(0.5, nrow(intersection_tibble) + 0.5), expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = set_size_tibble$y, labels = set_size_tibble$set, expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  set_size_plot <- set_size_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$percent, y = .data$y)) +
    ggplot2::geom_col(width = 0.65, fill = "gray23", orientation = "y") +
    ggplot2::scale_x_reverse(labels = format_QC_percent, expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::scale_y_continuous(breaks = set_size_tibble$y, labels = NULL, expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(x = "Set size (% of input)", y = NULL) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  (patchwork::plot_spacer() + bar_plot + patchwork::plot_layout(widths = c(0.35, 1))) /
    (set_size_plot + matrix_plot + patchwork::plot_layout(widths = c(0.35, 1))) +
    patchwork::plot_layout(heights = c(0.65, 0.35)) +
    patchwork::plot_annotation(title = stringr::str_glue("Excluded barcodes by QC filter type (input barcodes: {n_total})"))
}

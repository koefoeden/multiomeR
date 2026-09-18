#' Measurement-support presets for peak-gene hypotheses
peak_gene_filter_settings <- function() {
  tibble::tibble(
    filter = c("lenient", "moderate", "strict"),
    RNA_count = c(5, 10, 10), ATAC_count = c(3, 5, 5),
    min_aggregates = c(6L, 6L, 10L), aggregate_fraction = c(0.1, 0.1, 0.2),
    min_donors = c(2L, 2L, 3L), aggregates_per_donor = c(2L, 2L, 3L)
  )
}

#' Cache selected candidates and retention diagnostics for all three presets
#'
#' Support is determined separately for RNA and ATAC. Supported donors must
#' overlap, but qualifying aggregates need not coincide. All observations remain
#' in each retained regression. No association statistic enters the filter.
filter_peak_gene_candidate_pairs <- function(
  aggregate_matrices, normalized_aggregate_matrices, candidate_pairs_tibble, filter
) {
  settings <- peak_gene_filter_settings()
  if (length(filter) != 1L || is.na(filter) || !filter %in% settings$filter) {
    stop("peak_gene_correlation_filter must be lenient, moderate, or strict")
  }
  branch <- prepare_peak_gene_correlation_branch(normalized_aggregate_matrices, candidate_pairs_tibble)
  pairs <- branch$candidate_pairs
  if (!is.na(branch$skipped_reason)) pairs <- pairs[0, ]
  diagnostics <- settings |>
    dplyr::mutate(cell_group = branch$cell_group, chr = branch$chr,
      selected = .data$filter == !!filter, n_candidates = nrow(pairs),
      n_retained = 0L, n_low_RNA_support = 0L, n_low_ATAC_support = 0L,
      n_low_donor_support = 0L, skipped_reason = branch$skipped_reason)
  if (!nrow(pairs)) return(list(candidate_pairs = pairs, diagnostics = diagnostics))
  depth <- aggregate_matrices$aggregate_depth_tibble
  RNA <- as.matrix(aggregate_matrices$GEX_counts)
  ATAC <- as.matrix(aggregate_matrices$ATAC_counts)
  stopifnot(identical(colnames(RNA), depth$aggregate_id),
    identical(colnames(ATAC), depth$aggregate_id),
    all(is.finite(depth$GEX_depth) & depth$GEX_depth > 0),
    all(is.finite(depth$ATAC_depth) & depth$ATAC_depth > 0))
  gi <- match(pairs$gene_matrix_feature, rownames(RNA))
  pi <- match(pairs$peak, rownames(ATAC))
  stopifnot(!anyNA(gi), !anyNA(pi))
  donor_indices <- split(seq_len(nrow(depth)), depth$donor_id)
  for (k in seq_len(nrow(settings))) {
    cfg <- settings[k, ]
    gene_support <- sweep(RNA, 2,
      pmax(2, cfg$RNA_count * depth$GEX_depth / stats::median(depth$GEX_depth)), ">=")
    peak_support <- sweep(ATAC, 2,
      pmax(2, cfg$ATAC_count * depth$ATAC_depth / stats::median(depth$ATAC_depth)), ">=")
    minimum <- max(cfg$min_aggregates, ceiling(cfg$aggregate_fraction * nrow(depth)))
    low_gene <- rowSums(gene_support)[gi] < minimum
    low_peak <- rowSums(peak_support)[pi] < minimum
    common_donors <- integer(nrow(pairs))
    for (indices in donor_indices) {
      gene_d <- rowSums(gene_support[, indices, drop = FALSE]) >= cfg$aggregates_per_donor
      peak_d <- rowSums(peak_support[, indices, drop = FALSE]) >= cfg$aggregates_per_donor
      common_donors <- common_donors + as.integer(gene_d[gi] & peak_d[pi])
    }
    low_donor <- common_donors < cfg$min_donors
    keep <- !low_gene & !low_peak & !low_donor
    diagnostics$n_retained[[k]] <- sum(keep)
    diagnostics$n_low_RNA_support[[k]] <- sum(low_gene)
    diagnostics$n_low_ATAC_support[[k]] <- sum(low_peak)
    diagnostics$n_low_donor_support[[k]] <- sum(low_donor)
    if (cfg$filter == filter) selected_pairs <- pairs[keep, ]
  }
  list(candidate_pairs = selected_pairs, diagnostics = diagnostics)
}

#' Compare hypothesis retention under each support preset
plot_peak_gene_filter_retention <- function(diagnostics, filter) {
  plot_data <- diagnostics |>
    dplyr::summarise(n_candidates = sum(.data$n_candidates),
      n_retained = sum(.data$n_retained), .by = c("cell_group", "filter")) |>
    dplyr::filter(.data$n_candidates > 0) |>
    dplyr::mutate(retained_fraction = .data$n_retained / .data$n_candidates,
      filter = factor(.data$filter, levels = peak_gene_filter_settings()$filter))
  if (!nrow(plot_data)) return(make_empty_peak_gene_correlation_plot("No eligible peak-gene pairs before support filtering"))
  settings <- peak_gene_filter_settings()
  methods <- vapply(seq_len(nrow(settings)), function(i) {
    x <- settings[i, ]
    paste0(x$filter, ": RNA ", x$RNA_count, "/ATAC ", x$ATAC_count,
      "; ≥max(", x$min_aggregates, ", ", 100 * x$aggregate_fraction,
      "% of aggregates); ≥", x$min_donors, " shared donors with ≥",
      x$aggregates_per_donor, " supported aggregates per feature/donor")
  }, character(1))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$retained_fraction,
    y = .data$cell_group, fill = .data$filter)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(title = "Peak-gene pairs retained by measurement-support filters",
      subtitle = stringr::str_wrap(paste0("Configured filter: ", filter,
        ". Compare count and donor coverage across cell types; retention does not establish an association, and few donors can limit eligibility."), 110),
      caption = stringr::str_wrap(paste0("Denominator: pairs passing distance, detection and aggregate eligibility before support filtering. Counts scale with aggregate library depth relative to its cell-type median, with a two-count raw floor. ",
        paste(methods, collapse = ". "),
        ". Both features must qualify; their supported aggregates need not coincide. All observations remain in retained models. Failed fits within the retained family still count toward BH correction."), 140),
      x = "Fraction of eligible pairs retained", y = NULL, fill = "Filter") +
    ggplot2::theme(legend.position = "top")
}

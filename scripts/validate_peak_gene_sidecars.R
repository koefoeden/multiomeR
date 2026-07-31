#!/usr/bin/env Rscript

source("R/peak_gene_correlation_helpers.R")

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

expect_identical <- function(observed, expected, label) {
  if (!identical(observed, expected)) {
    fail(label, ": observed result differs from the legacy result")
  }
  invisible(TRUE)
}

normalized_branch <- list(
  cell_group = "B cell",
  chr = "chr1",
  GEX_norm = matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("ENSG1", "ENSG2"), c("agg1", "agg2", "agg3"))
  ),
  ATAC_norm = matrix(
    c(7, 8, 9, 10, 11, 12),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("chr1:100-200", "chr1:300-400"), c("agg1", "agg2", "agg3"))
  ),
  aggregate_depth_tibble = tibble::tibble(
    aggregate_id = c("agg1", "agg2", "agg3"),
    n_cells = c(10, 11, 12),
    GEX_depth = c(1000, 1100, 1200),
    ATAC_depth = c(2000, 2100, 2200)
  )
)

top_links <- tibble::tibble(
  scatter_plot_name = c("link_1", "link_other_chr"),
  cell_group = c("T cell", "T cell"),
  chr = c("chr1", "chr2"),
  gene_matrix_feature = c("ENSG1", "ENSG1"),
  peak = c("chr1:100-200", "chr2:100-200"),
  TargetGeneID = c("ENSG1", "ENSG1"),
  TargetGene = c("GENE1", "GENE1"),
  correlation = c(0.8, 0.7),
  FDR = c(0.01, 0.02),
  rank_in_cell_group = c(1L, 2L),
  rank_for_gene = c(1L, 1L)
)

observed_values <- extract_peak_gene_correlation_top_link_aggregate_values(
  normalized_aggregate_matrices = normalized_branch,
  top_links_tibble = top_links
)

expected_values <- dplyr::bind_cols(
  tibble::tibble(
    scatter_plot_name = rep("link_1", 3),
    primary_cell_group = rep("T cell", 3),
    cell_group = rep("B cell", 3),
    chr = rep("chr1", 3),
    peak = rep("chr1:100-200", 3),
    TargetGeneID = rep("ENSG1", 3),
    TargetGene = rep("GENE1", 3),
    correlation = rep(0.8, 3),
    FDR = rep(0.01, 3),
    rank_in_cell_group = rep(1L, 3),
    rank_for_gene = rep(1L, 3)
  ),
  normalized_branch$aggregate_depth_tibble,
  tibble::tibble(
    gene_expression_logCPM = c(1, 2, 3),
    peak_accessibility_logCPM = c(7, 8, 9)
  )
)

expect_identical(observed_values, expected_values, "top-link aggregate extraction")

empty_values <- extract_peak_gene_correlation_top_link_aggregate_values(
  normalized_aggregate_matrices = normalized_branch,
  top_links_tibble = dplyr::filter(top_links, .data$chr == "chr2")
)
if (nrow(empty_values) != 0L || !identical(names(empty_values), names(expected_values))) {
  fail("empty top-link extraction does not preserve the compact output schema")
}

results <- tibble::tibble(
  cell_group = c("B cell", "B cell", "B cell", "T cell", "T cell", "T cell"),
  correlation = c(-0.2, 0.12, NA, 0.33, 0.34, -0.1),
  FDR = c(0.2, 0.01, 0.03, 0.04, 0.06, 0.01),
  isSelfPromoter = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
  distance = c(-5100, 12000, 100, 6000, 9000, -100)
)

bin_width <- 0.025
legacy_histogram <- results |>
  dplyr::filter(!is.na(.data$correlation)) |>
  dplyr::mutate(
    correlation_bin = pmin(
      1 - bin_width,
      pmax(-1, floor((.data$correlation + 1) / bin_width) * bin_width - 1)
    ),
    correlation_mid = .data$correlation_bin + bin_width / 2
  ) |>
  dplyr::count(.data$cell_group, .data$correlation_mid, name = "n_pairs")
histogram <- summarize_peak_gene_correlation_histogram(results, bin_width)
expect_identical(histogram, legacy_histogram, "histogram summary")

legacy_support <- results |>
  dplyr::summarise(
    tested_pairs = dplyr::n(),
    FDR_significant_pairs = sum(.data$FDR < 0.05, na.rm = TRUE),
    positive_non_promoter_links = sum(
      .data$correlation > 0 & .data$FDR < 0.05 & !.data$isSelfPromoter,
      na.rm = TRUE
    ),
    .by = "cell_group"
  ) |>
  tidyr::pivot_longer(-"cell_group", names_to = "metric", values_to = "n") |>
  dplyr::mutate(n_for_plot = pmax(.data$n, 1))
support <- summarize_peak_gene_correlation_support_counts(results)
expect_identical(support, legacy_support, "support-count summary")

legacy_distance <- results |>
  dplyr::filter(!is.na(.data$correlation), !.data$isSelfPromoter) |>
  dplyr::mutate(
    abs_distance_bin = pmin(250000, floor(abs(.data$distance) / 5000) * 5000)
  ) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    median_correlation = stats::median(.data$correlation, na.rm = TRUE),
    significant_fraction = mean(.data$FDR < 0.05, na.rm = TRUE),
    .by = c("cell_group", "abs_distance_bin")
  )
distance <- summarize_peak_gene_correlation_by_distance(results)
expect_identical(distance, legacy_distance, "distance summary")

plots <- list(
  plot_peak_gene_correlation_histogram(histogram, bin_width),
  plot_peak_gene_correlation_support_counts(support),
  plot_peak_gene_correlation_by_distance(distance)
)
plot_data <- list(histogram, support, distance)
for (index in seq_along(plots)) {
  if (!inherits(plots[[index]], "ggplot")) {
    fail("plot helper ", index, " did not return a ggplot")
  }
  expect_identical(plots[[index]]$data, plot_data[[index]], paste0("plot data ", index))
}

cat("peak-gene sidecar validation ok\n")

source_project_file("R/peak_gene_correlation_helpers.R")

make_normalized_branch <- function() {
  list(
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
      donor_id = c("donor1", "donor2", "donor3"),
      state_bin = c("state1", "state1", "state2"),
      n_cells = c(10, 11, 12),
      GEX_depth = c(1000, 1100, 1200),
      ATAC_depth = c(2000, 2100, 2200)
    )
  )
}

make_top_links <- function() {
  tibble::tibble(
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
}

testthat::test_that("top-link aggregate values match the legacy result", {
  normalized_branch <- make_normalized_branch()
  observed_values <- extract_peak_gene_correlation_top_link_aggregate_values(
    normalized_aggregate_matrices = normalized_branch,
    top_links_tibble = make_top_links()
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

  testthat::expect_identical(observed_values, expected_values)
})

testthat::test_that("empty top-link extraction preserves the compact schema", {
  normalized_branch <- make_normalized_branch()
  empty_values <- extract_peak_gene_correlation_top_link_aggregate_values(
    normalized_aggregate_matrices = normalized_branch,
    top_links_tibble = dplyr::filter(make_top_links(), .data$chr == "chr2")
  )
  expected_names <- names(extract_peak_gene_correlation_top_link_aggregate_values(
    normalized_aggregate_matrices = normalized_branch,
    top_links_tibble = dplyr::filter(make_top_links(), .data$chr == "chr1")
  ))

  testthat::expect_equal(nrow(empty_values), 0L)
  testthat::expect_identical(names(empty_values), expected_names)
})

make_peak_gene_results <- function() {
  tibble::tibble(
    cell_group = c("B cell", "B cell", "B cell", "T cell", "T cell", "T cell"),
    correlation = c(-0.2, 0.12, NA, 0.33, 0.34, -0.1),
    FDR = c(0.2, 0.01, 0.03, 0.04, 0.06, 0.01),
    isSelfPromoter = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
    isTargetGeneBody = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
    distance = c(-5100, 12000, 100, 6000, 9000, -100)
  )
}

bin_width <- 0.025
testthat::test_that("the histogram summary matches the legacy result", {
  histogram <- summarize_peak_gene_correlation_histogram(
    make_peak_gene_results(),
    bin_width
  )
  expected <- tibble::tibble(
    cell_group = c("B cell", "B cell", "T cell", "T cell"),
    correlation_mid = c(-0.1875, 0.1125, -0.0875, 0.3375),
    n_pairs = c(1L, 1L, 1L, 2L)
  )

  testthat::expect_equal(histogram, expected)
})

testthat::test_that("support counts match the legacy result", {
  support <- summarize_peak_gene_correlation_support_counts(make_peak_gene_results())
  expected <- tibble::tibble(
    cell_group = rep(c("B cell", "T cell"), each = 3),
    metric = rep(
      c("tested_pairs", "FDR_significant_pairs", "candidate_enhancer_links"),
      2
    ),
    n = c(3L, 2L, 0L, 3L, 2L, 1L),
    n_for_plot = c(3, 2, 1, 3, 2, 1)
  )

  testthat::expect_equal(support, expected)
})

testthat::test_that("distance summaries match the legacy result", {
  distance <- summarize_peak_gene_correlation_by_distance(make_peak_gene_results())
  expected <- tibble::tibble(
    cell_group = c("B cell", "B cell", "T cell"),
    abs_distance_bin = c(5000, 10000, 5000),
    n_pairs = c(1L, 1L, 2L),
    median_correlation = c(-0.2, 0.12, 0.335),
    significant_fraction = c(0, 1, 0.5)
  )

  testthat::expect_equal(distance, expected)
})

testthat::test_that("sidecar plot helpers return ggplots with the supplied data", {
  results <- make_peak_gene_results()
  plot_data <- list(
    summarize_peak_gene_correlation_histogram(results, bin_width),
    summarize_peak_gene_correlation_support_counts(results),
    summarize_peak_gene_correlation_by_distance(results)
  )
  plots <- list(
    plot_peak_gene_correlation_histogram(plot_data[[1]], bin_width),
    plot_peak_gene_correlation_support_counts(plot_data[[2]]),
    plot_peak_gene_correlation_by_distance(plot_data[[3]])
  )

  for (index in seq_along(plots)) {
    testthat::expect_s3_class(plots[[index]], "ggplot")
    testthat::expect_identical(plots[[index]]$data, plot_data[[index]])
  }
})

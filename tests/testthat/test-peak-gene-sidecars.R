source_project_file("R/peak_gene_correlation_helpers.R")
source_project_file("R/peak_gene_finemapping_helpers.R")

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
    is_analyzable_link = c(TRUE, TRUE),
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
      is_analyzable_link = rep(TRUE, 3),
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
      peak_accessibility_logCPM = c(7, 8, 9),
      gene_expression_residual = c(0, 0, 0),
      peak_accessibility_residual = c(0, 0, 0)
    )
  )

  testthat::expect_equal(observed_values, expected_values, tolerance = 1e-10)
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

# Conditional inference must work with one donor and retain the full-model
# leverage correction; merely running the pipeline cannot verify these numbers.
make_peak_gene_scoring_case <- function(n_donors = 1L) {
  set.seed(17)
  n <- 60L
  donor <- rep(seq_len(n_donors), length.out = n)
  state <- rep(seq_len(20), length.out = n)
  peak <- stats::rnorm(n)
  depth <- exp(stats::rnorm(n, 8, 0.4))
  gene <- 5 + 0.6 * peak + donor + 0.3 * log(depth) +
    stats::rnorm(n, sd = 0.2 + abs(peak))
  list(
    cell_group = "example", chr = "chr1",
    GEX_norm = matrix(gene, 1, dimnames = list("gene", NULL)),
    ATAC_norm = rbind(peak = peak, constant = rep(1, n)),
    aggregate_depth_tibble = tibble::tibble(
      donor_id = as.character(donor), state_bin = as.character(state),
      GEX_depth = depth, ATAC_depth = rep(10000, n)
    ),
    gene_detected_frac = c(gene = 1), peak_accessible_frac = c(peak = 1, constant = 1),
    mean_gene_expression = c(gene = mean(gene)),
    mean_peak_accessibility = c(peak = mean(peak), constant = 1)
  )
}

testthat::test_that("one- and two-donor HC3 statistics match a full regression", {
  pairs <- tibble::tibble(chr = "chr1", gene_matrix_feature = "gene", peak = c("peak", "constant"))
  for (n_donors in c(1L, 2L)) {
    branch <- make_peak_gene_scoring_case(n_donors)
    observed <- score_peak_gene_correlations_for_cell_group(branch, pairs)
    design <- make_peak_gene_correlation_design_matrix(branch$aggregate_depth_tibble)
    gene <- as.numeric(branch$GEX_norm)
    peak <- as.numeric(branch$ATAC_norm["peak", ])
    fit <- stats::lm(gene ~ design + peak - 1)
    reference_SE <- sqrt(diag(sandwich::vcovHC(fit, type = "HC3")))[["peak"]]
    reference_p <- 2 * stats::pt(abs(stats::coef(fit)[["peak"]] / reference_SE),
      df = stats::df.residual(fit), lower.tail = FALSE)
    row <- observed[observed$peak == "peak", ]
    testthat::expect_equal(row$coefficient, unname(stats::coef(fit)[["peak"]]), tolerance = 1e-10)
    testthat::expect_equal(row$association_SE, reference_SE, tolerance = 1e-10)
    testthat::expect_equal(row$nominal_pvalue, reference_p, tolerance = 1e-10)
    testthat::expect_equal(row$n_informative_donors, n_donors)
    testthat::expect_true(is.na(observed$correlation[observed$peak == "constant"]))
    testthat::expect_true(is.na(observed$nominal_pvalue[observed$peak == "constant"]))
    testthat::expect_false(any(grepl("state", colnames(design))))
  }
})

testthat::test_that("a single donor retains non-overlapping state pseudobulks", {
  set.seed(1)
  barcodes <- paste0("cell", seq_len(1000))
  metadata <- tibble::tibble(barcode_w_prefix = barcodes, donor_id = "donor1",
    nCount_RNA = 1000, nCount_ATAC = 2000)
  embedding <- matrix(stats::rnorm(1000 * 3), 1000, dimnames = list(barcodes, paste0("LSI_", 2:4)))
  record <- make_peak_gene_correlation_donor_state_record(
    tibble::tibble(cell_group = "example", barcodes = list(barcodes)), metadata, embedding)
  testthat::expect_true(is.na(record$diagnostics$skipped_reason))
  testthat::expect_equal(record$diagnostics$n_donors, 1L)
  testthat::expect_gte(nrow(record$aggregates), 10L)
  testthat::expect_equal(anyDuplicated(unlist(record$aggregates$barcodes)), 0L)
})

testthat::test_that("compact FDR breakpoints reproduce global BH values within chromosome slices", {
  set.seed(42)
  results <- tibble::tibble(
    cell_group = rep(c("A", "B", "no_tests"), each = 1000),
    chr = rep(c("chr1", "chr2"), 1500),
    nominal_pvalue = c(
      c(0, 1, NA, rep(.001, 10), stats::runif(987)),
      c(rep(.02, 100), stats::runif(900)^4),
      rep(NA_real_, 1000)
    )
  ) |>
    dplyr::mutate(FDR = stats::p.adjust(.data$nominal_pvalue, "BH"), .by = "cell_group")
  reference <- make_peak_gene_finemapping_reference(
    results, tibble::tibble(TargetGeneID = c("g", "g"), gene_matrix_feature = c("G", "G"))
  )
  testthat::expect_lt(nrow(reference$FDR), sum(is.finite(results$nominal_pvalue)))
  testthat::expect_equal(nrow(reference$gene_features), 1L)
  branches <- results |> dplyr::group_by(.data$cell_group, .data$chr) |> dplyr::group_split()
  for (branch in branches) {
    observed <- restore_peak_gene_correlation_FDR(dplyr::select(branch, -"FDR"), reference$FDR)
    testthat::expect_identical(observed$FDR, branch$FDR)
  }
  empty <- restore_peak_gene_correlation_FDR(dplyr::select(results[0, ], -"FDR"), reference$FDR)
  testthat::expect_identical(empty$FDR, numeric())
})

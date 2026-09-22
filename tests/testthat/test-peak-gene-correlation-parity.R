source_project_file("R/peak_gene_correlation_helpers.R")
source_project_file("R/peak_gene_finemapping_helpers.R")

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

testthat::test_that("HC3 statistics match a full regression with one to six donors", {
  require_reference_version("sandwich", "3.1.1")
  pairs <- tibble::tibble(chr = "chr1", gene_matrix_feature = "gene", peak = c("peak", "constant"))
  for (n_donors in c(1L, 2L, 6L)) {
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
    testthat::expect_true(is.na(observed$correlation[observed$peak == "constant"]))
    testthat::expect_true(is.na(observed$nominal_pvalue[observed$peak == "constant"]))
  }
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
  branches <- results |> dplyr::group_by(.data$cell_group, .data$chr) |> dplyr::group_split()
  for (branch in branches) {
    observed <- restore_peak_gene_correlation_FDR(dplyr::select(branch, -"FDR"), reference$FDR)
    testthat::expect_identical(observed$FDR, branch$FDR)
  }
})

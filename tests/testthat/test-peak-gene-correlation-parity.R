source_project_file("R/peak_gene_finemapping_helpers.R")

testthat::test_that("compact FDR breakpoints reproduce global BH values within chromosome slices", {
  set.seed(42)
  results <- tibble::tibble(
    cell_group = rep(c("A", "B", "no_tests"), each = 1000),
    chr = rep(c("chr1", "chr2"), 1500),
    hierarchical_pvalue = c(
      c(0, 1, NA, rep(.001, 10), stats::runif(987)),
      c(rep(.02, 100), stats::runif(900)^4),
      rep(NA_real_, 1000)
    )
  ) |>
    dplyr::mutate(hierarchical_FDR = stats::p.adjust(.data$hierarchical_pvalue, "BH", n = dplyr::n()), .by = "cell_group")
  reference <- make_peak_gene_finemapping_reference(results)
  branches <- results |> dplyr::group_by(.data$cell_group, .data$chr) |> dplyr::group_split()
  for (branch in branches) {
    observed <- restore_peak_gene_correlation_FDR(dplyr::select(branch, -"hierarchical_FDR"), reference)
    testthat::expect_identical(observed$hierarchical_FDR, branch$hierarchical_FDR)
  }
})

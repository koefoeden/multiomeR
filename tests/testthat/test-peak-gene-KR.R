source_project_file("R/peak_gene_correlation_helpers.R")
source_project_file("R/peak_gene_hierarchical_helpers.R")
source_project_file("R/peak_gene_KR_helpers.R")

testthat::test_that("compiled batches preserve donor-slope inference and isolate bad rows", {
  kernel <- load_peak_gene_KR_kernel(file.path(multiomeR_project_root, "src/peak_gene_KR.cpp"))
  set.seed(1701)
  donor <- rep(letters[1:6], each = 8)
  peak <- stats::rnorm(length(donor))
  depth <- data.frame(donor_id = donor,
    GEX_depth = stats::runif(length(donor), 500, 2000),
    ATAC_depth = stats::runif(length(donor), 500, 2000))
  genes <- rbind(
    0.4 * peak + stats::rnorm(length(donor)),
    rep(c(-0.3, 0.1, 0.4, 0.7, 1, 1.3), each = 8) * peak + stats::rnorm(length(donor)),
    stats::rnorm(length(donor))
  )
  reference <- dplyr::bind_rows(lapply(seq_len(nrow(genes)), \(i)
    fit_peak_gene_hierarchical_association(genes[i, ], peak, depth)))
  compiled <- fit_peak_gene_hierarchical_batch(genes, peak, depth, kernel)
  testthat::expect_identical(compiled$hierarchical_status, reference$hierarchical_status)
  testthat::expect_equal(compiled$hierarchical_pvalue, reference$hierarchical_pvalue, tolerance = 1e-8)
  testthat::expect_equal(compiled$hierarchical_df, reference$hierarchical_df, tolerance = 1e-8)
  reversed <- fit_peak_gene_hierarchical_batch(genes[3:1, ], peak, depth, kernel)
  testthat::expect_equal(reversed$hierarchical_pvalue[3:1], compiled$hierarchical_pvalue, tolerance = 1e-10)
  uninformative <- fit_peak_gene_hierarchical_batch(genes, rep(0, length(peak)), depth, kernel)
  testthat::expect_true(all(uninformative$hierarchical_status == "insufficient_donor_information"))
  X <- cbind(1, peak)
  Z <- stats::model.matrix(~ 0 + factor(donor)) * peak
  covariance <- solve(crossprod(X))
  batch <- kernel(X, Z, list(covariance, covariance), c(0.4, 0.4), c(0, -1), c(1, 1))
  testthat::expect_true(is.finite(batch$p[[1]]))
  testthat::expect_true(is.na(batch$p[[2]]))
  testthat::expect_match(batch$diagnostic[[2]], "Invalid")
  testthat::expect_error(kernel(X, Z, list(covariance), numeric(), 0, 1), "dimensions")
})

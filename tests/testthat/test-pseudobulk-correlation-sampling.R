source_project_file("R/pseudobulk_helpers.R")

make_correlation_sampling_fixture <- function() {
  withr::with_seed(19L, {
    block <- rep(seq_len(8L), each = 3L)
    exposure <- rep(c(-1, 0, 1), 8L)
    design <- stats::model.matrix(~ factor(exposure))
    family_effects <- matrix(stats::rnorm(200L * 8L, sd = 0.5), 200L, 8L)
    means <- 40 * exp(family_effects[, block])
    counts <- matrix(stats::rpois(length(means), means), nrow(means)) + 1L
    dimnames(counts) <- list(paste0("feature_", seq_len(200L)), paste0("sample_", seq_len(24L)))
    rownames(design) <- colnames(counts)
    list(counts = counts, design = design, block = block)
  })
}

testthat::test_that("unlimited correlation fitting preserves edgeR including sparse-count DF", {
  input <- make_correlation_sampling_fixture()
  input$counts[1:5, seq.int(1L, 24L, by = 3L)] <- 0L
  reference <- do.call(edgeR::voomLmFit, input)
  unlimited <- do.call(fit_psbulk_voom, input)
  all_features <- do.call(fit_psbulk_voom, c(input, list(correlation_max_features = 200L)))
  sampled <- do.call(fit_psbulk_voom, c(input, list(correlation_max_features = 60L)))

  testthat::expect_equal(unlimited, reference, tolerance = 1e-12)
  testthat::expect_equal(all_features, reference, tolerance = 1e-12)
  testthat::expect_identical(sampled$df.residual, reference$df.residual)
  testthat::expect_gt(length(unique(reference$df.residual)), 1L)
})

testthat::test_that("sampling changes only correlation estimation and restores RNG and namespaces", {
  input <- make_correlation_sampling_fixture()
  original_function <- edgeR::voomLmFit
  original_estimator <- get("duplicateCorrelation", envir = environment(original_function))
  withr::local_seed(53L)
  seed_before <- .Random.seed
  fit <- do.call(fit_psbulk_voom, c(input, list(correlation_max_features = 60L)))
  testthat::expect_identical(.Random.seed, seed_before)
  rows <- withr::with_seed(732L, sort(sample.int(nrow(input$counts), 60L)))
  correlation <- limma::duplicateCorrelation(
    fit$EList$E[rows, , drop = FALSE], design = input$design, block = input$block,
    weights = fit$EList$weights[rows, , drop = FALSE]
  )$consensus.correlation

  testthat::expect_equal(fit$correlation, correlation, tolerance = 1e-12)
  testthat::expect_identical(rownames(fit$coefficients), rownames(input$counts))
  testthat::expect_identical(dim(fit$EList$weights), dim(input$counts))
  testthat::expect_identical(edgeR::voomLmFit, original_function)
  testthat::expect_identical(
    get("duplicateCorrelation", envir = environment(edgeR::voomLmFit)), original_estimator
  )
})

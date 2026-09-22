testthat::test_that("profiled donor-slope REML agrees with independently optimized fits", {
  native <- new.env(parent = baseenv())
  Rcpp::sourceCpp(file.path(multiomeR_project_root, "src/peak_gene_REML.cpp"), env = native)
  set.seed(1702)
  donor <- factor(rep(seq_len(6), each = 8))
  x <- stats::rnorm(length(donor))
  x <- x - ave(x, donor)
  design <- stats::model.matrix(~ 0 + donor)
  X <- cbind(design, x)
  Z <- design * x
  Y <- cbind(stats::rnorm(length(x)) + 0.5 * x,
    stats::rnorm(length(x)) + x * rep(seq(-2, 2, length.out = 6), each = 8))
  result <- native$peak_gene_REML_batch_cpp(X, Z, Y)
  for (i in seq_len(ncol(Y))) {
    fits <- lapply(c(0, 1, 5), function(start) suppressMessages(lme4::lmer(
      y ~ 0 + design + x + (0 + x | donor), data = data.frame(y = Y[, i], x, donor),
      REML = TRUE, start = list(theta = start),
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-10)))))
    reference <- fits[[which.max(vapply(fits, function(fit) as.numeric(stats::logLik(fit)), numeric(1)))]]
    testthat::expect_equal(result$REML_deviance[[i]], -2 * as.numeric(stats::logLik(reference)), tolerance = 1e-8)
    testthat::expect_equal(result$beta[, i], unname(lme4::fixef(reference)), tolerance = 1e-6)
    testthat::expect_equal(result$residual_variance[[i]], stats::sigma(reference)^2, tolerance = 1e-6)
    testthat::expect_equal(result$slope_variance[[i]], as.numeric(lme4::VarCorr(reference)[[1]]), tolerance = 1e-6)
    testthat::expect_equal(result$covariance[[i]], unname(as.matrix(stats::vcov(reference))), tolerance = 1e-6)
  }
  reversed <- native$peak_gene_REML_batch_cpp(X, Z, Y[, 2:1])
  testthat::expect_equal(reversed$beta[, 2:1], result$beta, tolerance = 1e-12)
  invalid <- native$peak_gene_REML_batch_cpp(X, Z, cbind(Y, 0))
  testthat::expect_equal(invalid$beta[, 1:2], result$beta, tolerance = 1e-12)
  testthat::expect_match(invalid$diagnostic[[3]], "No residual response variation")
  testthat::expect_error(native$peak_gene_REML_batch_cpp(cbind(X, x), Z, Y), "rank deficient")
})

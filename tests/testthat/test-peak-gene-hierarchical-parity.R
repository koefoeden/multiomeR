source_project_file("R/peak_gene_correlation_helpers.R")
source_project_file("R/peak_gene_hierarchical_helpers.R")

# Reference donor-slope fit: lme4 REML with pbkrtest Kenward-Roger inference.
#
# Donor fixed intercepts and log-depth covariates match the existing adjustment
# space. A random slope models donor heterogeneity; Kenward-Roger inference
# tests the average slope. Diagnostic p-values are retained when inference is
# unreliable. A zero variance estimate alone does not invalidate a fit.
fit_peak_gene_hierarchical_association <- function(gene, peak, aggregate_depth_tibble) {
  result <- tibble::tibble(
    hierarchical_coefficient = NA_real_,
    hierarchical_slope_SD = NA_real_,
    hierarchical_df = NA_real_,
    hierarchical_raw_pvalue = NA_real_,
    hierarchical_pvalue = NA_real_,
    hierarchical_singular = NA,
    hierarchical_status = "insufficient_donor_information",
    hierarchical_diagnostic = ""
  )
  donor <- factor(aggregate_depth_tibble$donor_id)
  x <- peak - ave(peak, donor)
  donor_ss <- rowsum(matrix(x^2, ncol = 1), donor)[, 1]
  # This checks identifiability of slope heterogeneity, not donor significance.
  if (sum(donor_ss > 1e-12 * max(sum(x^2), 1)) < 2L) {
    result$hierarchical_diagnostic <- "Within-donor peak variation in fewer than two donors"
    return(result)
  }
  design <- make_peak_gene_correlation_design_matrix(aggregate_depth_tibble)
  data <- data.frame(y = gene, x = x, donor = donor)
  diagnostics <- character()
  tryCatch(withCallingHandlers({
    fit <- suppressMessages(lme4::lmer(
      y ~ 0 + design + x + (0 + x | donor), data = data, REML = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa")
    ))
    result$hierarchical_coefficient <- unname(lme4::fixef(fit)["x"])
    result$hierarchical_slope_SD <- unname(attr(lme4::VarCorr(fit)$donor, "stddev"))
    result$hierarchical_singular <- lme4::isSingular(fit)
    restriction <- matrix(as.numeric(names(lme4::fixef(fit)) == "x"), nrow = 1)
    test <- pbkrtest::KRmodcomp(fit, restriction)
    result$hierarchical_df <- test$stats$ddf
    result$hierarchical_raw_pvalue <- test$stats$p.value
    # df < 1 is a conservative prototype reliability flag, not a validated
    # calibration boundary. Never interpret those extreme approximations.
    result$hierarchical_status <- if (
      !is.finite(test$stats$ddf) || test$stats$ddf < 1
    ) "insufficient_donor_information" else if (
      length(diagnostics) || !is.finite(test$stats$p.value) ||
      test$stats$p.value < 0 || test$stats$p.value > 1
    ) "numerical_failure" else "estimable"
    if (result$hierarchical_status == "estimable") {
      result$hierarchical_pvalue <- test$stats$p.value
    }
  }, warning = function(w) {
    diagnostics <<- c(diagnostics, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    result$hierarchical_status <<- "numerical_failure"
    diagnostics <<- c(diagnostics, conditionMessage(e))
  })
  result$hierarchical_diagnostic <- paste(unique(diagnostics), collapse = "; ")
  result
}

make_hierarchical_scan_case <- function() {
  set.seed(1703)
  n_donors <- 8L
  n_states <- 6L
  donor <- rep(sprintf("donor%d", seq_len(n_donors)), each = n_states)
  n <- length(donor)
  depth <- tibble::tibble(donor_id = donor, state_bin = as.character(rep(seq_len(n_states), n_donors)),
    GEX_depth = exp(stats::rnorm(n, 8, 0.3)), ATAC_depth = exp(stats::rnorm(n, 9, 0.3)))
  peaks <- rbind(common = stats::rnorm(n), second = stats::rnorm(n),
    one_donor = ifelse(donor == "donor1", stats::rnorm(n), 0))
  donor_effect <- stats::rnorm(n_donors)[factor(donor)]
  slopes <- list(null = rep(0, n_donors), shared = rep(0.5, n_donors),
    heterogeneous = seq(-0.5, 1.5, length.out = n_donors), strong = seq(0.8, 1.2, length.out = n_donors))
  genes <- t(vapply(slopes, function(slope) {
    donor_effect + slope[factor(donor)] * peaks["common", ] + 0.3 * peaks["second", ] +
      0.2 * log(depth$GEX_depth) + stats::rnorm(n, sd = 0.5)
  }, numeric(n)))
  list(
    matrices = list(cell_group = "fixture", chr = "chr1", GEX_norm = genes, ATAC_norm = peaks,
      aggregate_depth_tibble = depth,
      gene_detected_frac = stats::setNames(rep(1, nrow(genes)), rownames(genes)),
      peak_accessible_frac = stats::setNames(rep(1, nrow(peaks)), rownames(peaks)),
      mean_gene_expression = rowMeans(genes), mean_peak_accessibility = rowMeans(peaks)),
    pairs = tidyr::expand_grid(gene_matrix_feature = rownames(genes), peak = rownames(peaks)) |>
      dplyr::mutate(chr = "chr1", TargetGeneID = .data$gene_matrix_feature)
  )
}

testthat::test_that("the compiled REML kernel matches multi-start lme4 fits", {
  require_reference_version("lme4", "2.0.1")
  kernel <- load_peak_gene_kernel(file.path(multiomeR_project_root, "src/peak_gene_REML.cpp"), "peak_gene_REML_batch_cpp")
  set.seed(1702)
  donor <- factor(rep(seq_len(6), each = 8))
  x <- stats::rnorm(length(donor))
  x <- x - ave(x, donor)
  design <- stats::model.matrix(~ 0 + donor)
  X <- cbind(design, x)
  Z <- design * x
  Y <- cbind(stats::rnorm(length(x)) + 0.5 * x,
    stats::rnorm(length(x)) + x * rep(seq(-2, 2, length.out = 6), each = 8))
  result <- kernel(X, Z, Y)
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
})

testthat::test_that("the hierarchical scan matches lme4 with pbkrtest Kenward-Roger inference", {
  require_reference_version("lme4", "2.0.1")
  require_reference_version("pbkrtest", "0.5.5")
  case <- make_hierarchical_scan_case()
  observed <- score_peak_gene_hierarchical_associations(case$matrices, case$pairs,
    REML_source_file = file.path(multiomeR_project_root, "src/peak_gene_REML.cpp"),
    KR_source_file = file.path(multiomeR_project_root, "src/peak_gene_KR.cpp"))
  reference <- dplyr::bind_rows(lapply(seq_len(nrow(observed)), function(i) fit_peak_gene_hierarchical_association(
    case$matrices$GEX_norm[observed$TargetGeneID[[i]], ], case$matrices$ATAC_norm[observed$peak[[i]], ],
    case$matrices$aggregate_depth_tibble)))

  testthat::expect_identical(observed$hierarchical_status, reference$hierarchical_status)
  testthat::expect_setequal(observed$hierarchical_status, c("estimable", "insufficient_donor_information"))
  testthat::expect_true(any(observed$hierarchical_singular, na.rm = TRUE))
  testthat::expect_true(any(!observed$hierarchical_singular, na.rm = TRUE))
  for (column in c("hierarchical_coefficient", "hierarchical_slope_SD", "hierarchical_df",
    "hierarchical_raw_pvalue", "hierarchical_pvalue")) {
    testthat::expect_equal(observed[[column]], reference[[column]], tolerance = 1e-6, info = column)
  }
})

source_project_file("R/peak_gene_correlation_helpers.R")
source_project_file("R/peak_gene_KR_helpers.R")

# Fit an exploratory donor-varying peak-gene slope
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

# Apply the compiled correction to fresh models sharing one peak and design
#
# The final fixed-effect coefficient must be the peak slope named `x`.
# Variance estimates are extracted separately for each fitted response.
compute_peak_gene_KR_batch <- function(fits, kernel) {
  X <- lme4::getME(fits[[1]], "X")
  Z <- as.matrix(t(lme4::getME(fits[[1]], "Zt")))
  stopifnot(tail(colnames(X), 1) == "x")
  for (fit in fits) {
    stopifnot(identical(lme4::getME(fit, "X"), X),
      identical(as.matrix(t(lme4::getME(fit, "Zt"))), Z),
      length(lme4::VarCorr(fit)) == 1L,
      length(lme4::VarCorr(fit)[[1]]) == 1L,
      lme4::isREML(fit))
  }
  kernel(
    X, Z,
    lapply(fits, \(fit) as.matrix(stats::vcov(fit))),
    vapply(fits, \(fit) unname(lme4::fixef(fit)["x"]), numeric(1)),
    vapply(fits, \(fit) as.numeric(lme4::VarCorr(fit)[[1]]), numeric(1)),
    vapply(fits, \(fit) stats::sigma(fit)^2, numeric(1))
  )
}

# Fit fresh donor-slope models and correct all genes for one peak in a batch
fit_peak_gene_hierarchical_batch <- function(gene_matrix, peak, aggregate_depth_tibble, kernel) {
  donor <- factor(aggregate_depth_tibble$donor_id)
  x <- peak - ave(peak, donor)
  donor_ss <- rowsum(matrix(x^2, ncol = 1), donor)[, 1]
  result <- fit_peak_gene_hierarchical_association(0, 0, data.frame(donor_id = "empty"))[
    rep(1L, nrow(gene_matrix)), ]
  if (sum(donor_ss > 1e-12 * max(sum(x^2), 1)) < 2L) return(result)
  design <- make_peak_gene_correlation_design_matrix(aggregate_depth_tibble)
  fits <- lapply(seq_len(nrow(gene_matrix)), \(i) {
    diagnostics <- character()
    fit <- tryCatch(withCallingHandlers(
      suppressMessages(lme4::lmer(
        y ~ 0 + design + x + (0 + x | donor),
        data = data.frame(y = as.numeric(gene_matrix[i, ]), x = x, donor = donor),
        REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa")
      )), warning = function(w) {
        diagnostics <<- c(diagnostics, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        diagnostics <<- c(diagnostics, conditionMessage(e))
        NULL
      })
    result$hierarchical_diagnostic[[i]] <<- paste(unique(diagnostics), collapse = "; ")
    result$hierarchical_status[[i]] <<- "numerical_failure"
    fit
  })
  indices <- which(!vapply(fits, is.null, logical(1)))
  if (!length(indices)) return(result)
  correction <- compute_peak_gene_KR_batch(fits[indices], kernel)
  for (j in seq_along(indices)) {
    i <- indices[[j]]
    fit <- fits[[i]]
    result$hierarchical_coefficient[[i]] <- unname(lme4::fixef(fit)["x"])
    result$hierarchical_slope_SD[[i]] <- unname(attr(lme4::VarCorr(fit)$donor, "stddev"))
    result$hierarchical_singular[[i]] <- lme4::isSingular(fit)
    result$hierarchical_df[[i]] <- correction$df[[j]]
    result$hierarchical_raw_pvalue[[i]] <- correction$p[[j]]
    diagnostics <- c(result$hierarchical_diagnostic[[i]], correction$diagnostic[[j]])
    result$hierarchical_diagnostic[[i]] <- paste(diagnostics[nzchar(diagnostics)], collapse = "; ")
    status <- if (nzchar(correction$diagnostic[[j]])) {
      "numerical_failure"
    } else if (!is.finite(correction$df[[j]]) || correction$df[[j]] < 1) {
      "insufficient_donor_information"
    } else if (nzchar(result$hierarchical_diagnostic[[i]]) || !is.finite(correction$p[[j]])) {
      "numerical_failure"
    } else "estimable"
    result$hierarchical_status[[i]] <- status
    if (status == "estimable") result$hierarchical_pvalue[[i]] <- correction$p[[j]]
  }
  result
}

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

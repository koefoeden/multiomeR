#' Fit an exploratory donor-varying peak-gene slope
#'
#' Donor fixed intercepts and log-depth covariates match the existing adjustment
#' space. A random slope models donor heterogeneity; Kenward-Roger inference
#' tests the average slope. Diagnostic p-values are retained when inference is
#' unreliable. A zero variance estimate alone does not invalidate a fit.
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

#' Load the profiled REML kernel once per source revision in each worker
load_peak_gene_REML_kernel <- function(native_source_file) {
  environment <- new.env(parent = baseenv())
  Rcpp::sourceCpp(native_source_file, env = environment,
    cacheDir = file.path(tempdir(), "peak_gene_REML"), showOutput = FALSE)
  environment$peak_gene_REML_batch_cpp
}

#' Scan every eligible pair in one cell-type/chromosome branch
#'
#' Measurement-support candidates, detection and aggregate eligibility match the HC3 scan. No
#' correlation, HC3 p-value, promoter or top-N screen is applied before fitting.
#' Dense matrices and nuisance design are reused across peak batches.
score_peak_gene_hierarchical_associations <- function(
  normalized_aggregate_matrices, candidate_pairs_tibble,
  REML_source_file, KR_source_file
) {
  matrices <- normalized_aggregate_matrices
  branch <- prepare_peak_gene_correlation_branch(matrices, candidate_pairs_tibble)
  pairs <- branch$candidate_pairs
  if (!is.na(branch$skipped_reason)) pairs <- pairs[0, ]
  result <- data.frame(
    cell_group = rep(branch$cell_group, nrow(pairs)),
    peak = pairs$peak, TargetGeneID = pairs$TargetGeneID,
    hierarchical_coefficient = rep(NA_real_, nrow(pairs)), hierarchical_slope_SD = rep(NA_real_, nrow(pairs)),
    hierarchical_df = rep(NA_real_, nrow(pairs)), hierarchical_raw_pvalue = rep(NA_real_, nrow(pairs)),
    hierarchical_pvalue = rep(NA_real_, nrow(pairs)), hierarchical_singular = rep(NA, nrow(pairs)),
    hierarchical_status = rep("insufficient_donor_information", nrow(pairs)),
    hierarchical_diagnostic = rep("Within-donor peak variation in fewer than two donors", nrow(pairs))
  )
  if (!nrow(pairs)) return(tibble::as_tibble(result))
  donor <- factor(matrices$aggregate_depth_tibble$donor_id)
  if (nlevels(donor) < 2L) return(tibble::as_tibble(result))
  fit_kernel <- load_peak_gene_REML_kernel(REML_source_file)
  correction_kernel <- load_peak_gene_KR_kernel(KR_source_file)
  donor_indices <- split(seq_along(donor), donor)
  donor_matrix <- stats::model.matrix(~ 0 + donor)
  design <- branch$design
  design_qr <- qr(design)
  GEX <- as.matrix(matrices$GEX_norm)
  ATAC <- as.matrix(matrices$ATAC_norm)
  groups <- split(seq_len(nrow(pairs)), pairs$peak)
  for (indices in groups) {
    x <- as.numeric(ATAC[pairs$peak[[indices[[1]]]], ])
    x <- x - ave(x, donor)
    donor_ss <- vapply(donor_indices, \(i) sum(x[i]^2), numeric(1))
    if (sum(donor_ss > 1e-12 * max(sum(x^2), 1)) < 2L) next
    if (sum(qr.resid(design_qr, x)^2) <= 1e-12 * max(sum(x^2), 1)) {
      result$hierarchical_diagnostic[indices] <- "No peak variation after nuisance adjustment"
      next
    }
    X <- cbind(design, x)
    Z <- donor_matrix * x
    fit <- fit_kernel(X, Z, t(GEX[pairs$gene_matrix_feature[indices], , drop = FALSE]))
    correction <- correction_kernel(X, Z, fit$covariance, fit$beta[ncol(X), ],
      fit$slope_variance, fit$residual_variance)
    diagnostic <- ifelse(nzchar(fit$diagnostic), fit$diagnostic, correction$diagnostic)
    status <- ifelse(nzchar(diagnostic) | !is.finite(correction$p) |
      correction$p < 0 | correction$p > 1, "numerical_failure",
      ifelse(!is.finite(correction$df) | correction$df < 1,
        "insufficient_donor_information", "estimable"))
    result$hierarchical_coefficient[indices] <- fit$beta[ncol(X), ]
    result$hierarchical_slope_SD[indices] <- sqrt(fit$slope_variance)
    result$hierarchical_singular[indices] <- sqrt(fit$variance_ratio) < 1e-4
    result$hierarchical_df[indices] <- correction$df
    result$hierarchical_raw_pvalue[indices] <- correction$p
    result$hierarchical_pvalue[indices] <- ifelse(status == "estimable", correction$p, NA_real_)
    result$hierarchical_status[indices] <- status
    result$hierarchical_diagnostic[indices] <- diagnostic
  }
  tibble::as_tibble(result)
}

#' Adjust across the complete candidate family within each cell type
#'
#' Unreliable tests retain NA p/FDR but count in the BH family size, so numerical
#' exclusions do not make the correction less stringent. Direction and promoter
#' filtering happen only after this adjustment.
finalize_peak_gene_hierarchical_results <- function(results_tibbles) {
  result <- data.table::rbindlist(results_tibbles, use.names = TRUE)
  result[, hierarchical_FDR := stats::p.adjust(hierarchical_pvalue, method = "BH", n = .N), by = cell_group]
  data.table::setDF(result)
  tibble::as_tibble(result)
}

#' Select positive nonpromoter hierarchical links without an HC3 dependency
select_peak_gene_hierarchical_top_links <- function(
  hierarchical_results, candidate_pairs_tibble, n_per_cell_group
) {
  selected <- hierarchical_results |>
    dplyr::filter(.data$hierarchical_status == "estimable", .data$hierarchical_coefficient > 0) |>
    dplyr::semi_join(dplyr::filter(candidate_pairs_tibble, !.data$isSelfPromoter),
      by = c("peak", "TargetGeneID")) |>
    dplyr::arrange(.data$cell_group, .data$hierarchical_pvalue, .data$TargetGeneID, .data$peak) |>
    dplyr::slice_head(n = n_per_cell_group, by = "cell_group") |>
    dplyr::left_join(candidate_pairs_tibble, by = c("peak", "TargetGeneID"),
      relationship = "many-to-one")
  make_peak_gene_correlation_top_links(selected, candidate_pairs_tibble, n_per_cell_group)
}

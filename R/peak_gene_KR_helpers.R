#' Load the peak-gene-specific batched Kenward-Roger kernel
#'
#' Each worker compiles at most once per source revision into its own temporary
#' cache, avoiding concurrent writes to shared build artifacts.
load_peak_gene_KR_kernel <- function(native_source_file) {
  environment <- new.env(parent = baseenv())
  Rcpp::sourceCpp(
    native_source_file, env = environment,
    cacheDir = file.path(tempdir(), "peak_gene_KR"), showOutput = FALSE
  )
  environment$peak_gene_KR_batch_cpp
}

#' Apply the compiled correction to fresh models sharing one peak and design
#'
#' The final fixed-effect coefficient must be the peak slope named `x`.
#' Variance estimates are extracted separately for each fitted response.
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

#' Fit fresh donor-slope models and correct all genes for one peak in a batch
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

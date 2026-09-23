#' Load the profiled REML kernel once per source revision in each worker
load_peak_gene_REML_kernel <- function(native_source_file) {
  environment <- new.env(parent = baseenv())
  Rcpp::sourceCpp(native_source_file, env = environment,
    cacheDir = file.path(tempdir(), "peak_gene_REML"), showOutput = FALSE)
  environment$peak_gene_REML_batch_cpp
}

#' Scan every eligible pair in one cell-type/chromosome branch
#'
#' Measurement-support candidates, detection and aggregate eligibility come from
#' the branch preparation. No correlation, promoter or top-N screen is applied
#' before fitting. Dense matrices and nuisance design are reused across peak
#' batches. The descriptive `correlation` is the Pearson correlation of the
#' donor- and depth-adjusted peak and gene values.
score_peak_gene_hierarchical_associations <- function(
  normalized_aggregate_matrices, candidate_pairs_tibble,
  REML_source_file, KR_source_file
) {
  matrices <- normalized_aggregate_matrices
  branch <- prepare_peak_gene_correlation_branch(matrices, candidate_pairs_tibble)
  pairs <- branch$candidate_pairs
  if (!is.na(branch$skipped_reason)) pairs <- pairs[0, ]
  n_pairs <- nrow(pairs)
  result <- data.frame(
    cell_group = rep(branch$cell_group, n_pairs), pairs,
    n_aggregates = rep(branch$n_aggregates, n_pairs), n_donors = rep(branch$n_donors, n_pairs),
    residual_df = rep(branch$residual_df, n_pairs),
    gene_detected_frac = unname(as.numeric(matrices$gene_detected_frac[pairs$gene_matrix_feature])),
    peak_accessible_frac = unname(as.numeric(matrices$peak_accessible_frac[pairs$peak])),
    correlation = rep(NA_real_, n_pairs),
    hierarchical_coefficient = rep(NA_real_, n_pairs), hierarchical_slope_SD = rep(NA_real_, n_pairs),
    hierarchical_df = rep(NA_real_, n_pairs), hierarchical_raw_pvalue = rep(NA_real_, n_pairs),
    hierarchical_pvalue = rep(NA_real_, n_pairs), hierarchical_singular = rep(NA, n_pairs),
    hierarchical_status = rep("insufficient_donor_information", n_pairs),
    hierarchical_diagnostic = rep("Within-donor peak variation in fewer than two donors", n_pairs),
    check.names = FALSE
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
  GEX_residual <- qr.resid(design_qr, t(GEX))
  GEX_residual_ss <- colSums(GEX_residual^2)
  groups <- split(seq_len(nrow(pairs)), pairs$peak)
  for (indices in groups) {
    x <- as.numeric(ATAC[pairs$peak[[indices[[1]]]], ])
    x <- x - ave(x, donor)
    x_residual <- qr.resid(design_qr, x)
    x_residual_ss <- sum(x_residual^2)
    genes <- pairs$gene_matrix_feature[indices]
    correlation <- crossprod(GEX_residual[, genes, drop = FALSE], x_residual)[, 1] /
      sqrt(x_residual_ss * GEX_residual_ss[genes])
    correlation[!is.finite(correlation) | GEX_residual_ss[genes] <= 1e-12 | x_residual_ss <= 1e-12] <- NA_real_
    result$correlation[indices] <- pmax(pmin(correlation, 1), -1)
    donor_ss <- vapply(donor_indices, \(i) sum(x[i]^2), numeric(1))
    if (sum(donor_ss > 1e-12 * max(sum(x^2), 1)) < 2L) next
    if (x_residual_ss <= 1e-12 * max(sum(x^2), 1)) {
      result$hierarchical_diagnostic[indices] <- "No peak variation after nuisance adjustment"
      next
    }
    X <- cbind(design, x)
    Z <- donor_matrix * x
    fit <- fit_kernel(X, Z, t(GEX[genes, , drop = FALSE]))
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

#' Candidate enhancer links: positive, reliable, FDR-significant nonpromoter pairs.
is_peak_gene_link <- function(results_tibble) {
  results_tibble$hierarchical_status == "estimable" & results_tibble$hierarchical_coefficient > 0 &
    results_tibble$hierarchical_FDR < 0.05 & !results_tibble$isSelfPromoter
}

make_peak_gene_correlation_links <- function(results_tibble) {
  results_tibble[which(is_peak_gene_link(results_tibble)), ] |>
    dplyr::arrange(.data$cell_group, .data$hierarchical_pvalue, .data$TargetGeneID, .data$peak) |>
    dplyr::mutate(rank_in_cell_group = dplyr::row_number(), .by = "cell_group")
}

#' Select the lowest reliable p-values among positive nonpromoter pairs
#'
#' Returns one placeholder row when no link is analyzable, so downstream plot
#' targets can still branch over the result.
select_peak_gene_hierarchical_top_links <- function(hierarchical_results, n_per_cell_group) {
  top_links <- hierarchical_results |>
    dplyr::filter(.data$hierarchical_status == "estimable", .data$hierarchical_coefficient > 0,
      !.data$isSelfPromoter) |>
    dplyr::arrange(.data$cell_group, .data$hierarchical_pvalue, .data$TargetGeneID, .data$peak) |>
    dplyr::slice_head(n = n_per_cell_group, by = "cell_group") |>
    dplyr::mutate(rank_in_cell_group = dplyr::row_number(), .by = "cell_group") |>
    dplyr::mutate(
      scatter_plot_name = paste(
        make.names(.data$cell_group),
        sprintf("rank%03d", .data$rank_in_cell_group),
        make.names(.data$TargetGene),
        .data$chr,
        sep = "_"
      ),
      is_analyzable_link = TRUE
    )
  if (nrow(top_links) > 0L) {
    return(top_links)
  }
  top_links[NA_integer_, , drop = FALSE] |>
    dplyr::mutate(
      cell_group = "__no_analyzable_link__",
      scatter_plot_name = "no_analyzable_peak_gene_link",
      is_analyzable_link = FALSE
    )
}

#' Load a peak-gene Rcpp kernel once per source revision in each worker
#'
#' Each worker compiles into its own temporary cache, avoiding concurrent writes
#' to shared build artifacts.
load_peak_gene_kernel <- function(native_source_file, kernel_name) {
  environment <- new.env(parent = baseenv())
  Rcpp::sourceCpp(native_source_file, env = environment,
    cacheDir = file.path(tempdir(), kernel_name), showOutput = FALSE)
  environment[[kernel_name]]
}

#' Sample each gene's background peaks from its branch
#'
#' A gene that varies along a cell type's chromatin states correlates with many
#' peaks regardless of distance. Up to `n_per_gene` of the branch's candidate
#' peaks more than `min_distance` from the gene's TSS share this non-specific
#' covariation but no cis-regulatory relationship, so they calibrate the gene's
#' candidate pairs.
#' Like the background peaks of Signac's LinkPeaks() (Stuart et al. 2021, doi:10.1038/s41592-021-01282-5).
add_peak_gene_background_pairs <- function(pairs, n_per_gene = 200L, min_distance = 1e7) {
  peaks <- dplyr::distinct(pairs, .data$peak, .data$chr, .data$start, .data$end, .data$peak_center)
  genes <- dplyr::distinct(pairs, .data$TargetGeneID, .data$TargetGene, .data$gene_matrix_feature,
    .data$TargetGeneTSS, .data$TargetGeneStart, .data$TargetGeneEnd)
  draws <- lapply(genes$TargetGeneTSS, \(TSS) {
    eligible <- which(abs(peaks$peak_center - TSS) > min_distance)
    eligible[sample.int(length(eligible), min(n_per_gene, length(eligible)))]
  })
  background <- dplyr::bind_cols(peaks[unlist(draws), ], genes[rep(seq_len(nrow(genes)), lengths(draws)), ]) |>
    dplyr::mutate(distance = as.integer(.data$peak_center - .data$TargetGeneTSS), isSelfPromoter = FALSE,
      isTargetGeneBody = FALSE, link_class = "background")
  dplyr::bind_rows(pairs, background[names(pairs)])
}

#' Scan every eligible pair in one cell-type/chromosome branch
#'
#' Measurement-support candidates, detection and aggregate eligibility come from
#' the branch preparation. No correlation, promoter or top-N screen is applied
#' before fitting. Dense matrices and nuisance design are reused across peak
#' batches. The descriptive `correlation` is the Pearson correlation of the
#' nuisance-adjusted peak and gene values. Each gene's background pairs
#' (`link_class == "background"`) are scored with the same model, and
#' `hierarchical_gene_z` standardizes every pair's signed statistic by the mean
#' and SD of its gene's background statistics (leave-one-out for background
#' rows; missing with fewer than ten). Background rows keep only the columns the
#' calibration needs. Random slopes need `min_donors` donors with within-donor
#' peak variation.
score_peak_gene_hierarchical_associations <- function(
  normalized_aggregate_matrices, candidate_pairs_tibble,
  REML_source_file, KR_source_file, min_donors = 8L
) {
  matrices <- normalized_aggregate_matrices
  branch <- prepare_peak_gene_correlation_branch(matrices, candidate_pairs_tibble)
  pairs <- branch$candidate_pairs
  if (!is.na(branch$skipped_reason)) pairs <- pairs[0, ]
  pairs <- add_peak_gene_background_pairs(pairs)
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
    hierarchical_diagnostic = rep(paste("Within-donor peak variation in fewer than", min_donors, "donors"), n_pairs),
    hierarchical_gene_z = rep(NA_real_, n_pairs), check.names = FALSE
  )
  if (!nrow(pairs)) return(tibble::as_tibble(result))
  donor <- factor(matrices$aggregate_depth_tibble$donor_id)
  if (nlevels(donor) < min_donors) return(tibble::as_tibble(result))
  fit_kernel <- load_peak_gene_kernel(REML_source_file, "peak_gene_REML_batch_cpp")
  correction_kernel <- load_peak_gene_kernel(KR_source_file, "peak_gene_KR_batch_cpp")
  donor_index <- as.integer(donor)
  donor_counts <- tabulate(donor_index, nlevels(donor))
  donor_matrix <- stats::model.matrix(~ 0 + donor)
  design <- branch$design
  design_qr <- qr(design)
  GEX <- as.matrix(matrices$GEX_norm)
  ATAC <- as.matrix(matrices$ATAC_norm)
  GEX_residual <- qr.resid(design_qr, t(GEX))
  GEX_residual_ss <- colSums(GEX_residual^2)
  # Fill plain vectors per peak; assigning into data-frame columns copies them each time.
  estimates <- as.list(result[c("correlation", "hierarchical_coefficient", "hierarchical_slope_SD",
    "hierarchical_df", "hierarchical_raw_pvalue", "hierarchical_pvalue", "hierarchical_singular",
    "hierarchical_status", "hierarchical_diagnostic")])
  groups <- split(seq_len(nrow(pairs)), pairs$peak)
  for (indices in groups) {
    x <- as.numeric(ATAC[pairs$peak[[indices[[1]]]], ])
    x <- x - (rowsum(x, donor_index)[, 1] / donor_counts)[donor_index]
    x_residual <- qr.resid(design_qr, x)
    x_residual_ss <- sum(x_residual^2)
    genes <- pairs$gene_matrix_feature[indices]
    correlation <- crossprod(GEX_residual[, genes, drop = FALSE], x_residual)[, 1] /
      sqrt(x_residual_ss * GEX_residual_ss[genes])
    correlation[!is.finite(correlation) | GEX_residual_ss[genes] <= 1e-12 | x_residual_ss <= 1e-12] <- NA_real_
    estimates$correlation[indices] <- pmax(pmin(correlation, 1), -1)
    donor_ss <- rowsum(x^2, donor_index)[, 1]
    if (sum(donor_ss > 1e-12 * max(sum(x^2), 1)) < min_donors) next
    if (x_residual_ss <= 1e-12 * max(sum(x^2), 1)) {
      estimates$hierarchical_diagnostic[indices] <- "No peak variation after nuisance adjustment"
      next
    }
    # The donor-varying slope makes donors the units of inference; Kenward-Roger
    # (Kenward & Roger 1997) adjusts its SE and df for the few donors.
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
    estimates$hierarchical_coefficient[indices] <- fit$beta[ncol(X), ]
    estimates$hierarchical_slope_SD[indices] <- sqrt(fit$slope_variance)
    estimates$hierarchical_singular[indices] <- sqrt(fit$variance_ratio) < 1e-4
    estimates$hierarchical_df[indices] <- correction$df
    estimates$hierarchical_raw_pvalue[indices] <- correction$p
    estimates$hierarchical_pvalue[indices] <- ifelse(status == "estimable", correction$p, NA_real_)
    estimates$hierarchical_status[indices] <- status
    estimates$hierarchical_diagnostic[indices] <- diagnostic
  }
  result[names(estimates)] <- estimates
  statistic <- sign(result$hierarchical_coefficient) *
    stats::qnorm(pmax(result$hierarchical_pvalue, .Machine$double.xmin) / 2, lower.tail = FALSE)
  is_background <- pairs$link_class == "background"
  in_null <- is_background & is.finite(statistic)
  null_statistic <- ifelse(in_null, statistic, 0)
  gene_sum <- function(x) stats::ave(x, pairs$gene_matrix_feature, FUN = sum)
  n <- gene_sum(as.numeric(in_null)) - in_null
  gene_mean <- (gene_sum(null_statistic) - null_statistic) / n
  gene_SD <- sqrt(pmax(gene_sum(null_statistic^2) - null_statistic^2 - n * gene_mean^2, 0) / (n - 1))
  result$hierarchical_gene_z <- ifelse(n >= 10, (statistic - gene_mean) / gene_SD, NA_real_)
  dplyr::bind_rows(tibble::as_tibble(result[!is_background, ]), tibble::as_tibble(result[is_background, c("cell_group", "chr",
    "gene_matrix_feature", "link_class", "peak_accessible_frac", "hierarchical_coefficient", "hierarchical_pvalue",
    "hierarchical_status", "hierarchical_gene_z")]))
}

#' Calibrate against each gene's background pairs and adjust within each cell type
#'
#' The gene-standardized statistics are rescaled by the cell type's background SD
#' within quintiles of peak detection, because more accessible peaks have wider
#' null distributions. A pair's empirical P-value is the share of the cell type's
#' rescaled background statistics at least as extreme. Unreliable tests retain
#' NA p/FDR but count in the BH family size, so numerical exclusions do not make
#' the correction less stringent. Direction and promoter filtering happen only
#' after this adjustment; background rows are dropped.
finalize_peak_gene_hierarchical_results <- function(results_tibbles) {
  result <- data.table::rbindlist(results_tibbles, use.names = TRUE, fill = TRUE)
  result[, is_background := link_class == "background" & !is.na(hierarchical_gene_z)]
  result[, hierarchical_background_pvalue := if (any(is_background)) {
    quintile <- findInterval(peak_accessible_frac, stats::quantile(peak_accessible_frac[is_background], 1:4 / 5))
    scale <- tapply(hierarchical_gene_z[is_background], quintile[is_background], stats::sd)
    z <- abs(hierarchical_gene_z / scale[as.character(quintile)])
    background <- sort(z[is_background])
    (length(background) - findInterval(z, background, left.open = TRUE) + 1) / (length(background) + 1)
  } else NA_real_, by = cell_group]
  result[, n_background_pairs := sum(is_background), by = cell_group]
  result <- result[link_class != "background"][, is_background := NULL]
  result[, hierarchical_FDR := stats::p.adjust(hierarchical_background_pvalue, method = "BH", n = .N), by = cell_group]
  data.table::setDF(result)
  tibble::as_tibble(result)
}

#' Candidate enhancer links: positive, reliable, FDR-significant nonpromoter
#' pairs that exceed their gene's background.
is_peak_gene_link <- function(results_tibble) {
  results_tibble$hierarchical_status == "estimable" & results_tibble$hierarchical_coefficient > 0 &
    results_tibble$hierarchical_gene_z > 0 & results_tibble$hierarchical_FDR < 0.05 & !results_tibble$isSelfPromoter
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

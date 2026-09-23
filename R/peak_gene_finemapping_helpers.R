#' Build compact FDR breakpoints for chromosome-level peak gene fine-mapping
#'
#' BH-adjusted p-values are monotone in the p-value. Retaining the largest
#' p-value for each distinct FDR reconstructs the tested values without sending
#' the complete finalized result table to every worker.
#'
#' @param finalized_results_tibble Combined hierarchical association results.
#' @return Tibble of FDR breakpoints per cell group.
#' @keywords internal

make_peak_gene_finemapping_reference <- function(finalized_results_tibble) {
  finalized_results_tibble |>
    dplyr::select("cell_group", "hierarchical_pvalue", "hierarchical_FDR") |>
    dplyr::filter(is.finite(.data$hierarchical_pvalue), is.finite(.data$hierarchical_FDR)) |>
    dplyr::summarise(hierarchical_pvalue = max(.data$hierarchical_pvalue), .by = c("cell_group", "hierarchical_FDR")) |>
    dplyr::arrange(.data$cell_group, .data$hierarchical_pvalue)
}

#' Restore within-cell-group FDR values on one chromosome result
#'
#' @param results_tibble Hierarchical result from one cell-group/chromosome branch.
#' @param FDR_reference_tibble Breakpoints from `make_peak_gene_finemapping_reference()`.
#' @return Input result with its within-cell-group `hierarchical_FDR` restored.
#' @keywords internal

restore_peak_gene_correlation_FDR <- function(results_tibble, FDR_reference_tibble) {
  FDR <- rep(NA_real_, nrow(results_tibble))
  if (nrow(results_tibble) > 0L) {
    reference <- FDR_reference_tibble[FDR_reference_tibble$cell_group == results_tibble$cell_group[[1]], ]
    valid <- is.finite(results_tibble$hierarchical_pvalue)
    index <- findInterval(results_tibble$hierarchical_pvalue[valid], reference$hierarchical_pvalue, left.open = TRUE) + 1L
    FDR[valid] <- reference$hierarchical_FDR[index]
  }
  dplyr::mutate(results_tibble, hierarchical_FDR = FDR)
}

#' Empty peak gene fine-mapping result
#'
#' @return A typed empty tibble for compact conditional peak prioritization.
#' @keywords internal

empty_peak_gene_correlation_finemapping_tibble <- function() {
  tibble::tibble(
    cell_group = character(),
    chr = character(),
    TargetGeneID = character(),
    TargetGene = character(),
    peak = character(),
    correlation = numeric(),
    FDR = numeric(),
    susie_PIP = numeric(),
    credible_set = integer(),
    in_credible_set = logical(),
    n_candidate_peaks = integer(),
    n_aggregates = integer(),
    n_donors = integer(),
    converged = logical(),
    finemapping_method = character()
  )
}

#' Prioritize conditionally supported peak-gene links with SuSiE
#'
#' Fine-map detected cis peaks for genes that contain at least one candidate
#' enhancer link. This is multipeak prioritization, not causal proof: the SuSiE
#' likelihood does not itself model donor clustering.
#'
#' @param normalized_aggregate_matrices Donor-state pseudobulk matrices.
#' @param branch_results Hierarchical results of the same branch with
#'   `hierarchical_FDR` restored.
#' @param max_genes Maximum screened genes fine-mapped in one branch.
#' @param max_candidate_peaks Maximum peaks retained per gene, prioritized by
#'   absolute donor-adjusted correlation when necessary.
#' @param L Maximum SuSiE single-effect components.
#' @return Compact peak-level PIP and credible-set records.
#' @keywords internal

finemap_peak_gene_correlations_for_branch <- function(
  normalized_aggregate_matrices,
  branch_results,
  max_genes = 50L,
  max_candidate_peaks = 500L,
  L = 10L
) {
  screened_genes <- branch_results[which(is_peak_gene_link(branch_results)), ] |>
    dplyr::summarise(
      best_FDR = min(.data$hierarchical_FDR),
      .by = c("TargetGeneID", "TargetGene")
    ) |>
    dplyr::slice_min(.data$best_FDR, n = max_genes, with_ties = FALSE)
  if (nrow(screened_genes) == 0L) {
    return(empty_peak_gene_correlation_finemapping_tibble())
  }

  GEX_norm <- normalized_aggregate_matrices$GEX_norm
  ATAC_norm <- normalized_aggregate_matrices$ATAC_norm
  design <- make_peak_gene_correlation_design_matrix(
    normalized_aggregate_matrices$aggregate_depth_tibble
  )
  design_qr <- qr(design)
  n_aggregates <- ncol(GEX_norm)
  n_donors <- dplyr::n_distinct(
    normalized_aggregate_matrices$aggregate_depth_tibble$donor_id
  )

  purrr::map_dfr(seq_len(nrow(screened_genes)), \(gene_index) {
    gene_id <- screened_genes$TargetGeneID[[gene_index]]
    gene_results <- branch_results |>
      dplyr::filter(.data$TargetGeneID == !!gene_id) |>
      dplyr::arrange(dplyr::desc(abs(.data$correlation))) |>
      dplyr::slice_head(n = max_candidate_peaks)
    if (nrow(gene_results) < 2L) {
      return(empty_peak_gene_correlation_finemapping_tibble())
    }

    peak_residual <- residualize_peak_gene_correlation_matrix(
      ATAC_norm[gene_results$peak, , drop = FALSE],
      design
    )
    gene_residual <- as.numeric(qr.resid(
      design_qr,
      as.numeric(GEX_norm[gene_results$gene_matrix_feature[[1]], ])
    ))
    variable_peaks <- rowSums(peak_residual^2) > .Machine$double.eps
    peak_residual <- peak_residual[variable_peaks, , drop = FALSE]
    gene_results <- gene_results[variable_peaks, , drop = FALSE]
    if (nrow(gene_results) < 2L || sum(gene_residual^2) <= .Machine$double.eps) {
      return(empty_peak_gene_correlation_finemapping_tibble())
    }

    fit <- susieR::susie(
      X = t(peak_residual),
      y = gene_residual,
      L = min(as.integer(L), nrow(peak_residual)),
      intercept = FALSE,
      standardize = TRUE,
      estimate_residual_variance = TRUE,
      max_iter = 100L,
      verbose = FALSE
    )
    credible_sets <- susieR::susie_get_cs(fit, coverage = 0.95)$cs
    credible_set <- rep(NA_integer_, nrow(gene_results))
    if (length(credible_sets) > 0L) {
      for (set_index in seq_along(credible_sets)) {
        peak_indices <- credible_sets[[set_index]]
        unassigned_indices <- peak_indices[is.na(credible_set[peak_indices])]
        credible_set[unassigned_indices] <- set_index
      }
    }

    prioritized <- gene_results |>
      dplyr::transmute(
        cell_group = .data$cell_group,
        chr = .data$chr,
        TargetGeneID = .data$TargetGeneID,
        TargetGene = .data$TargetGene,
        peak = .data$peak,
        correlation = as.numeric(.data$correlation),
        FDR = .data$hierarchical_FDR,
        susie_PIP = as.numeric(fit$pip),
        credible_set = credible_set,
        in_credible_set = !is.na(.data$credible_set),
        n_candidate_peaks = nrow(gene_results),
        n_aggregates = n_aggregates,
        n_donors = n_donors,
        converged = isTRUE(fit$converged),
        finemapping_method = "SuSiE on donor/depth residuals"
      )
    retained <- prioritized |>
      dplyr::filter(.data$susie_PIP >= 0.01 | .data$in_credible_set)
    if (nrow(retained) == 0L) {
      retained <- prioritized |>
        dplyr::slice_max(.data$susie_PIP, n = 1L, with_ties = FALSE)
    }
    retained
  })
}

#' Filter credible set variants
#'
#' Filter and optionally reweight credible-set variants before peak overlap.
#'
#' @param credible_set_GRanges GRanges object containing credible set GRanges coordinates and metadata.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param ... Additional arguments forwarded to the variant weighting function.
#' @return A GRanges object containing retained variants, with any weighting
#'   function side effects applied to metadata columns.
#' @keywords internal

filter_credible_set_variants <- function(credible_set_GRanges, posterior_probability_cutoff = NULL, posterior_probability_weighting_function = NULL, ...) {
  if (!is.null(posterior_probability_cutoff)) {
    credible_set_GRanges <- S4Vectors::subset(credible_set_GRanges, posteriorProbability > posterior_probability_cutoff)
  }

  if (!is.null(posterior_probability_weighting_function)) {
    credible_set_GRanges <- posterior_probability_weighting_function(credible_set_GRanges, ...)
  }

  credible_set_GRanges
}

#' Sum variant weights in peaks
#'
#' Sum overlapping variant weights for each peak range.
#'
#' @param variant_GRanges GRanges object containing variant GRanges coordinates and metadata.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param weight_col Column in `variant_tibble` to draw as vertical PIP/weight values.
#' @return Named numeric vector with one value per peak; peaks without variants
#'   receive zero.
#' @keywords internal

sum_variant_weights_in_peaks <- function(variant_GRanges, peak_ranges, weight_col = "posteriorProbability") {
  if (!weight_col %in% names(GenomicRanges::mcols(variant_GRanges))) {
    stop("variant_GRanges is missing weight column: ", weight_col)
  }

  peak_variant_overlaps_hits <- GenomicRanges::findOverlaps(
    query = peak_ranges,
    subject = variant_GRanges
  )

  out <- numeric(length(peak_ranges))
  names(out) <- get_peak_names_from_GRanges(peak_ranges)
  if (length(peak_variant_overlaps_hits) == 0) {
    return(out)
  }

  peak_idx <- S4Vectors::queryHits(peak_variant_overlaps_hits)
  variant_idx <- S4Vectors::subjectHits(peak_variant_overlaps_hits)
  variant_weights <- GenomicRanges::mcols(variant_GRanges)[[weight_col]][variant_idx]
  variant_weights[is.na(variant_weights)] <- 0

  out <- Matrix::sparseMatrix(
    i = peak_idx,
    j = rep.int(1L, length(peak_idx)),
    x = variant_weights,
    dims = c(length(peak_ranges), 1)
  ) |>
    Matrix::rowSums()
  names(out) <- get_peak_names_from_GRanges(peak_ranges)
  out
}

#' Get summed posterior probabilities per peak
#'
#' Convert credible-set variant PIPs into peak-level GWAS weights.
#'
#' @param credible_set_GRanges GRanges object containing credible set GRanges coordinates and metadata.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param GWAS_ID Configured GWAS label used in target names, plots, and Open Targets joins.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param weight_transform Weight post-processing mode. `cap_1` caps summed peak
#'   weights at 1; `sum` leaves summed weights unchanged.
#' @param ... Additional arguments forwarded to `filter_credible_set_variants()`.
#' @return Named numeric vector of peak weights, with names matching peak ranges.
#'   Errors if no credible-set variants overlap peaks for the GWAS.
#' @keywords internal

get_summed_posterior_probabilities_per_peak <- function(
  credible_set_GRanges,
  peak_ranges,
  GWAS_ID,
  posterior_probability_cutoff = NULL,
  posterior_probability_weighting_function = NULL,
  weight_transform = "cap_1",
  ...
) {
  credible_set_GRanges <- credible_set_GRanges |>
    filter_credible_set_variants(
      posterior_probability_cutoff = posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function,
      ...
    )

  posterior_probability_sums_per_peak <- sum_variant_weights_in_peaks(
    variant_GRanges = credible_set_GRanges,
    peak_ranges = peak_ranges,
    weight_col = "posteriorProbability"
  )

  if (sum(posterior_probability_sums_per_peak > 0) == 0) {
    stop(sprintf("No peaks overlap with credible-set variants for GWAS trait: %s - skipped.", GWAS_ID))
  }

  if (identical(weight_transform, "cap_1")) {
    n_capped <- sum(posterior_probability_sums_per_peak > 1)
    if (n_capped > 0) {
      message(sprintf(
        "Capped %s GWAS_chromVAR peak weight(s) above 1 for %s; max uncapped weight was %.3f.",
        n_capped,
        GWAS_ID,
        max(posterior_probability_sums_per_peak)
      ))
    }
    posterior_probability_sums_per_peak <- pmin(posterior_probability_sums_per_peak, 1)
  } else if (!identical(weight_transform, "sum")) {
    stop("Unsupported weight_transform: ", weight_transform)
  }

  posterior_probability_sums_per_peak
}

#' Get GWAS chromVAR peak weight record
#'
#' Build one peak-weight record for GWAS chromVAR scoring.
#'
#' @param GWAS_input_record Single GWAS branch record containing the study, finemapping method, and weighting mode.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param weight_transform Optional function or scalar transform applied to variant weights before aggregation.
#' @param ... Additional arguments forwarded to peak-weight construction helpers.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_GWAS_chromVAR_peak_weight_record <- function(
  GWAS_input_record,
  peak_ranges,
  posterior_probability_cutoff = NULL,
  posterior_probability_weighting_function = NULL,
  weight_transform = "cap_1",
  ...
) {
  GWAS_ID <- GWAS_input_record$GWAS_ID
  list(
    GWAS_ID = GWAS_ID,
    peak_weights_vec = get_summed_posterior_probabilities_per_peak(
      credible_set_GRanges = GWAS_input_record$credible_set_GRanges,
      peak_ranges = peak_ranges,
      GWAS_ID = GWAS_ID,
      posterior_probability_cutoff = posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function,
      weight_transform = weight_transform,
      ...
    )
  )
}

get_GWAS_chromVAR_peak_weight_summary_tibble <- function(peak_weight_records) {
  total_peaks <- length(peak_weight_records[[1]]$peak_weights_vec)

  peak_weight_records |>
    purrr::map_dfr(\(peak_weight_record) {
      tibble::tibble(
        GWAS_ID = peak_weight_record$GWAS_ID,
        total_capped_posteriorProbability = sum(peak_weight_record$peak_weights_vec),
        frac_overlapped_peaks = sum(peak_weight_record$peak_weights_vec > 0) / total_peaks
      )
    })
}

get_GWAS_chromVAR_peak_weight_matrix <- function(peak_weight_records, RSE_ATAC) {
  peak_weights_matrix <- peak_weight_records |>
    purrr::map("peak_weights_vec") |>
    do.call(what = cbind)
  colnames(peak_weights_matrix) <- purrr::map_chr(peak_weight_records, "GWAS_ID")

  peak_weights_matrix |>
    align_peak_weights_to_RSE(RSE_ATAC = RSE_ATAC) |>
    Matrix::Matrix(sparse = TRUE)
}

align_peak_weights_to_RSE <- function(peak_weights_matrix, RSE_ATAC) {
  peak_names <- get_peak_names_from_GRanges(SummarizedExperiment::rowRanges(RSE_ATAC))
  if (is.null(rownames(peak_weights_matrix))) {
    if (nrow(peak_weights_matrix) != length(peak_names)) {
      stop("Unnamed peak_weights_matrix must have the same number of rows as RSE_ATAC peaks.")
    }
    rownames(peak_weights_matrix) <- peak_names
  }

  missing_peaks <- setdiff(peak_names, rownames(peak_weights_matrix))
  if (length(missing_peaks) > 0) {
    stop("Missing peak weight rows for ", length(missing_peaks), " ATAC peak(s).")
  }

  peak_weights_matrix[peak_names, , drop = FALSE]
}

#' Get GWAS chromVAR z-score chunk record
#'
#' Compute one GWAS chromVAR z-score vector for one cell chunk.
#'
#' @param peak_weight_record One GWAS record containing `GWAS_ID` and named
#'   peak weights.
#' @param chunk_context_record Cell-chunk record containing counts, background,
#'   chunk ID, and cell names.
#' @param RSE_ATAC RangedSummarizedExperiment for ATAC peaks, with row ranges aligned to peak-level matrices.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_GWAS_chromVAR_z_score_chunk_record <- function(peak_weight_record, chunk_context_record, RSE_ATAC) {
  peak_weights_matrix <- matrix(peak_weight_record$peak_weights_vec, ncol = 1)
  rownames(peak_weights_matrix) <- names(peak_weight_record$peak_weights_vec)
  colnames(peak_weights_matrix) <- peak_weight_record$GWAS_ID
  peak_weights_matrix <- align_peak_weights_to_RSE(
    peak_weights_matrix = peak_weights_matrix,
    RSE_ATAC = RSE_ATAC
  )
  peak_weights_matrix <- Matrix::Matrix(peak_weights_matrix, sparse = TRUE)

  chunk_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = chunk_context_record$counts),
    rowRanges = SummarizedExperiment::rowRanges(RSE_ATAC)
  )
  SummarizedExperiment::rowData(chunk_obj) <- SummarizedExperiment::rowData(RSE_ATAC)

  z_score_matrix <- betterChromVAR::computeDeviationsAnalytic(
    object = chunk_obj,
    background = chunk_context_record$background,
    annotations = peak_weights_matrix,
    verbose = FALSE,
    retSE = FALSE,
    compute = "z"
  )$z

  z_score_vec <- z_score_matrix[peak_weight_record$GWAS_ID, ]
  names(z_score_vec) <- chunk_context_record$cell_names
  tibble::tibble(
    GWAS_ID = peak_weight_record$GWAS_ID,
    chunk_id = chunk_context_record$chunk_id,
    z_score_vec = list(z_score_vec)
  )
}

combine_GWAS_chromVAR_z_score_chunk_records <- function(chromVAR_z_score_chunk_records) {
  chromVAR_z_score_chunk_records <- dplyr::arrange(chromVAR_z_score_chunk_records, chunk_id)
  list(
    GWAS_ID = chromVAR_z_score_chunk_records$GWAS_ID[[1]],
    z_score_vec = unlist(chromVAR_z_score_chunk_records$z_score_vec, use.names = TRUE)
  )
}

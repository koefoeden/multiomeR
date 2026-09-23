#' Map credible-set variant weights to ATAC peaks
#'
#' Sum the weights of the variants overlapping each peak, cap each peak weight
#' at 1 for the betterChromVAR continuous-annotation contract, and allocate the
#' capped weight back to the variants in proportion to their weights.
#'
#' @param GWAS_ID Configured GWAS label.
#' @param variant_GRanges Retained credible-set variants.
#' @param peak_ranges Ordered ATAC peak ranges.
#' @param weight_col Variant weight column.
#' @return One row per peak-variant overlap.
#' @keywords internal

get_GWAS_chromVAR_peak_variant_weight_tibble <- function(
  GWAS_ID,
  variant_GRanges,
  peak_ranges,
  weight_col = "posteriorProbability"
) {
  overlap_hits <- GenomicRanges::findOverlaps(peak_ranges, variant_GRanges)
  if (length(overlap_hits) == 0L) {
    stop("No credible-set variants overlap ATAC peaks for GWAS_ID: ", GWAS_ID)
  }

  peak_idx <- S4Vectors::queryHits(overlap_hits)
  variant_tibble <- GenomicRanges::mcols(variant_GRanges)[S4Vectors::subjectHits(overlap_hits), , drop = FALSE] |>
    as.data.frame() |>
    tibble::as_tibble()

  tibble::tibble(
    peak_name = get_peak_names_from_GRanges(peak_ranges)[peak_idx],
    peak_chromosome = as.character(GenomicRanges::seqnames(peak_ranges))[peak_idx],
    peak_start = GenomicRanges::start(peak_ranges)[peak_idx],
    peak_end = GenomicRanges::end(peak_ranges)[peak_idx]
  ) |>
    dplyr::bind_cols(variant_tibble) |>
    dplyr::mutate(
      GWAS_ID = .env$GWAS_ID,
      uncapped_peak_weight = sum(.data[[weight_col]]),
      peak_weight = pmin(.data$uncapped_peak_weight, 1),
      peak_weight_scale = .data$peak_weight / .data$uncapped_peak_weight,
      peak_variant_weight = .data[[weight_col]] * .data$peak_weight_scale,
      .by = peak_name
    ) |>
    dplyr::relocate(
      GWAS_ID,
      peak_name,
      peak_chromosome,
      peak_start,
      peak_end,
      peak_weight,
      uncapped_peak_weight,
      peak_variant_weight,
      peak_weight_scale
    )
}

#' Get GWAS chromVAR peak weight record
#'
#' @param peak_variant_weight_tibble Peak-variant weights of one GWAS.
#' @param peak_ranges Ordered ATAC peak ranges.
#' @return List with `GWAS_ID` and the capped weight of every peak, named by peak.
#' @keywords internal

get_GWAS_chromVAR_peak_weight_record <- function(peak_variant_weight_tibble, peak_ranges) {
  peak_weight_tibble <- dplyr::distinct(peak_variant_weight_tibble, peak_name, peak_weight)
  peak_weights_vec <- stats::setNames(numeric(length(peak_ranges)), get_peak_names_from_GRanges(peak_ranges))
  peak_weights_vec[peak_weight_tibble$peak_name] <- peak_weight_tibble$peak_weight
  list(GWAS_ID = peak_variant_weight_tibble$GWAS_ID[[1]], peak_weights_vec = peak_weights_vec)
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

# Peak weights are named from the rowRanges of the chromVAR object they annotate.
get_GWAS_chromVAR_peak_weight_matrix <- function(peak_weight_records) {
  peak_weights_matrix <- peak_weight_records |>
    purrr::map("peak_weights_vec") |>
    do.call(what = cbind)
  colnames(peak_weights_matrix) <- purrr::map_chr(peak_weight_records, "GWAS_ID")
  Matrix::Matrix(peak_weights_matrix, sparse = TRUE)
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
  peak_weights_matrix <- Matrix::Matrix(peak_weights_matrix, sparse = TRUE)

  chunk_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = chunk_context_record$counts),
    rowRanges = SummarizedExperiment::rowRanges(RSE_ATAC)
  )

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

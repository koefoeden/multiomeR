#' Infer automatic absolute-effect weighting support for one GWAS
#'
#' Prefer variant-level betas when every credible-set variant has one. Otherwise,
#' use one beta per locus when every locus has at least one usable beta. Inputs
#' with incomplete locus coverage are skipped rather than partially represented.
#'
#' @param GWAS_input_record One normalized GWAS input record.
#' @param posterior_probability_cutoff Raw PIP threshold applied before effect
#'   weighting.
#' @return One-row tibble describing eligibility, route, and beta coverage.
#' @keywords internal

infer_GWAS_absolute_effect_weighting <- function(
  GWAS_input_record,
  posterior_probability_cutoff = NULL
) {
  GWAS_ID <- GWAS_input_record$GWAS_ID
  variant_tibble <- S4Vectors::mcols(GWAS_input_record$credible_set_GRanges) |>
    as.data.frame() |>
    tibble::as_tibble()
  required_columns <- c("studyLocusId", "posteriorProbability", "beta")
  missing_columns <- setdiff(required_columns, names(variant_tibble))

  if (length(missing_columns) > 0L) {
    return(tibble::tibble(
      GWAS_ID = GWAS_ID,
      studyId = GWAS_input_record$studyId,
      finemappingMethod = GWAS_input_record$finemappingMethod,
      eligible = FALSE,
      effect_weighting_route = NA_character_,
      status = paste0("missing_columns:", paste(missing_columns, collapse = ",")),
      n_variants = nrow(variant_tibble),
      n_variants_with_beta = 0L,
      n_loci = dplyr::n_distinct(variant_tibble$studyLocusId),
      n_loci_with_beta = 0L
    ))
  }

  raw_PIP <- if ("posteriorProbability_raw" %in% names(variant_tibble)) {
    variant_tibble$posteriorProbability_raw
  } else {
    variant_tibble$posteriorProbability
  }
  finite_beta <- is.finite(variant_tibble$beta)
  retained <- is.finite(raw_PIP) & raw_PIP >= 0
  if (!is.null(posterior_probability_cutoff)) {
    retained <- retained & raw_PIP > posterior_probability_cutoff
  }
  locus_has_beta <- tapply(
    finite_beta,
    variant_tibble$studyLocusId,
    any
  )
  all_variants_have_beta <- all(finite_beta)
  all_loci_have_beta <- all(locus_has_beta)

  route <- dplyr::case_when(
    all_variants_have_beta ~ "variant_absolute_effect",
    all_loci_have_beta ~ "locus_absolute_effect",
    TRUE ~ NA_character_
  )
  positive_effect_available <- any(
    retained & finite_beta & abs(variant_tibble$beta) > 0
  )
  eligible <- !is.na(route) && positive_effect_available
  status <- dplyr::case_when(
    !all_loci_have_beta ~ "skipped_incomplete_locus_beta_coverage",
    !positive_effect_available ~ "skipped_no_positive_retained_effect",
    TRUE ~ "eligible"
  )

  tibble::tibble(
    GWAS_ID = GWAS_ID,
    studyId = GWAS_input_record$studyId,
    finemappingMethod = GWAS_input_record$finemappingMethod,
    eligible = eligible,
    effect_weighting_route = route,
    status = status,
    n_variants = nrow(variant_tibble),
    n_variants_with_beta = sum(finite_beta),
    n_loci = length(locus_has_beta),
    n_loci_with_beta = sum(locus_has_beta)
  )
}

#' Summarize automatically inferred absolute-effect weighting routes
#'
#' @param GWAS_input_records List of normalized GWAS input records.
#' @param posterior_probability_cutoff Raw PIP threshold applied before effect
#'   weighting.
#' @return One row per raw-PIP GWAS input with eligibility and coverage.
#' @keywords internal

get_GWAS_absolute_effect_weighting_status_tibble <- function(
  GWAS_input_records,
  posterior_probability_cutoff = NULL
) {
  GWAS_input_records |>
    purrr::keep(\(record) identical(record$variant_weighting_mode, "raw_PIP")) |>
    purrr::map_dfr(
      infer_GWAS_absolute_effect_weighting,
      posterior_probability_cutoff = posterior_probability_cutoff
    )
}

#' Build absolute-effect-weighted credible-set ranges
#'
#' @param GWAS_input_record One normalized GWAS input record.
#' @param effect_weighting_route Inferred variant- or locus-level beta route.
#' @param posterior_probability_cutoff Raw PIP threshold applied before effect
#'   weighting.
#' @return Credible-set ranges carrying PIP-mass-preserving
#'   `absolute_effect_weight`.
#' @keywords internal

get_GWAS_absolute_effect_variant_GRanges <- function(
  GWAS_input_record,
  effect_weighting_route,
  posterior_probability_cutoff = NULL
) {
  variant_GRanges <- GWAS_input_record$credible_set_GRanges
  variant_tibble <- S4Vectors::mcols(variant_GRanges) |>
    as.data.frame() |>
    tibble::as_tibble()
  raw_PIP <- if ("posteriorProbability_raw" %in% names(variant_tibble)) {
    variant_tibble$posteriorProbability_raw
  } else {
    variant_tibble$posteriorProbability
  }

  effect_beta <- if (identical(effect_weighting_route, "variant_absolute_effect")) {
    variant_tibble$beta
  } else if (identical(effect_weighting_route, "locus_absolute_effect")) {
    locus_beta_tibble <- variant_tibble |>
      dplyr::mutate(raw_PIP = raw_PIP) |>
      dplyr::filter(is.finite(.data$beta)) |>
      dplyr::arrange(dplyr::desc(.data$raw_PIP)) |>
      dplyr::slice_head(n = 1L, by = "studyLocusId") |>
      dplyr::select(studyLocusId, effect_beta = beta)
    locus_beta_tibble$effect_beta[
      match(variant_tibble$studyLocusId, locus_beta_tibble$studyLocusId)
    ]
  } else {
    stop("Unsupported effect_weighting_route: ", effect_weighting_route)
  }

  keep <- is.finite(raw_PIP) & raw_PIP >= 0 & is.finite(effect_beta)
  if (!is.null(posterior_probability_cutoff)) {
    keep <- keep & raw_PIP > posterior_probability_cutoff
  }
  absolute_effect_weight <- raw_PIP * abs(effect_beta)
  keep <- keep & absolute_effect_weight > 0
  absolute_effect_weight <- absolute_effect_weight[keep]
  if (length(absolute_effect_weight) == 0L) {
    stop("No positive absolute-effect variant weights for GWAS_ID: ", GWAS_input_record$GWAS_ID)
  }
  retained_raw_PIP <- raw_PIP[keep]
  absolute_effect_scale <-
    sum(absolute_effect_weight) / sum(retained_raw_PIP)

  variant_GRanges <- variant_GRanges[keep]
  S4Vectors::mcols(variant_GRanges)$posteriorProbability_raw <- retained_raw_PIP
  S4Vectors::mcols(variant_GRanges)$effect_beta <- effect_beta[keep]
  S4Vectors::mcols(variant_GRanges)$absolute_effect_scale <-
    absolute_effect_scale
  S4Vectors::mcols(variant_GRanges)$absolute_effect_weight <-
    absolute_effect_weight / absolute_effect_scale
  variant_GRanges
}

#' Build peak weights for every automatically eligible GWAS
#'
#' @param GWAS_input_records List of normalized GWAS input records.
#' @param weighting_status_tibble Output of
#'   `get_GWAS_absolute_effect_weighting_status_tibble()`.
#' @param peak_ranges Ordered ATAC peak ranges.
#' @param posterior_probability_cutoff Raw PIP threshold applied before effect
#'   weighting.
#' @return List of eligible GWAS peak-weight records. Variant weights preserve
#'   retained raw-PIP mass before peak aggregation; summed peak annotations are
#'   capped at 1 for the betterChromVAR continuous-annotation contract.
#' @keywords internal

get_GWAS_absolute_effect_peak_weight_records <- function(
  GWAS_input_records,
  weighting_status_tibble,
  peak_ranges,
  posterior_probability_cutoff = NULL
) {
  eligible_status <- weighting_status_tibble |>
    dplyr::filter(.data$eligible)
  input_records <- GWAS_input_records |>
    purrr::keep(\(record) record$GWAS_ID %in% eligible_status$GWAS_ID)

  purrr::map(input_records, \(record) {
    route <- eligible_status$effect_weighting_route[
      match(record$GWAS_ID, eligible_status$GWAS_ID)
    ]
    variant_GRanges <- get_GWAS_absolute_effect_variant_GRanges(
      GWAS_input_record = record,
      effect_weighting_route = route,
      posterior_probability_cutoff = posterior_probability_cutoff
    )
    uncapped_peak_weights_vec <- sum_variant_weights_in_peaks(
      variant_GRanges = variant_GRanges,
      peak_ranges = peak_ranges,
      weight_col = "absolute_effect_weight"
    )
    peak_weights_vec <- pmin(uncapped_peak_weights_vec, 1)
    if (!any(peak_weights_vec > 0)) {
      stop("No absolute-effect-weighted variants overlap ATAC peaks for GWAS_ID: ", record$GWAS_ID)
    }
    list(
      GWAS_ID = record$GWAS_ID,
      effect_weighting_route = route,
      n_capped_peaks = sum(uncapped_peak_weights_vec > 1),
      peak_weights_vec = peak_weights_vec
    )
  })
}

#' Compute the absolute-effect-weighted cell-type chromVAR heatmap table
#'
#' @param psbulk_ATAC_data_matrix Peak-by-cell-type pseudobulk counts.
#' @param chromVAR_obj Template chromVAR object with aligned peak ranges.
#' @param annotation_matrix Peak-by-GWAS absolute-effect weight matrix.
#' @param background_record Precomputed pseudobulk betterChromVAR background.
#' @param GWAS_inputs_tibble Metadata for eligible weighted GWAS inputs.
#' @param cell_type_support_tibble Optional nuclei and ATAC-depth support.
#' @return One row per eligible GWAS and cell type with deviation, relative
#'   deviation, z score, one-sided P value, BH-adjusted P value, and labels.
#' @keywords internal

get_GWAS_absolute_effect_chromVAR_deviation_tibble <- function(
  psbulk_ATAC_data_matrix,
  chromVAR_obj,
  annotation_matrix,
  background_record,
  GWAS_inputs_tibble,
  cell_type_support_tibble = NULL
) {
  GWAS_IDs <- GWAS_inputs_tibble$GWAS_ID
  annotation_matrix <- annotation_matrix[, GWAS_IDs, drop = FALSE]
  chromVAR_dev <- compute_pseudobulk_chromVAR_deviation_SE(
    psbulk_ATAC_data_matrix = psbulk_ATAC_data_matrix,
    chromVAR_obj = chromVAR_obj,
    annotation_matrix = annotation_matrix,
    background_record = background_record
  )
  deviation_matrix <- SummarizedExperiment::assay(chromVAR_dev, "deviations")
  z_matrix <- SummarizedExperiment::assay(chromVAR_dev, "z")

  deviation_tibble <- as.data.frame(
    as.table(deviation_matrix),
    responseName = "deviation"
  ) |>
    tibble::as_tibble() |>
    dplyr::rename(GWAS_ID = Var1, cluster = Var2) |>
    dplyr::mutate(
      GWAS_ID = as.character(.data$GWAS_ID),
      cluster = as.character(.data$cluster)
    )
  z_tibble <- as.data.frame(as.table(z_matrix), responseName = "z") |>
    tibble::as_tibble() |>
    dplyr::rename(GWAS_ID = Var1, cluster = Var2) |>
    dplyr::mutate(
      GWAS_ID = as.character(.data$GWAS_ID),
      cluster = as.character(.data$cluster)
    )

  out <- deviation_tibble |>
    dplyr::left_join(
      z_tibble,
      by = c("GWAS_ID", "cluster"),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      relative_deviation = scale_within_vector(.data$deviation),
      .by = GWAS_ID
    ) |>
    dplyr::left_join(
      GWAS_inputs_tibble |>
        dplyr::select(
          GWAS_ID,
          Category,
          variant_weighting_mode,
          effect_weighting_route
        ),
      by = "GWAS_ID",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      z_p = stats::pnorm(.data$z, lower.tail = FALSE),
      z_q = stats::p.adjust(.data$z_p, method = "BH"),
      support_label = dplyr::case_when(
        .data$z >= 3 ~ "**",
        .data$z >= 2 ~ "*",
        .default = ""
      )
    )

  if (!is.null(cell_type_support_tibble)) {
    out <- out |>
      dplyr::left_join(cell_type_support_tibble, by = "cluster")
  }

  out |>
    dplyr::relocate(
      GWAS_ID,
      Category,
      variant_weighting_mode,
      effect_weighting_route,
      cluster,
      deviation,
      relative_deviation,
      z,
      z_p,
      z_q,
      support_label
    )
}

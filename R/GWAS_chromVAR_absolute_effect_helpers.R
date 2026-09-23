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

infer_GWAS_absolute_effect_weighting <- function(GWAS_input_record, posterior_probability_cutoff) {
  variant_tibble <- S4Vectors::mcols(GWAS_input_record$credible_set_GRanges) |>
    as.data.frame() |>
    tibble::as_tibble()

  raw_PIP <- variant_tibble$posteriorProbability
  finite_beta <- is.finite(variant_tibble$beta)
  retained <- is.finite(raw_PIP) & raw_PIP > posterior_probability_cutoff
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
    GWAS_ID = GWAS_input_record$GWAS_ID,
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
#' @return One row per GWAS input with eligibility and coverage.
#' @keywords internal

get_GWAS_absolute_effect_weighting_status_tibble <- function(GWAS_input_records, posterior_probability_cutoff) {
  GWAS_input_records |>
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
  posterior_probability_cutoff
) {
  variant_GRanges <- GWAS_input_record$credible_set_GRanges
  variant_tibble <- S4Vectors::mcols(variant_GRanges) |>
    as.data.frame() |>
    tibble::as_tibble()
  raw_PIP <- variant_tibble$posteriorProbability

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

  keep <- is.finite(raw_PIP) & raw_PIP > posterior_probability_cutoff & is.finite(effect_beta)
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
  S4Vectors::mcols(variant_GRanges)$effect_beta <- effect_beta[keep]
  S4Vectors::mcols(variant_GRanges)$absolute_effect_scale <-
    absolute_effect_scale
  S4Vectors::mcols(variant_GRanges)$absolute_effect_weight <-
    absolute_effect_weight / absolute_effect_scale
  variant_GRanges
}

#' Map absolute-effect variant weights to ATAC peaks for every eligible GWAS
#'
#' @param GWAS_input_records List of normalized GWAS input records.
#' @param weighting_status_tibble Output of
#'   `get_GWAS_absolute_effect_weighting_status_tibble()`.
#' @param peak_ranges Ordered ATAC peak ranges.
#' @param posterior_probability_cutoff Raw PIP threshold applied before effect
#'   weighting.
#' @return One peak-variant weight tibble per eligible GWAS. Variant weights
#'   preserve retained raw-PIP mass before peak aggregation.
#' @keywords internal

get_GWAS_absolute_effect_peak_variant_weight_records <- function(
  GWAS_input_records,
  weighting_status_tibble,
  peak_ranges,
  posterior_probability_cutoff
) {
  eligible <- dplyr::filter(weighting_status_tibble, .data$eligible)
  GWAS_input_records |>
    purrr::keep(\(record) record$GWAS_ID %in% eligible$GWAS_ID) |>
    purrr::map(\(record) {
      get_GWAS_chromVAR_peak_variant_weight_tibble(
        GWAS_ID = record$GWAS_ID,
        variant_GRanges = get_GWAS_absolute_effect_variant_GRanges(
          GWAS_input_record = record,
          effect_weighting_route = eligible$effect_weighting_route[match(record$GWAS_ID, eligible$GWAS_ID)],
          posterior_probability_cutoff = posterior_probability_cutoff
        ),
        peak_ranges = peak_ranges,
        weight_col = "absolute_effect_weight"
      )
    })
}

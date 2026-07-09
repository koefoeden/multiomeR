#' Normalize open targets variant weighting mode
#'
#' Normalize configured Open Targets variant-weighting mode values.
#'
#' @param mode Configured Open Targets variant weighting mode; supported values normalize to `raw_PIP` or `locus_lead_effect`.
#' @param GWAS_ID Configured GWAS label used in target names, plots, and Open Targets joins.
#' @return One of `raw_PIP` or `locus_lead_effect`; unsupported non-empty values
#'   error with the GWAS ID for diagnosis.
#' @keywords internal

normalize_open_targets_variant_weighting_mode <- function(mode, GWAS_ID = NA_character_) {
  if (is.null(mode) || length(mode) == 0 || isFALSE(mode)) {
    return("raw_PIP")
  }
  if (length(mode) != 1L) {
    stop("variant_weighting_mode must be length 1 for GWAS_ID: ", GWAS_ID)
  }
  if (is.na(mode)) {
    return("raw_PIP")
  }

  mode_chr <- stringr::str_to_lower(stringr::str_trim(as.character(mode)))
  if (mode_chr %in% c("", "null", "false", "none", "raw_pip")) {
    return("raw_PIP")
  }
  if (identical(mode_chr, "locus_lead_effect")) {
    return("locus_lead_effect")
  }

  stop(
    "Unsupported variant_weighting_mode for GWAS_ID=",
    GWAS_ID,
    ": ",
    mode,
    ". Allowed values are NULL, FALSE, '', raw_PIP, none, and locus_lead_effect."
  )
}

#' Get open targets locus effect record
#'
#' Derive a locus-level effect multiplier from the lead credible-set variant.
#'
#' @param locus_tibble One Open Targets credible-set locus tibble with PIP, variant, and effect-size fields.
#' @return One-row tibble with GWAS ID, study locus ID, raw effect multiplier,
#'   effect source, and lead variant ID.
#' @keywords internal

get_open_targets_locus_effect_record <- function(locus_tibble) {
  GWAS_ID <- unique(locus_tibble$GWAS_ID)
  studyLocusId <- unique(locus_tibble$studyLocusId)
  finemappingMethod <- paste(sort(unique(locus_tibble$finemappingMethod)), collapse = ", ")
  if (length(GWAS_ID) != 1L || length(studyLocusId) != 1L) {
    stop("Expected one GWAS_ID and one studyLocusId per Open Targets locus weighting group.")
  }

  if (identical(locus_tibble$finemappingMethod[[1]], "SuSiE-inf")) {
    lead_row <- locus_tibble |>
      dplyr::filter(!is.na(.data$beta)) |>
      dplyr::arrange(dplyr::desc(.data$posteriorProbability)) |>
      dplyr::slice_head(n = 1)

    if (nrow(lead_row) == 1L) {
      return(tibble::tibble(
        GWAS_ID = GWAS_ID,
        studyLocusId = studyLocusId,
        locus_effect_raw = abs(lead_row$beta[[1]]),
        locus_effect_source = "susie_highest_PIP_beta",
        lead_variantId = lead_row$variantId[[1]]
      ))
    }
  } else {
    lead_row <- locus_tibble |>
      dplyr::mutate(
        has_pvalue = !is.na(.data$pValueMantissa) & !is.na(.data$pValueExponent),
        has_beta_standardError = !is.na(.data$beta) & !is.na(.data$standardError),
        has_effect_evidence = .data$has_pvalue | .data$has_beta_standardError
      ) |>
      dplyr::filter(.data$has_effect_evidence) |>
      dplyr::arrange(dplyr::desc(.data$posteriorProbability)) |>
      dplyr::slice_head(n = 1)

    if (nrow(lead_row) == 1L) {
      has_pvalue <- !is.na(lead_row$pValueMantissa[[1]]) && !is.na(lead_row$pValueExponent[[1]])
      log_p_value <- if (has_pvalue && lead_row$pValueMantissa[[1]] > 0) {
        log(lead_row$pValueMantissa[[1]]) + lead_row$pValueExponent[[1]] * log(10)
      } else {
        NA_real_
      }
      z_abs <- if (!is.na(log_p_value) && log_p_value <= 0) {
        stats::qnorm(log_p_value - log(2), lower.tail = FALSE, log.p = TRUE)
      } else {
        NA_real_
      }

      if (!is.na(lead_row$beta[[1]]) && !is.na(lead_row$standardError[[1]])) {
        locus_effect_raw <- pmax(abs(lead_row$beta[[1]]) - lead_row$standardError[[1]], 0)
        locus_effect_source <- "beta_standardError"
      } else if (!is.na(z_abs) && !is.na(lead_row$standardError[[1]])) {
        locus_effect_raw <- pmax(z_abs * lead_row$standardError[[1]] - lead_row$standardError[[1]], 0)
        locus_effect_source <- "pvalue_standardError"
      } else if (!is.na(z_abs)) {
        locus_effect_raw <- z_abs
        locus_effect_source <- "pvalue_z"
      } else {
        locus_effect_raw <- NA_real_
        locus_effect_source <- NA_character_
      }

      if (!is.na(locus_effect_raw) && is.finite(locus_effect_raw)) {
        return(tibble::tibble(
          GWAS_ID = GWAS_ID,
          studyLocusId = studyLocusId,
          locus_effect_raw = locus_effect_raw,
          locus_effect_source = locus_effect_source,
          lead_variantId = lead_row$variantId[[1]]
        ))
      }
    }
  }

  stop(
    "Could not derive locus_lead_effect multiplier for GWAS_ID=",
    GWAS_ID,
    ", studyLocusId=",
    studyLocusId,
    ", finemappingMethod=",
    finemappingMethod,
    ". Missing usable beta, standardError, pValueMantissa/pValueExponent, or logBF evidence."
  )
}

#' Weight open targets posterior probability by locus effect
#'
#' Apply configured locus-effect scaling to Open Targets posterior probabilities.
#'
#' @param variants_tibble Open Targets variant tibble with posterior probabilities and optional locus-effect weighting fields.
#' @return Variant tibble with raw PIP preserved, weighting-mode metadata, locus
#'   multiplier fields, and updated non-negative `posteriorProbability`.
#' @keywords internal

weight_open_targets_posterior_probability_by_locus_effect <- function(variants_tibble) {
  allowed_modes <- c("raw_PIP", "locus_lead_effect")
  unknown_modes <- setdiff(unique(variants_tibble$variant_weighting_mode), allowed_modes)
  if (length(unknown_modes) > 0) {
    stop("Unsupported internal variant_weighting_mode value(s): ", paste(unknown_modes, collapse = ", "))
  }

  variants_tibble <- variants_tibble |>
    dplyr::mutate(
      posteriorProbability_raw = .data$posteriorProbability,
      posteriorProbability_weighting_mode = .data$variant_weighting_mode,
      locus_effect_multiplier = 1,
      locus_effect_source = dplyr::if_else(.data$variant_weighting_mode == "raw_PIP", "raw_PIP", NA_character_),
      lead_variantId = NA_character_
    )

  weighted_tibble <- variants_tibble |>
    dplyr::filter(.data$variant_weighting_mode == "locus_lead_effect")
  if (nrow(weighted_tibble) == 0) {
    return(variants_tibble)
  }

  locus_effect_tibble <- weighted_tibble |>
    dplyr::group_by(GWAS_ID, studyLocusId) |>
    dplyr::group_split() |>
    purrr::map_dfr(get_open_targets_locus_effect_record)

  multiplier_tibble <- locus_effect_tibble |>
    dplyr::mutate(
      positive_effect = is.finite(.data$locus_effect_raw) & .data$locus_effect_raw > 0,
      effect_cap = stats::quantile(.data$locus_effect_raw[.data$positive_effect], 0.95, na.rm = TRUE),
      locus_effect_capped = pmin(.data$locus_effect_raw, .data$effect_cap),
      positive_capped = is.finite(.data$locus_effect_capped) & .data$locus_effect_capped > 0,
      positive_median = stats::median(.data$locus_effect_capped[.data$positive_capped], na.rm = TRUE),
      .by = GWAS_ID
    )

  invalid_medians <- multiplier_tibble |>
    dplyr::filter(is.na(.data$positive_median) | !is.finite(.data$positive_median) | .data$positive_median <= 0) |>
    dplyr::distinct(GWAS_ID)
  if (nrow(invalid_medians) > 0) {
    stop(
      "Could not normalize locus_lead_effect multipliers for GWAS_ID(s): ",
      paste(invalid_medians$GWAS_ID, collapse = ", "),
      ". Positive finite locus effects had missing, zero, or non-finite median."
    )
  }

  multiplier_tibble <- multiplier_tibble |>
    dplyr::transmute(
      GWAS_ID = .data$GWAS_ID,
      studyLocusId = .data$studyLocusId,
      locus_effect_multiplier = .data$locus_effect_capped / .data$positive_median,
      locus_effect_source = .data$locus_effect_source,
      lead_variantId = .data$lead_variantId
    )

  out <- variants_tibble |>
    dplyr::select(-dplyr::all_of(c("locus_effect_multiplier", "locus_effect_source", "lead_variantId"))) |>
    dplyr::left_join(multiplier_tibble, by = c("GWAS_ID", "studyLocusId")) |>
    dplyr::mutate(
      locus_effect_multiplier = dplyr::coalesce(.data$locus_effect_multiplier, 1),
      locus_effect_source = dplyr::coalesce(.data$locus_effect_source, "raw_PIP"),
      posteriorProbability = .data$posteriorProbability_raw * .data$locus_effect_multiplier
    )

  if (any(is.na(out$posteriorProbability) | !is.finite(out$posteriorProbability) | out$posteriorProbability < 0)) {
    stop("Open Targets variant weighting produced missing, non-finite, or negative posteriorProbability values.")
  }

  out
}

#' Get open targets credible set variants tibble
#'
#' Extract configured Open Targets credible-set variants for GWAS-gene chromVAR.
#'
#' @param GWAS_inputs_tibble GWAS configuration tibble with labels, Open Targets study IDs, finemapping methods, and weighting modes.
#' @param open_targets_credible_set_dataset_path Local directory or dataset path
#'   for the Open Targets credible-set Parquet dataset.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_open_targets_credible_set_variants_tibble <- function(GWAS_inputs_tibble, open_targets_credible_set_dataset_path) {
  if (nrow(GWAS_inputs_tibble) == 0) {
    return(tibble::tibble(
      GWAS_ID = character(),
      open_targets_release = character(),
      studyId = character(),
      studyLocusId = character(),
      credibleSetIndex = integer(),
      finemappingMethod = character(),
      variant_weighting_mode = character(),
      confidence = character(),
      variantId = character(),
      chromosome = character(),
      position = integer(),
      posteriorProbability = numeric(),
      posteriorProbability_raw = numeric(),
      posteriorProbability_weighting_mode = character(),
      locus_effect_multiplier = numeric(),
      locus_effect_source = character(),
      lead_variantId = character(),
      logBF = numeric(),
      pValueMantissa = numeric(),
      pValueExponent = integer(),
      beta = numeric(),
      standardError = numeric(),
      r2Overall = numeric(),
      is95CredibleSet = logical(),
      is99CredibleSet = logical(),
      locusStart = integer(),
      locusEnd = integer()
    ))
  }

  study_ids <- unique(GWAS_inputs_tibble$studyId)
  credible_set_tibble <- arrow::open_dataset(open_targets_credible_set_dataset_path) |>
    dplyr::filter(studyType == "gwas", studyId %in% study_ids) |>
    dplyr::select(
      studyId,
      studyLocusId,
      credibleSetIndex,
      finemappingMethod,
      confidence,
      locusStart,
      locusEnd,
      locus
    ) |>
    dplyr::collect() |>
    dplyr::inner_join(
      GWAS_inputs_tibble |>
        dplyr::select(GWAS_ID, studyId, finemappingMethod, open_targets_release, variant_weighting_mode),
      by = c("studyId", "finemappingMethod"),
      relationship = "many-to-many"
    )

  if (nrow(credible_set_tibble) == 0) {
    stop("No Open Targets credible sets found for the configured GWAS input rows.")
  }

  variants_tibble <- credible_set_tibble |>
    tidyr::unnest(locus)
  variant_parts <- stringr::str_match(variants_tibble$variantId, "^(.+)_([0-9]+)_([^_]+)_([^_]+)$")

  variants_tibble <- variants_tibble |>
    dplyr::mutate(
      chromosome = variant_parts[, 2],
      position = suppressWarnings(as.integer(variant_parts[, 3]))
    ) |>
    dplyr::filter(stringr::str_detect(chromosome, "^([0-9]+|X|Y|MT)$")) |>
    dplyr::select(
      GWAS_ID,
      open_targets_release,
      studyId,
      studyLocusId,
      credibleSetIndex,
      finemappingMethod,
      variant_weighting_mode,
      confidence,
      variantId,
      chromosome,
      position,
      posteriorProbability,
      logBF,
      pValueMantissa,
      pValueExponent,
      beta,
      standardError,
      r2Overall,
      is95CredibleSet,
      is99CredibleSet,
      locusStart,
      locusEnd
    )

  if (any(is.na(variants_tibble$position))) {
    stop("Open Targets variantId coordinate parsing produced missing positions.")
  }
  if (any(is.na(variants_tibble$posteriorProbability))) {
    stop("Open Targets credible sets contain missing posteriorProbability values.")
  }
  if (any(!variants_tibble$is95CredibleSet)) {
    stop("Open Targets credible-set table contains variants outside the 95% credible set.")
  }

  variants_tibble <- weight_open_targets_posterior_probability_by_locus_effect(variants_tibble)

  duplicated_variants <- variants_tibble |>
    dplyr::count(GWAS_ID, studyLocusId, variantId) |>
    dplyr::filter(n > 1)
  if (nrow(duplicated_variants) > 0) {
    stop("Open Targets method-only filtering duplicated GWAS_ID/studyLocusId/variantId rows. Review confidence values before processing.")
  }

  variants_tibble
}

get_open_targets_credible_set_GRanges <- function(variants_tibble) {
  variants_tibble |>
    dplyr::mutate(seqnames = stringr::str_c("chr", chromosome), start = position, end = position) |>
    GenomicRanges::makeGRangesFromDataFrame(
      seqnames.field = "seqnames",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = TRUE
    ) |>
    unname()
}

get_GWAS_chromVAR_input_record <- function(GWAS_input_tibble, open_targets_credible_set_dataset_path) {
  if (nrow(GWAS_input_tibble) != 1) {
    stop("GWAS_input_tibble must contain exactly one GWAS row.")
  }

  variants_tibble <- get_open_targets_credible_set_variants_tibble(
    GWAS_inputs_tibble = GWAS_input_tibble,
    open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
  )

  list(
    GWAS_ID = GWAS_input_tibble$GWAS_ID[[1]],
    studyId = GWAS_input_tibble$studyId[[1]],
    finemappingMethod = GWAS_input_tibble$finemappingMethod[[1]],
    variant_weighting_mode = GWAS_input_tibble$variant_weighting_mode[[1]],
    credible_set_GRanges = get_open_targets_credible_set_GRanges(variants_tibble)
  )
}

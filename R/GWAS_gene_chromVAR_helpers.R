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

get_GWAS_gene_chromVAR_inputs_tibble <- function(GWAS_analysis_inputs_tibble, selected_GWAS_IDs) {
  if (is.null(selected_GWAS_IDs) || length(selected_GWAS_IDs) == 0) {
    return(GWAS_analysis_inputs_tibble[0, ])
  }

  missing_GWAS_IDs <- setdiff(selected_GWAS_IDs, GWAS_analysis_inputs_tibble$GWAS_ID)
  if (length(missing_GWAS_IDs) > 0) {
    stop("Configured gene-level GWAS chromVAR GWAS_ID(s) missing from GWAS_analysis_inputs_tibble: ", paste(missing_GWAS_IDs, collapse = ", "))
  }

  GWAS_analysis_inputs_tibble |>
    dplyr::filter(GWAS_ID %in% selected_GWAS_IDs) |>
    dplyr::arrange(match(GWAS_ID, selected_GWAS_IDs))
}

#' Get GWAS gene chromVAR L2G tibble
#'
#' Retrieve Open Targets locus-to-gene scores for selected credible-set loci.
#'
#' @param GWAS_gene_chromVAR_credible_set_variants_tibble Credible-set variant
#'   tibble containing selected `studyLocusId`, `studyId`, `GWAS_ID`, and release.
#' @param open_targets_gwas_credible_sets_evidence_dataset_path Local directory
#'   or dataset path for Open Targets GWAS credible-set evidence.
#' @return Distinct GWAS/study-locus/target rows with `L2G_score >= 0.05`.
#' @keywords internal

get_GWAS_gene_chromVAR_L2G_tibble <- function(GWAS_gene_chromVAR_credible_set_variants_tibble, open_targets_gwas_credible_sets_evidence_dataset_path) {
  if (nrow(GWAS_gene_chromVAR_credible_set_variants_tibble) == 0) {
    return(tibble::tibble(
      GWAS_ID = character(),
      studyId = character(),
      studyLocusId = character(),
      targetId = character(),
      L2G_score = numeric(),
      open_targets_release = character()
    ))
  }

  selected_loci_tibble <- GWAS_gene_chromVAR_credible_set_variants_tibble |>
    dplyr::distinct(GWAS_ID, studyId, studyLocusId, open_targets_release)

  evidence_tibble <- arrow::open_dataset(open_targets_gwas_credible_sets_evidence_dataset_path) |>
    dplyr::filter(studyLocusId %in% selected_loci_tibble$studyLocusId) |>
    dplyr::select(studyLocusId, targetId, resourceScore, score) |>
    dplyr::collect()

  if (nrow(evidence_tibble) == 0) {
    return(tibble::tibble(
      GWAS_ID = character(),
      studyId = character(),
      studyLocusId = character(),
      targetId = character(),
      L2G_score = numeric(),
      open_targets_release = character()
    ))
  }

  score_delta <- abs(evidence_tibble$score - evidence_tibble$resourceScore)
  if (any(score_delta > sqrt(.Machine$double.eps), na.rm = TRUE)) {
    stop("Open Targets evidence_gwas_credible_sets score and resourceScore are not equivalent.")
  }

  evidence_tibble |>
    dplyr::filter(!is.na(score), score >= 0.05) |>
    dplyr::transmute(studyLocusId, targetId, L2G_score = score) |>
    dplyr::inner_join(selected_loci_tibble, by = "studyLocusId") |>
    dplyr::select(GWAS_ID, studyId, studyLocusId, targetId, L2G_score, open_targets_release) |>
    dplyr::distinct()
}

#' Get empty GWAS gene chromVAR record
#'
#' Return an empty GWAS-gene chromVAR record with valid matrix/tibble schema.
#'
#' @param GWAS_ID Configured GWAS label used in target names, plots, and Open Targets joins.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param open_targets_release Open Targets release label propagated into empty
#'   feature metadata.
#' @param variant_weighting_mode Normalized variant weighting mode, currently `raw_PIP` or `locus_lead_effect`.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_empty_GWAS_gene_chromVAR_record <- function(GWAS_ID, peak_ranges, open_targets_release, variant_weighting_mode = "raw_PIP") {
  targetId <- character()
  peak_names <- get_peak_names_from_GRanges(peak_ranges)
  peak_gene_weight_matrix <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    dims = c(length(peak_names), 0),
    dimnames = list(peak_names, targetId)
  )

  list(
    GWAS_ID = GWAS_ID,
    peak_gene_weight_matrix = peak_gene_weight_matrix,
    gene_support_tibble = tibble::tibble(
      GWAS_ID = character(),
      variant_weighting_mode = character(),
      targetId = character(),
      n_loci = integer(),
      n_variants = integer(),
      n_peaks = integer(),
      sum_PIP_x_L2G = numeric(),
      sum_peak_gene_weight = numeric(),
      sum_peak_gene_weight_raw = numeric(),
      effect_weighting_lift = numeric(),
      max_L2G_score = numeric(),
      open_targets_release = character()
    )
  )
}

#' Get GWAS gene chromVAR peak weight record
#'
#' Build peak-by-gene weights from credible-set variants and L2G evidence.
#'
#' @param GWAS_input_tibble One configured GWAS row containing GWAS ID, Open
#'   Targets release, and variant-weighting mode.
#' @param GWAS_gene_chromVAR_credible_set_variants_tibble Variant-level
#'   credible-set table with PIP and locus metadata for the configured GWAS.
#' @param GWAS_gene_chromVAR_L2G_tibble Locus-to-gene evidence table with
#'   target IDs and L2G scores for the selected loci.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param ... Additional arguments passed to variant weighting or overlap helpers.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_GWAS_gene_chromVAR_peak_weight_record <- function(
  GWAS_input_tibble,
  GWAS_gene_chromVAR_credible_set_variants_tibble,
  GWAS_gene_chromVAR_L2G_tibble,
  peak_ranges,
  posterior_probability_cutoff = NULL,
  posterior_probability_weighting_function = NULL,
  ...
) {
  if (nrow(GWAS_input_tibble) != 1) {
    stop("GWAS_input_tibble must contain exactly one GWAS row.")
  }

  GWAS_ID <- GWAS_input_tibble$GWAS_ID[[1]]
  open_targets_release <- GWAS_input_tibble$open_targets_release[[1]]
  variant_weighting_mode <- GWAS_input_tibble$variant_weighting_mode[[1]]
  peak_names <- get_peak_names_from_GRanges(peak_ranges)

  variants_tibble <- GWAS_gene_chromVAR_credible_set_variants_tibble |>
    dplyr::filter(GWAS_ID == .env$GWAS_ID)
  L2G_tibble <- GWAS_gene_chromVAR_L2G_tibble |>
    dplyr::filter(GWAS_ID == .env$GWAS_ID)

  if (nrow(variants_tibble) == 0 || nrow(L2G_tibble) == 0) {
    return(get_empty_GWAS_gene_chromVAR_record(GWAS_ID, peak_ranges, open_targets_release, variant_weighting_mode))
  }

  credible_set_GRanges <- variants_tibble |>
    get_open_targets_credible_set_GRanges() |>
    filter_credible_set_variants(
      posterior_probability_cutoff = posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function,
      ...
    )

  peak_variant_overlaps_hits <- GenomicRanges::findOverlaps(
    query = peak_ranges,
    subject = credible_set_GRanges
  )
  if (length(peak_variant_overlaps_hits) == 0) {
    return(get_empty_GWAS_gene_chromVAR_record(GWAS_ID, peak_ranges, open_targets_release, variant_weighting_mode))
  }

  variant_metadata_tibble <- GenomicRanges::mcols(credible_set_GRanges) |>
    as.data.frame() |>
    tibble::as_tibble()
  if (!"posteriorProbability_raw" %in% names(variant_metadata_tibble)) {
    variant_metadata_tibble$posteriorProbability_raw <- variant_metadata_tibble$posteriorProbability
  }
  variant_metadata_tibble <- variant_metadata_tibble |>
    dplyr::select(studyLocusId, variantId, posteriorProbability, posteriorProbability_raw)

  peak_variant_hits_tibble <- tibble::tibble(
    peak_idx = S4Vectors::queryHits(peak_variant_overlaps_hits),
    variant_idx = S4Vectors::subjectHits(peak_variant_overlaps_hits)
  )
  peak_variant_tibble <- peak_variant_hits_tibble |>
    dplyr::bind_cols(variant_metadata_tibble[peak_variant_hits_tibble$variant_idx, ]) |>
    dplyr::mutate(peak_name = peak_names[peak_idx])

  peak_locus_tibble <- peak_variant_tibble |>
    dplyr::summarise(
      capped_peak_locus_PIP = pmin(sum(posteriorProbability, na.rm = TRUE), 1),
      capped_peak_locus_PIP_raw = pmin(sum(posteriorProbability_raw, na.rm = TRUE), 1),
      .by = c(peak_idx, peak_name, studyLocusId)
    ) |>
    dplyr::inner_join(
      L2G_tibble |> dplyr::select(studyLocusId, targetId, L2G_score),
      by = "studyLocusId",
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      peak_gene_weight = capped_peak_locus_PIP * L2G_score,
      peak_gene_weight_raw = capped_peak_locus_PIP_raw * L2G_score
    )

  if (nrow(peak_locus_tibble) == 0) {
    return(get_empty_GWAS_gene_chromVAR_record(GWAS_ID, peak_ranges, open_targets_release, variant_weighting_mode))
  }

  peak_gene_tibble <- peak_locus_tibble |>
    dplyr::summarise(
      peak_gene_weight = sum(peak_gene_weight, na.rm = TRUE),
      peak_gene_weight_raw = sum(peak_gene_weight_raw, na.rm = TRUE),
      .by = c(peak_idx, peak_name, targetId)
    )

  target_ids <- sort(unique(peak_gene_tibble$targetId))
  peak_gene_weight_matrix <- Matrix::sparseMatrix(
    i = peak_gene_tibble$peak_idx,
    j = match(peak_gene_tibble$targetId, target_ids),
    x = peak_gene_tibble$peak_gene_weight,
    dims = c(length(peak_names), length(target_ids)),
    dimnames = list(peak_names, target_ids)
  )

  variant_L2G_tibble <- peak_variant_tibble |>
    dplyr::inner_join(
      L2G_tibble |> dplyr::select(studyLocusId, targetId, L2G_score),
      by = "studyLocusId",
      relationship = "many-to-many"
    )

  gene_support_tibble <- peak_locus_tibble |>
    dplyr::summarise(
      n_loci = dplyr::n_distinct(studyLocusId),
      n_peaks = dplyr::n_distinct(peak_name),
      sum_peak_gene_weight = sum(peak_gene_weight, na.rm = TRUE),
      sum_peak_gene_weight_raw = sum(peak_gene_weight_raw, na.rm = TRUE),
      max_L2G_score = max(L2G_score, na.rm = TRUE),
      .by = targetId
    ) |>
    dplyr::left_join(
      variant_L2G_tibble |>
        dplyr::summarise(
          n_variants = dplyr::n_distinct(variantId),
          sum_PIP_x_L2G = sum(posteriorProbability * L2G_score, na.rm = TRUE),
          .by = targetId
        ),
      by = "targetId"
    ) |>
    dplyr::mutate(
      GWAS_ID = GWAS_ID,
      variant_weighting_mode = variant_weighting_mode,
      effect_weighting_lift = dplyr::if_else(
        sum_peak_gene_weight_raw > 0,
        sum_peak_gene_weight / sum_peak_gene_weight_raw,
        NA_real_
      ),
      open_targets_release = open_targets_release,
      .before = targetId
    ) |>
    dplyr::select(
      GWAS_ID,
      variant_weighting_mode,
      targetId,
      n_loci,
      n_variants,
      n_peaks,
      sum_PIP_x_L2G,
      sum_peak_gene_weight,
      sum_peak_gene_weight_raw,
      effect_weighting_lift,
      max_L2G_score,
      open_targets_release
    ) |>
    dplyr::arrange(GWAS_ID, targetId)

  list(
    GWAS_ID = GWAS_ID,
    peak_gene_weight_matrix = peak_gene_weight_matrix,
    gene_support_tibble = gene_support_tibble
  )
}

#' Get GWAS gene chromVAR feature metadata tibble
#'
#' Summarize retained GWAS-gene chromVAR features and support filters.
#'
#' @param GWAS_gene_chromVAR_peak_weight_records List of per-GWAS peak/gene
#'   weight records containing `gene_support_tibble`.
#' @param min_n_peaks Minimum number of weighted peaks required for a GWAS/model feature to be kept.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_GWAS_gene_chromVAR_feature_metadata_tibble <- function(GWAS_gene_chromVAR_peak_weight_records, min_n_peaks = 3) {
  if (length(GWAS_gene_chromVAR_peak_weight_records) == 0) {
    return(tibble::tibble(
      feature_id = character(),
      GWAS_ID = character(),
      variant_weighting_mode = character(),
      targetId = character(),
      n_loci = integer(),
      n_variants = integer(),
      n_peaks = integer(),
      sum_PIP_x_L2G = numeric(),
      sum_peak_gene_weight = numeric(),
      sum_peak_gene_weight_raw = numeric(),
      effect_weighting_lift = numeric(),
      max_L2G_score = numeric(),
      open_targets_release = character(),
      support_filter_pass = logical()
    ))
  }

  GWAS_gene_chromVAR_peak_weight_records |>
    purrr::map_dfr("gene_support_tibble") |>
    dplyr::mutate(
      feature_id = stringr::str_c(GWAS_ID, targetId, sep = "__"),
      support_filter_pass = n_loci >= 1 & n_peaks >= min_n_peaks & sum_PIP_x_L2G > 0,
      .before = GWAS_ID
    ) |>
    dplyr::arrange(GWAS_ID, targetId)
}

#' Get open targets target metadata tibble
#'
#' Retrieve Open Targets target labels and biotypes for Ensembl target IDs.
#'
#' @param target_ids Character vector of Open Targets target IDs to retrieve from the target dataset.
#' @param open_targets_target_dataset_path Local directory or dataset path for
#'   the Open Targets target metadata Parquet dataset.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_open_targets_target_metadata_tibble <- function(target_ids, open_targets_target_dataset_path) {
  target_ids <- unique(stringr::str_remove(target_ids, "\\.[0-9]+$"))
  if (length(target_ids) == 0) {
    return(tibble::tibble(
      targetId = character(),
      gene_symbol = character(),
      gene_name = character(),
      gene_biotype = character(),
      gene_label = character()
    ))
  }

  target_metadata_tibble <- arrow::open_dataset(open_targets_target_dataset_path) |>
    dplyr::filter(id %in% target_ids) |>
    dplyr::select(id, approvedSymbol, approvedName, biotype) |>
    dplyr::collect()

  duplicated_target_ids <- target_metadata_tibble |>
    dplyr::count(id) |>
    dplyr::filter(n > 1)
  if (nrow(duplicated_target_ids) > 0) {
    stop("Open Targets target dataset contains duplicated target IDs: ", paste(duplicated_target_ids$id, collapse = ", "))
  }

  missing_target_ids <- setdiff(target_ids, target_metadata_tibble$id)
  if (length(missing_target_ids) > 0) {
    stop("Open Targets target metadata is missing target IDs used by L2G evidence: ", paste(missing_target_ids, collapse = ", "))
  }

  target_metadata_tibble <- target_metadata_tibble |>
    dplyr::transmute(
      targetId = id,
      gene_symbol = dplyr::na_if(stringr::str_squish(approvedSymbol), ""),
      gene_name = dplyr::na_if(stringr::str_squish(approvedName), ""),
      gene_biotype = dplyr::na_if(stringr::str_squish(biotype), "")
    ) |>
    dplyr::mutate(gene_symbol = dplyr::coalesce(gene_symbol, targetId))

  duplicated_symbols <- target_metadata_tibble |>
    dplyr::count(gene_symbol) |>
    dplyr::filter(n > 1) |>
    dplyr::pull(gene_symbol)

  target_metadata_tibble |>
    dplyr::mutate(
      gene_label = dplyr::if_else(
        gene_symbol %in% duplicated_symbols,
        stringr::str_c(gene_symbol, " (", targetId, ")"),
        gene_symbol
      )
    ) |>
    dplyr::arrange(match(targetId, target_ids))
}

add_open_targets_target_metadata <- function(GWAS_gene_chromVAR_feature_metadata_tibble, open_targets_target_dataset_path) {
  target_metadata_tibble <- get_open_targets_target_metadata_tibble(
    target_ids = GWAS_gene_chromVAR_feature_metadata_tibble$targetId,
    open_targets_target_dataset_path = open_targets_target_dataset_path
  )

  GWAS_gene_chromVAR_feature_metadata_tibble |>
    dplyr::left_join(target_metadata_tibble, by = "targetId") |>
    dplyr::relocate(gene_label, gene_symbol, gene_name, gene_biotype, .after = targetId)
}

#' Get open targets target locus tibble
#'
#' Retrieve genomic coordinates for one Open Targets target ID.
#'
#' @param target_id Open Targets target ID, usually an Ensembl gene ID.
#' @param open_targets_target_dataset_path Local directory or dataset path for
#'   the Open Targets target metadata Parquet dataset.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_open_targets_target_locus_tibble <- function(target_id, open_targets_target_dataset_path) {
  target_locus_tibble <- arrow::open_dataset(open_targets_target_dataset_path) |>
    dplyr::filter(id == .env$target_id) |>
    dplyr::select(id, approvedSymbol, approvedName, biotype, genomicLocation, tss) |>
    dplyr::collect()

  if (nrow(target_locus_tibble) != 1L) {
    stop("Expected one Open Targets target metadata row for ", target_id, ", found ", nrow(target_locus_tibble), ".")
  }

  target_locus_tibble |>
    dplyr::transmute(
      targetId = id,
      gene_symbol = dplyr::na_if(stringr::str_squish(approvedSymbol), ""),
      gene_name = dplyr::na_if(stringr::str_squish(approvedName), ""),
      gene_biotype = dplyr::na_if(stringr::str_squish(biotype), ""),
      chr = stringr::str_c("chr", genomicLocation$chromosome),
      start = genomicLocation$start,
      end = genomicLocation$end,
      strand = genomicLocation$strand,
      tss = tss
    ) |>
    dplyr::mutate(gene_symbol = dplyr::coalesce(gene_symbol, targetId))
}

#' Get GWAS gene chromVAR locus check tibble
#'
#' Assemble variant, L2G, and target-window data for one locus-check plot.
#'
#' @param psbulk_GWAS_gene_chromVAR_results_tibble Formatted pseudobulk results
#'   containing GWAS, target, and contrast hits.
#' @param GWAS_gene_chromVAR_credible_set_variants_tibble Credible-set variants
#'   used to draw PIP tracks.
#' @param GWAS_gene_chromVAR_L2G_tibble Locus-to-gene score tibble used to
#'   identify linked loci for the target.
#' @param open_targets_target_dataset_path Local directory or dataset path for
#'   Open Targets target metadata.
#' @param GWAS_ID Configured GWAS label used in target names, plots, and Open Targets joins.
#' @param target_id Open Targets target ID, usually an Ensembl gene ID.
#' @param contrast Model contrast label used to filter differential results or annotate plots.
#' @param flank Number of base pairs added upstream and downstream when building locus-check plotting windows.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_GWAS_gene_chromVAR_locus_check_tibble <- function(
  psbulk_GWAS_gene_chromVAR_results_tibble,
  GWAS_gene_chromVAR_credible_set_variants_tibble,
  GWAS_gene_chromVAR_L2G_tibble,
  open_targets_target_dataset_path,
  GWAS_ID,
  target_id,
  contrast,
  model = NULL,
  flank = 50000L
) {
  result_tibble <- psbulk_GWAS_gene_chromVAR_results_tibble |>
    dplyr::filter(
      .data$GWAS_ID == .env$GWAS_ID,
      .data$targetId == .env$target_id,
      .data$contrast == .env$contrast
    )
  if (!is.null(model)) {
    result_tibble <- result_tibble |>
      dplyr::filter(.data$model == .env$model)
  }

  if (nrow(result_tibble) != 1L) {
    stop(
      "Expected one GWAS-gene chromVAR result row for ",
      GWAS_ID,
      " / ",
      target_id,
      " / ",
      contrast,
      if (!is.null(model)) paste0(" / ", model) else "",
      ", found ",
      nrow(result_tibble),
      "."
    )
  }

  target_locus_tibble <- get_open_targets_target_locus_tibble(
    target_id = target_id,
    open_targets_target_dataset_path = open_targets_target_dataset_path
  )

  L2G_tibble <- GWAS_gene_chromVAR_L2G_tibble |>
    dplyr::filter(.data$GWAS_ID == .env$GWAS_ID, .data$targetId == .env$target_id)
  if (nrow(L2G_tibble) == 0L) {
    stop("No Open Targets L2G rows found for ", GWAS_ID, " / ", target_id, ".")
  }

  variant_tibble <- GWAS_gene_chromVAR_credible_set_variants_tibble |>
    dplyr::semi_join(
      L2G_tibble |> dplyr::select(studyLocusId),
      by = "studyLocusId"
    ) |>
    dplyr::filter(.data$GWAS_ID == .env$GWAS_ID) |>
    dplyr::mutate(chr = stringr::str_c("chr", .data$chromosome))
  if (nrow(variant_tibble) == 0L) {
    stop("No Open Targets credible-set variants found for ", GWAS_ID, " / ", target_id, ".")
  }

  locus_start <- pmax(
    1L,
    min(c(target_locus_tibble$start, target_locus_tibble$end, variant_tibble$position), na.rm = TRUE) - flank
  )
  locus_end <- max(c(target_locus_tibble$start, target_locus_tibble$end, variant_tibble$position), na.rm = TRUE) + flank

  result_tibble |>
    dplyr::select(GWAS_ID, variant_weighting_mode, model, contrast, targetId, gene_label, gene_symbol, logFC, PValue, FDR, n_loci, n_variants, n_peaks, sum_PIP_x_L2G, max_L2G_score) |>
    dplyr::mutate(
      gene_name = target_locus_tibble$gene_name[[1]],
      gene_biotype = target_locus_tibble$gene_biotype[[1]],
      chr = target_locus_tibble$chr[[1]],
      target_start = target_locus_tibble$start[[1]],
      target_end = target_locus_tibble$end[[1]],
      target_tss = target_locus_tibble$tss[[1]],
      locus_start = as.integer(locus_start),
      locus_end = as.integer(locus_end),
      L2G_score = max(L2G_tibble$L2G_score, na.rm = TRUE)
    )
}

empty_GWAS_gene_chromVAR_locus_check_tibble <- function(psbulk_GWAS_gene_chromVAR_results_tibble) {
  psbulk_GWAS_gene_chromVAR_results_tibble[0, , drop = FALSE] |>
    dplyr::mutate(
      chr = character(),
      target_start = integer(),
      target_end = integer(),
      target_tss = integer(),
      locus_start = integer(),
      locus_end = integer(),
      L2G_score = numeric(),
      locus_plot_rank = integer(),
      plot_name = character()
    )
}

make_GWAS_gene_chromVAR_locus_plot_name <- function(rank, GWAS_ID, gene_symbol, targetId, model, contrast) {
  label <- stringr::str_c(
    sprintf("%02d", rank),
    GWAS_ID,
    dplyr::coalesce(gene_symbol, targetId),
    model,
    contrast,
    sep = "__"
  )
  stringr::str_replace_all(label, "[^A-Za-z0-9_.-]+", "_")
}

get_top_GWAS_gene_chromVAR_locus_check_tibble <- function(
  psbulk_GWAS_gene_chromVAR_results_tibble,
  GWAS_gene_chromVAR_credible_set_variants_tibble,
  GWAS_gene_chromVAR_L2G_tibble,
  open_targets_target_dataset_path,
  n_top = 10,
  FDR_threshold = 0.05,
  flank = 50000L
) {
  n_top <- as.integer(n_top %||% 0L)
  FDR_threshold <- as.numeric(FDR_threshold %||% 0.05)
  flank <- as.integer(flank %||% 50000L)
  if (n_top <= 0L || nrow(psbulk_GWAS_gene_chromVAR_results_tibble) == 0L) {
    return(empty_GWAS_gene_chromVAR_locus_check_tibble(psbulk_GWAS_gene_chromVAR_results_tibble))
  }

  selected_tibble <- psbulk_GWAS_gene_chromVAR_results_tibble |>
    dplyr::filter(!is.na(.data$FDR), .data$FDR <= .env$FDR_threshold) |>
    dplyr::arrange(
      .data$FDR,
      .data$PValue,
      dplyr::desc(abs(.data$logFC)),
      dplyr::desc(.data$sum_peak_gene_weight),
      dplyr::desc(.data$max_L2G_score)
    ) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::mutate(locus_plot_rank = dplyr::row_number())

  if (nrow(selected_tibble) == 0L) {
    return(empty_GWAS_gene_chromVAR_locus_check_tibble(psbulk_GWAS_gene_chromVAR_results_tibble))
  }

  selected_tibble |>
    dplyr::select(locus_plot_rank, GWAS_ID, targetId, model, contrast) |>
    purrr::pmap_dfr(\(locus_plot_rank, GWAS_ID, targetId, model, contrast) {
      get_GWAS_gene_chromVAR_locus_check_tibble(
        psbulk_GWAS_gene_chromVAR_results_tibble = psbulk_GWAS_gene_chromVAR_results_tibble,
        GWAS_gene_chromVAR_credible_set_variants_tibble = GWAS_gene_chromVAR_credible_set_variants_tibble,
        GWAS_gene_chromVAR_L2G_tibble = GWAS_gene_chromVAR_L2G_tibble,
        open_targets_target_dataset_path = open_targets_target_dataset_path,
        GWAS_ID = GWAS_ID,
        target_id = targetId,
        contrast = contrast,
        model = model,
        flank = flank
      ) |>
        dplyr::mutate(
          locus_plot_rank = locus_plot_rank,
          plot_name = make_GWAS_gene_chromVAR_locus_plot_name(
            rank = locus_plot_rank,
            GWAS_ID = .data$GWAS_ID,
            gene_symbol = .data$gene_symbol,
            targetId = .data$targetId,
            model = .data$model,
            contrast = .data$contrast
          )
        )
    })
}

plot_GWAS_gene_chromVAR_locus_tracks <- function(
  locus_check_tibble,
  GWAS_gene_chromVAR_credible_set_variants_tibble,
  GWAS_gene_chromVAR_L2G_tibble,
  consensus_peak_GRanges,
  fragments,
  metadata_tibble,
  group_cells_by_col = "PCA_harmony_SNN_cluster_cell_type"
) {
  if (nrow(locus_check_tibble) == 0L) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }
  if (!group_cells_by_col %in% colnames(metadata_tibble)) {
    stop("metadata_tibble is missing locus-track grouping column: ", group_cells_by_col)
  }

  fragments <- BPCells::select_cells(fragments, metadata_tibble$barcode_w_prefix)
  fragment_cell_names <- BPCells::cellNames(fragments)
  metadata <- metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(
      .data$barcode_w_prefix %in% fragment_cell_names,
      !is.na(.data[[group_cells_by_col]])
    ) |>
    dplyr::arrange(match(.data$barcode_w_prefix, fragment_cell_names))
  if (nrow(metadata) == 0L) {
    stop("No cells remain for GWAS-gene chromVAR locus track plotting.")
  }
  fragments <- BPCells::select_cells(fragments, metadata$barcode_w_prefix)

  cell_read_counts <- if ("atac_fragments" %in% colnames(metadata)) metadata$atac_fragments else metadata$nCount_ATAC
  groups <- metadata[[group_cells_by_col]]
  locus_rows <- split(locus_check_tibble, locus_check_tibble$plot_name)

  purrr::map(locus_rows, \(locus_row) {
    region <- GenomicRanges::GRanges(
      seqnames = locus_row$chr[[1]],
      ranges = IRanges::IRanges(
        start = locus_row$locus_start[[1]],
        end = locus_row$locus_end[[1]]
      )
    )
    L2G_tibble <- GWAS_gene_chromVAR_L2G_tibble |>
      dplyr::filter(
        .data$GWAS_ID == locus_row$GWAS_ID[[1]],
        .data$targetId == locus_row$targetId[[1]]
      )
    variant_tibble <- GWAS_gene_chromVAR_credible_set_variants_tibble |>
      dplyr::filter(.data$GWAS_ID == locus_row$GWAS_ID[[1]]) |>
      dplyr::semi_join(
        L2G_tibble |> dplyr::select(studyLocusId),
        by = "studyLocusId"
      )

    coverage_tibble <- BPCells::trackplot_coverage(
      fragments = fragments,
      region = region,
      groups = groups,
      cell_read_counts = cell_read_counts,
      group_order = gtools::mixedsort(unique(as.character(groups))),
      bins = 500,
      return_data = TRUE
    )

    BPCells::trackplot_combine(
      tracks = c(
        list(
          make_open_targets_target_locus_track(locus_row, region),
          make_BPCells_ATAC_coverage_track_from_tibble(coverage_tibble, region),
          make_consensus_peak_locus_track(consensus_peak_GRanges, region)
        ),
        make_open_targets_variant_PIP_tracks(variant_tibble, region)
      ),
      title = stringr::str_glue(
        "{locus_row$GWAS_ID[[1]]} {locus_row$gene_symbol[[1]]}: {locus_row$model[[1]]} / {locus_row$contrast[[1]]}; ",
        "logFC={round(locus_row$logFC[[1]], 3)}, ",
        "FDR={format(locus_row$FDR[[1]], scientific = TRUE, digits = 3)}, ",
        "L2G={round(locus_row$L2G_score[[1]], 3)}"
      )
    )
  })
}

make_open_targets_target_locus_track <- function(locus_check_tibble, region) {
  target_tibble <- locus_check_tibble |>
    dplyr::transmute(
      chr = .data$chr,
      start = .data$target_start,
      end = .data$target_end,
      gene_label = .data$gene_symbol
    )

  BPCells::trackplot_genome_annotation(
    loci = target_tibble,
    region = region,
    label_by = "gene_label",
    track_label = "Open Targets gene"
  )
}

#' Make open targets variant PIP track
#'
#' Draw an Open Targets variant PIP or weighted-PIP track for a locus.
#'
#' @param variant_tibble Variant-level tibble for one locus or track, including position, variant ID, and PIP/weight columns.
#' @param region Genomic region accepted by BPCells trackplot helpers; normalized internally to a single plotting interval.
#' @param track_label Text label shown on the rendered BPCells track.
#' @param weight_col Column in `variant_tibble` to draw as vertical PIP/weight values.
#' @param y_label Y-axis label for the track, such as `PIP` or `PIP_raw`.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

make_open_targets_variant_PIP_track <- function(variant_tibble, region, track_label = "Open Targets PIP", weight_col = "posteriorProbability", y_label = "PIP") {
  region <- BPCells:::normalize_ranges(region)
  if (nrow(variant_tibble) == 0L) {
    return(BPCells:::trackplot_empty(region, track_label))
  }
  if (!weight_col %in% names(variant_tibble)) {
    stop("variant_tibble is missing PIP weight column: ", weight_col)
  }

  variant_tibble <- variant_tibble |>
    dplyr::mutate(
      PIP_track_value = .data[[weight_col]],
      rank = dplyr::min_rank(dplyr::desc(.data$PIP_track_value)),
      variant_label = dplyr::if_else(.data$rank <= 3L, .data$variantId, "")
    )
  ymax <- max(variant_tibble$PIP_track_value, na.rm = TRUE) * 1.1

  BPCells:::wrap_trackplot(
    ggplot2::ggplot(variant_tibble, ggplot2::aes(x = .data$position, y = .data$PIP_track_value)) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data$position, y = 0, yend = .data$PIP_track_value),
        linewidth = 0.3
      ) +
      ggplot2::geom_point(size = 1.5) +
      ggrepel::geom_text_repel(
        data = dplyr::filter(variant_tibble, .data$variant_label != ""),
        ggplot2::aes(label = .data$variant_label),
        size = 2,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels = scales::label_number()) +
      ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
      ggplot2::labs(x = "Genomic Position (bp)", y = y_label) +
      BPCells:::trackplot_theme(),
    ggplot2::unit(1.4, "null"),
    region = region
  )
}

make_open_targets_variant_PIP_tracks <- function(variant_tibble, region) {
  tracks <- list(make_open_targets_variant_PIP_track(variant_tibble, region))

  if (
    "posteriorProbability_raw" %in% names(variant_tibble) &&
      any(abs(variant_tibble$posteriorProbability - variant_tibble$posteriorProbability_raw) > sqrt(.Machine$double.eps), na.rm = TRUE)
  ) {
    tracks$PIP_raw <- make_open_targets_variant_PIP_track(
      variant_tibble = variant_tibble,
      region = region,
      track_label = "Open Targets PIP_raw",
      weight_col = "posteriorProbability_raw",
      y_label = "PIP_raw"
    )
  }

  tracks
}

#' Make consensus peak locus track
#'
#' Draw consensus ATAC peaks overlapping a locus-check region.
#'
#' @param consensus_peak_GRanges Consensus peak GRanges whose names identify ATAC features in downstream matrices and track plots.
#' @param region Genomic region accepted by BPCells trackplot helpers; normalized internally to a single plotting interval.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

make_consensus_peak_locus_track <- function(consensus_peak_GRanges, region) {
  peak_GRanges <- IRanges::subsetByOverlaps(consensus_peak_GRanges, region)
  if (length(peak_GRanges) == 0L) {
    return(BPCells:::trackplot_empty(region, "Consensus peaks"))
  }

  peak_tibble <- GenomicRanges::as.data.frame(peak_GRanges) |>
    tibble::as_tibble() |>
    dplyr::transmute(
      chr = as.character(.data$seqnames),
      start = .data$start,
      end = .data$end
    )

  BPCells::trackplot_genome_annotation(
    loci = peak_tibble,
    region = region,
    track_label = "Consensus peaks"
  )
}

#' Get GWAS gene chromVAR peak weight matrix
#'
#' Build the peak-by-GWAS-gene feature weight matrix for chromVAR.
#'
#' @param GWAS_gene_chromVAR_peak_weight_records List of per-GWAS records with
#'   sparse peak-gene weight matrices.
#' @param GWAS_gene_chromVAR_feature_metadata_tibble Feature metadata defining
#'   which `feature_id`s pass support filters and should be retained.
#' @param RSE_ATAC RangedSummarizedExperiment for ATAC peaks, with row ranges aligned to peak-level matrices.
#' @return A named matrix-like object with rows and columns aligned to the input feature/cell identifiers.
#' @keywords internal

get_GWAS_gene_chromVAR_peak_weight_matrix <- function(GWAS_gene_chromVAR_peak_weight_records, GWAS_gene_chromVAR_feature_metadata_tibble, RSE_ATAC) {
  peak_names <- get_peak_names_from_GRanges(SummarizedExperiment::rowRanges(RSE_ATAC))
  retained_feature_ids <- GWAS_gene_chromVAR_feature_metadata_tibble |>
    dplyr::filter(support_filter_pass) |>
    dplyr::pull(feature_id)

  if (length(GWAS_gene_chromVAR_peak_weight_records) == 0 || length(retained_feature_ids) == 0) {
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      dims = c(length(peak_names), 0),
      dimnames = list(peak_names, character())
    ))
  }

  peak_weight_matrices <- GWAS_gene_chromVAR_peak_weight_records |>
    purrr::map(\(record) {
      mat <- record$peak_gene_weight_matrix
      colnames(mat) <- stringr::str_c(record$GWAS_ID, colnames(mat), sep = "__")
      mat
    })

  Matrix::Matrix(do.call(cbind, peak_weight_matrices), sparse = TRUE) |>
    align_peak_weights_to_RSE(RSE_ATAC = RSE_ATAC) |>
    (\(mat) mat[, retained_feature_ids, drop = FALSE])()
}

#' Format psbulk GWAS gene chromVAR results
#'
#' Join pseudobulk GWAS-gene chromVAR results to feature and target metadata.
#'
#' @param psbulk_GWAS_gene_chromVAR_results_tibble Raw pseudobulk model results
#'   keyed by chromVAR feature ID.
#' @param GWAS_gene_chromVAR_feature_metadata_tibble Feature metadata produced by
#'   `get_GWAS_gene_chromVAR_feature_metadata_tibble()`.
#' @return Formatted result tibble with original feature ID, GWAS ID, target
#'   metadata, support summaries, and model statistics. Empty input returns the
#'   same schema.
#' @keywords internal

format_psbulk_GWAS_gene_chromVAR_results <- function(psbulk_GWAS_gene_chromVAR_results_tibble, GWAS_gene_chromVAR_feature_metadata_tibble) {
  empty_result <- tibble::tibble(
    feature_id = character(),
    feature_id_original = character(),
    GWAS_ID = character(),
    variant_weighting_mode = character(),
    targetId = character(),
    gene_label = character(),
    gene_symbol = character(),
    gene_name = character(),
    gene_biotype = character(),
    model = character(),
    contrast = character(),
    logFC = numeric(),
    AveExpr = numeric(),
    t = numeric(),
    PValue = numeric(),
    FDR = numeric(),
    cell_type_subset = character(),
    n_samples = integer(),
    n_donors = integer(),
    n_paired_donors = integer(),
    min_group_n_samples = integer(),
    analysis_type = character(),
    n_loci = integer(),
    n_variants = integer(),
    n_peaks = integer(),
    sum_PIP_x_L2G = numeric(),
    sum_peak_gene_weight = numeric(),
    sum_peak_gene_weight_raw = numeric(),
    effect_weighting_lift = numeric(),
    max_L2G_score = numeric(),
    FDR_within_GWAS_model_contrast = numeric(),
    open_targets_release = character()
  )

  combined_results_tibble <- psbulk_GWAS_gene_chromVAR_results_tibble |>
    dplyr::bind_rows()

  if (nrow(combined_results_tibble) == 0) {
    return(empty_result)
  }

  combined_results_tibble |>
    dplyr::inner_join(
      GWAS_gene_chromVAR_feature_metadata_tibble |>
        dplyr::filter(support_filter_pass) |>
        dplyr::select(-support_filter_pass),
      by = "feature_id"
    ) |>
    dplyr::mutate(
      feature_id_original = feature_id,
      FDR_within_GWAS_model_contrast = stats::p.adjust(PValue, method = "BH"),
      FDR = FDR_within_GWAS_model_contrast,
      feature_id = gene_label,
      .by = c(GWAS_ID, model, contrast)
    ) |>
    dplyr::relocate(feature_id, GWAS_ID, variant_weighting_mode, targetId, gene_label, gene_symbol, gene_name, gene_biotype, model, contrast)
}

#' Plot psbulk GWAS gene chromVAR volcano
#'
#' Render locus-check plots for significant GWAS-gene chromVAR hits.
#'
#' @param psbulk_GWAS_gene_chromVAR_result_tibble One significant formatted
#'   GWAS-gene chromVAR result row used to choose GWAS, target, and contrast.
#' @param label_max_overlaps Maximum labels allowed by ggrepel before labels are dropped to reduce overplotting.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_psbulk_GWAS_gene_chromVAR_volcano <- function(psbulk_GWAS_gene_chromVAR_result_tibble, label_max_overlaps = 10) {
  if (dplyr::n_distinct(psbulk_GWAS_gene_chromVAR_result_tibble$GWAS_ID) != 1L || dplyr::n_distinct(psbulk_GWAS_gene_chromVAR_result_tibble$contrast) != 1L) {
    stop("plot_psbulk_GWAS_gene_chromVAR_volcano() expects one GWAS_ID and one contrast.")
  }

  uses_effect_weighting <- identical(unique(psbulk_GWAS_gene_chromVAR_result_tibble$variant_weighting_mode), "locus_lead_effect")
  color_label <- "log1p(raw peak-gene support)"

  plot_tibble <- psbulk_GWAS_gene_chromVAR_result_tibble |>
    dplyr::mutate(
      log10PValue = -log10(pmax(.data$PValue, .Machine$double.xmin)),
      support_color = log1p(.data$sum_peak_gene_weight_raw),
      significant = .data$FDR < 0.05,
      base_priority_score =
        dplyr::min_rank(.data$PValue) +
          dplyr::min_rank(dplyr::desc(abs(.data$logFC)))
    ) |>
    dplyr::mutate(
      priority_score = if (.env$uses_effect_weighting) {
        .data$base_priority_score + dplyr::min_rank(dplyr::desc(.data$effect_weighting_lift))
      } else {
        .data$base_priority_score
      },
      priority_rank = dplyr::min_rank(.data$priority_score),
      priority_label = dplyr::if_else(
        !is.na(.data$gene_label) & .data$gene_label != "" & !is.na(.data$priority_rank),
        stringr::str_c(.data$gene_label, " (", .data$priority_rank, ")"),
        ""
      )
    ) |>
    dplyr::ungroup()

  ggplot2::ggplot(plot_tibble, ggplot2::aes(x = .data$logFC, y = .data$log10PValue)) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = .data$support_color,
        size = .data$effect_weighting_lift,
        shape = .data$significant
      ),
      alpha = 0.65
    ) +
    ggplot2::scale_color_gradientn(
      colors = c("#D9D9D9", "#74A9CF", "#045A8D", "#A50F15"),
      name = color_label
    ) +
    ggplot2::scale_size_continuous(range = c(0.8, 3.2), name = "effect-weighting lift") +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "FDR < 0.05") +
    ggrepel::geom_text_repel(
      data = dplyr::filter(plot_tibble, .data$priority_label != ""),
      ggplot2::aes(label = .data$priority_label),
      size = 2,
      min.segment.length = 0,
      max.overlaps = label_max_overlaps,
      box.padding = 0.25,
      point.padding = 0.1
    ) +
    ggplot2::scale_x_continuous(limits = symmetric_limits) +
    ggplot2::labs(
      title = stringr::str_c(plot_tibble$model[[1]], ": ", plot_tibble$GWAS_ID[[1]], " / ", plot_tibble$contrast[[1]]),
      subtitle = stringr::str_c(plot_tibble$variant_weighting_mode[[1]], "; color = ", color_label, "; size = effect-weighting lift"),
      x = "logFC",
      y = "-log10(PValue)"
    ) +
    ggplot2::theme(legend.position = "bottom")
}

plot_psbulk_GWAS_gene_chromVAR_volcanoes <- function(psbulk_GWAS_gene_chromVAR_results_tibble, model_name) {
  model_results_tibble <- psbulk_GWAS_gene_chromVAR_results_tibble |>
    dplyr::filter(model == .env$model_name)

  if (nrow(model_results_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot_tibbles <- model_results_tibble |>
    dplyr::group_split(GWAS_ID, contrast)
  names(plot_tibbles) <- purrr::map_chr(
    plot_tibbles,
    \(plot_tibble) stringr::str_c(plot_tibble$GWAS_ID[[1]], plot_tibble$contrast[[1]], sep = "__")
  )

  purrr::map(plot_tibbles, \(plot_tibble) {
    plot_psbulk_GWAS_gene_chromVAR_volcano(plot_tibble)
  })
}

#' Plot psbulk GWAS gene chromVAR QC
#'
#' Add Open Targets gene metadata to GWAS-gene chromVAR features.
#'
#' @param psbulk_GWAS_gene_chromVAR_results_tibble Formatted pseudobulk results
#'   used to rank and select locus-check plots.
#' @param GWAS_gene_chromVAR_feature_metadata_tibble Feature metadata with GWAS,
#'   target, support, and Open Targets release fields.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_psbulk_GWAS_gene_chromVAR_QC <- function(psbulk_GWAS_gene_chromVAR_results_tibble, GWAS_gene_chromVAR_feature_metadata_tibble) {
  if (nrow(GWAS_gene_chromVAR_feature_metadata_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  feature_count_tibble <- dplyr::bind_rows(
    GWAS_gene_chromVAR_feature_metadata_tibble |>
      dplyr::count(GWAS_ID, name = "n_target_genes") |>
      dplyr::mutate(filter_stage = "before support filter"),
    GWAS_gene_chromVAR_feature_metadata_tibble |>
      dplyr::filter(support_filter_pass) |>
      dplyr::count(GWAS_ID, name = "n_target_genes") |>
      dplyr::mutate(filter_stage = "after support filter")
  )

  feature_count_plot <- feature_count_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = GWAS_ID, y = n_target_genes, fill = filter_stage)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
    ggplot2::labs(x = NULL, y = "Target genes", fill = NULL, title = "GWAS-gene annotations before and after support filtering") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  n_peaks_plot <- GWAS_gene_chromVAR_feature_metadata_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = n_peaks, fill = support_filter_pass)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::facet_wrap(~GWAS_ID, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ggplot2::labs(x = "Peaks per GWAS-gene feature", y = "Features", fill = "Passes filter", title = "GWAS-gene peak support")

  pip_plot <- GWAS_gene_chromVAR_feature_metadata_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = sum_PIP_x_L2G, fill = support_filter_pass)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::facet_wrap(~GWAS_ID, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ggplot2::labs(x = "sum PIP x L2G", y = "Features", fill = "Passes filter", title = "GWAS-gene genetic support")

  plots <- list(
    feature_counts = feature_count_plot,
    n_peaks = n_peaks_plot,
    sum_PIP_x_L2G = pip_plot
  )

  if (nrow(psbulk_GWAS_gene_chromVAR_results_tibble) > 0) {
    plots$PValue_density <- psbulk_GWAS_gene_chromVAR_results_tibble |>
      plot_psbulk_DX_PValue_density()

    top_genes_tibble <- psbulk_GWAS_gene_chromVAR_results_tibble |>
      dplyr::mutate(one_vs_rest_contrast = stringr::str_detect(contrast, "_vs_rest$")) |>
      dplyr::filter(!one_vs_rest_contrast | logFC > 0) |>
      dplyr::arrange(FDR_within_GWAS_model_contrast, PValue, dplyr::desc(abs(logFC))) |>
      dplyr::slice_head(n = 25, by = c(GWAS_ID, model, contrast)) |>
      dplyr::select(-one_vs_rest_contrast)

    if (nrow(top_genes_tibble) > 0) {
      plots$top_gene_support <- top_genes_tibble |>
        dplyr::mutate(label = stringr::str_c(gene_label, contrast, sep = " / ")) |>
        tidyr::pivot_longer(
          cols = c(n_loci, n_variants, n_peaks, sum_PIP_x_L2G, sum_peak_gene_weight, max_L2G_score),
          names_to = "support_metric",
          values_to = "value"
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = value, y = stats::reorder(label, value))) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::facet_wrap(~support_metric, scales = "free_x") +
        ggplot2::labs(x = NULL, y = NULL, title = "Support metrics for top GWAS-gene results")
    }
  }

  plots
}

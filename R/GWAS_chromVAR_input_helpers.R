#' Classify a configured GWAS source
#'
#' @param sourceId Open Targets study accession or local Parquet filename.
#' @return One of `open_targets` or `local_file`.
#' @keywords internal

classify_GWAS_source <- function(sourceId) {
  if (length(sourceId) != 1L || is.na(sourceId) || !nzchar(sourceId)) {
    stop("sourceId must be one non-empty value.")
  }
  if (stringr::str_starts(sourceId, "GCST")) {
    if (!stringr::str_detect(sourceId, "^GCST[0-9]+$")) {
      stop("Malformed Open Targets sourceId: ", sourceId)
    }
    return("open_targets")
  }
  "local_file"
}

#' Resolve a local GWAS source file
#'
#' @param sourceId Configured source identifier.
#' @param sourceType Source type returned by `classify_GWAS_source()`.
#' @return A normalized local filename, or `character()` for Open Targets.
#' @keywords internal

resolve_GWAS_source_file <- function(sourceId, sourceType) {
  if (identical(sourceType, "open_targets")) {
    return(character())
  }
  if (!identical(sourceType, "local_file")) {
    stop("Unsupported GWAS sourceType: ", sourceType)
  }
  normalizePath(sourceId, mustWork = TRUE)
}

#' Normalize GWAS variant weighting mode
#'
#' Normalize configured GWAS variant-weighting mode values.
#'
#' @param mode Configured variant weighting mode; supported values normalize to
#'   `raw_PIP` or `locus_lead_effect`.
#' @param GWAS_ID Configured GWAS label used in target names and plots.
#' @return One of `raw_PIP` or `locus_lead_effect`; unsupported non-empty values
#'   error with the GWAS ID for diagnosis.
#' @keywords internal

normalize_GWAS_variant_weighting_mode <- function(mode, GWAS_ID = NA_character_) {
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

#' Resolve an Open Targets fine-mapping method
#'
#' @param GWAS_config_tibble One configured Open Targets GWAS row.
#' @param open_targets_credible_set_dataset_path Open Targets credible-set
#'   dataset path.
#' @return The input row with resolved study, method, release, and source
#'   metadata.
#' @keywords internal

resolve_open_targets_GWAS_input_tibble <- function(
  GWAS_config_tibble,
  open_targets_credible_set_dataset_path
) {
  if (nrow(GWAS_config_tibble) != 1L) {
    stop("GWAS_config_tibble must contain exactly one Open Targets row.")
  }
  study_id <- GWAS_config_tibble$sourceId[[1]]
  available_methods <- arrow::open_dataset(open_targets_credible_set_dataset_path) |>
    dplyr::filter(studyType == "gwas", studyId == .env$study_id) |>
    dplyr::select(finemappingMethod) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::pull(finemappingMethod)
  if (length(available_methods) == 0L) {
    stop("Configured Open Targets sourceId is missing from credible_set: ", study_id)
  }

  requested_method <- GWAS_config_tibble$requested_finemappingMethod[[1]]
  finemapping_method <- if (!identical(requested_method, "auto")) {
    if (!requested_method %in% available_methods) {
      stop(
        "Requested Open Targets finemappingMethod is unavailable for GWAS_ID=",
        GWAS_config_tibble$GWAS_ID[[1]],
        ", studyId=",
        study_id,
        ". Requested: ",
        requested_method,
        ". Available: ",
        paste(available_methods, collapse = ", ")
      )
    }
    requested_method
  } else {
    method_priority <- c("SuSie", "SuSiE-inf", "PICS")
    selected <- method_priority[method_priority %in% available_methods]
    if (length(selected) == 0L) {
      stop(
        "No supported Open Targets finemappingMethod found for GWAS_ID=",
        GWAS_config_tibble$GWAS_ID[[1]],
        ", studyId=",
        study_id,
        ". Available: ",
        paste(available_methods, collapse = ", ")
      )
    }
    selected[[1]]
  }
  release <- basename(dirname(open_targets_credible_set_dataset_path))

  GWAS_config_tibble |>
    dplyr::transmute(
      GWAS_ID,
      sourceId,
      sourceType,
      studyId = .env$study_id,
      finemappingMethod = .env$finemapping_method,
      variant_weighting_mode,
      sourceRelease = .env$release,
      open_targets_release = .env$release,
      credibleSetProbability = 0.95
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
      sourceId = character(),
      sourceType = character(),
      sourceRelease = character(),
      open_targets_release = character(),
      studyId = character(),
      studyLocusId = character(),
      credibleSetIndex = integer(),
      finemappingMethod = character(),
      variant_weighting_mode = character(),
      confidence = character(),
      credibleSetProbability = numeric(),
      variantId = character(),
      chromosome = character(),
      position = integer(),
      variantRepresentation = character(),
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
        dplyr::select(
          GWAS_ID,
          sourceId,
          sourceType,
          sourceRelease,
          open_targets_release,
          studyId,
          finemappingMethod,
          credibleSetProbability,
          variant_weighting_mode
        ),
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
      position = suppressWarnings(as.integer(variant_parts[, 3])),
      variantRepresentation = "reference_alt"
    ) |>
    dplyr::filter(stringr::str_detect(chromosome, "^([0-9]+|X|Y|MT)$")) |>
    dplyr::select(
      GWAS_ID,
      sourceId,
      sourceType,
      sourceRelease,
      open_targets_release,
      studyId,
      studyLocusId,
      credibleSetIndex,
      finemappingMethod,
      variant_weighting_mode,
      confidence,
      credibleSetProbability,
      variantId,
      chromosome,
      position,
      variantRepresentation,
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

#' Validate a published local fine-mapped GWAS file
#'
#' @param variants_tibble One self-contained fine-mapped GWAS table.
#' @param expected_schema_version Supported schema version.
#' @return The input invisibly, or an error for an invalid contract.
#' @keywords internal

validate_local_finemapped_GWAS_tibble <- function(
  variants_tibble,
  expected_schema_version = 1L
) {
  required_columns <- c(
    "schemaVersion", "studyId", "studyLocusId", "credibleSetIndex",
    "studyType", "finemappingMethod", "confidence",
    "credibleSetProbability", "genomeBuild", "variantId", "chromosome",
    "position", "variantRepresentation", "posteriorProbability", "logBF",
    "pValueMantissa", "pValueExponent", "beta", "standardError",
    "r2Overall", "is95CredibleSet", "is99CredibleSet", "locusStart",
    "locusEnd", "sampleSize", "sourceUrl", "sourceSha256",
    "sourcePublication", "sourceRelease"
  )
  missing_columns <- setdiff(required_columns, names(variants_tibble))
  if (length(missing_columns) > 0L) {
    stop(
      "Local fine-mapped GWAS file is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (nrow(variants_tibble) == 0L) {
    stop("Local fine-mapped GWAS file must contain at least one variant.")
  }

  singleton_columns <- c(
    "schemaVersion", "studyId", "studyType", "finemappingMethod",
    "confidence", "credibleSetProbability", "genomeBuild", "sampleSize",
    "sourceUrl", "sourceSha256", "sourcePublication", "sourceRelease"
  )
  singleton_values <- purrr::map(variants_tibble[singleton_columns], unique)
  invalid_singletons <- names(singleton_values)[
    purrr::map_lgl(singleton_values, \(value) length(value) != 1L || is.na(value[[1]]))
  ]
  if (length(invalid_singletons) > 0L) {
    stop(
      "Local fine-mapped GWAS columns must each contain one non-missing value: ",
      paste(invalid_singletons, collapse = ", ")
    )
  }
  if (!identical(as.integer(singleton_values$schemaVersion[[1]]), as.integer(expected_schema_version))) {
    stop("Unsupported local fine-mapped GWAS schemaVersion.")
  }
  if (!identical(singleton_values$studyType[[1]], "gwas")) {
    stop("Local fine-mapped GWAS studyType must be gwas.")
  }
  if (!identical(singleton_values$genomeBuild[[1]], "GRCh38")) {
    stop("Local fine-mapped GWAS variants must use GRCh38.")
  }
  credible_set_probability <- singleton_values$credibleSetProbability[[1]]
  if (
    !is.finite(credible_set_probability) ||
      credible_set_probability <= 0 ||
      credible_set_probability > 1
  ) {
    stop("credibleSetProbability must be finite and in (0, 1].")
  }
  if (
    any(!is.finite(variants_tibble$posteriorProbability)) ||
      any(!dplyr::between(variants_tibble$posteriorProbability, 0, 1))
  ) {
    stop("posteriorProbability must be finite and in [0, 1].")
  }
  if (
    any(is.na(variants_tibble$position)) ||
      any(variants_tibble$position < variants_tibble$locusStart) ||
      any(variants_tibble$position > variants_tibble$locusEnd)
  ) {
    stop("Local GWAS variants must fall within their declared locus bounds.")
  }
  if (anyDuplicated(variants_tibble[c("studyId", "studyLocusId", "variantId")])) {
    stop("Local fine-mapped GWAS variant keys are duplicated.")
  }
  if (!stringr::str_detect(singleton_values$sourceSha256[[1]], "^[0-9a-f]{64}$")) {
    stop("sourceSha256 must be a lowercase 64-character SHA-256 checksum.")
  }
  if (
    isTRUE(all.equal(credible_set_probability, 0.95)) &&
      any(!variants_tibble$is95CredibleSet)
  ) {
    stop("Local 95% credible-set file contains rows outside the 95% set.")
  }
  if (
    isTRUE(all.equal(credible_set_probability, 0.99)) &&
      any(!variants_tibble$is99CredibleSet)
  ) {
    stop("Local 99% credible-set file contains rows outside the 99% set.")
  }
  invisible(variants_tibble)
}

#' Read a local fine-mapped GWAS source
#'
#' @param GWAS_config_tibble One local-file GWAS configuration row.
#' @param source_tibble Local Parquet contents.
#' @return Canonical source-neutral credible-set variants.
#' @keywords internal

get_local_finemapped_GWAS_variants_tibble <- function(
  GWAS_config_tibble,
  source_tibble
) {
  if (nrow(GWAS_config_tibble) != 1L) {
    stop("Expected one local GWAS configuration row.")
  }
  if (!identical(GWAS_config_tibble$variant_weighting_mode[[1]], "raw_PIP")) {
    stop(
      "Local fine-mapped GWAS files currently support only raw_PIP weighting: ",
      GWAS_config_tibble$GWAS_ID[[1]]
    )
  }

  validate_local_finemapped_GWAS_tibble(source_tibble)
  source_tibble |>
    dplyr::transmute(
      GWAS_ID = GWAS_config_tibble$GWAS_ID[[1]],
      sourceId = GWAS_config_tibble$sourceId[[1]],
      sourceType = GWAS_config_tibble$sourceType[[1]],
      sourceRelease,
      open_targets_release = NA_character_,
      studyId,
      studyLocusId,
      credibleSetIndex,
      finemappingMethod,
      variant_weighting_mode = GWAS_config_tibble$variant_weighting_mode[[1]],
      confidence,
      credibleSetProbability,
      variantId,
      chromosome = as.character(chromosome),
      position = as.integer(position),
      variantRepresentation,
      posteriorProbability,
      logBF,
      pValueMantissa,
      pValueExponent,
      beta,
      standardError,
      r2Overall,
      is95CredibleSet,
      is99CredibleSet,
      locusStart = as.integer(locusStart),
      locusEnd = as.integer(locusEnd)
    ) |>
    weight_open_targets_posterior_probability_by_locus_effect()
}

#' Summarize a local fine-mapped GWAS source
#'
#' @param variants_tibble Validated local source table.
#' @param GWAS_config_tibble One local-file GWAS configuration row.
#' @return One source-neutral metadata row.
#' @keywords internal

get_local_finemapped_GWAS_metadata_tibble <- function(
  variants_tibble,
  GWAS_config_tibble
) {
  variants_tibble |>
    dplyr::summarise(
      GWAS_ID = GWAS_config_tibble$GWAS_ID[[1]],
      sourceId = GWAS_config_tibble$sourceId[[1]],
      sourceType = GWAS_config_tibble$sourceType[[1]],
      studyId = dplyr::first(studyId),
      finemappingMethod = dplyr::first(finemappingMethod),
      variant_weighting_mode = GWAS_config_tibble$variant_weighting_mode[[1]],
      confidence = dplyr::first(confidence),
      credibleSetProbability = dplyr::first(credibleSetProbability),
      sourceRelease = dplyr::first(sourceRelease),
      sample_size = as.numeric(dplyr::first(sampleSize)),
      n_credible_set_loci = dplyr::n_distinct(studyLocusId),
      ancestry_EUR = NA_real_,
      ancestry_EAS = NA_real_,
      ancestry_AFR = NA_real_,
      ancestry_AMR = NA_real_,
      ancestry_SAS = NA_real_,
      ancestry_OTH = NA_real_
    )
}

#' Build one source-neutral GWAS chromVAR input record
#'
#' @param GWAS_config_tibble One source-classified GWAS configuration row.
#' @param local_source_file Tracked local filename, or `character()` for Open
#'   Targets.
#' @param open_targets_study_dataset_path Open Targets study dataset.
#' @param open_targets_credible_set_dataset_path Open Targets credible-set
#'   dataset.
#' @return One record containing metadata and credible-set `GRanges`.
#' @keywords internal

get_GWAS_chromVAR_input_record <- function(
  GWAS_config_tibble,
  local_source_file,
  open_targets_study_dataset_path,
  open_targets_credible_set_dataset_path
) {
  if (nrow(GWAS_config_tibble) != 1L) {
    stop("GWAS_config_tibble must contain exactly one GWAS row.")
  }

  if (identical(GWAS_config_tibble$sourceType[[1]], "open_targets")) {
    GWAS_input_tibble <- resolve_open_targets_GWAS_input_tibble(
      GWAS_config_tibble = GWAS_config_tibble,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    )
    variants_tibble <- get_open_targets_credible_set_variants_tibble(
      GWAS_inputs_tibble = GWAS_input_tibble,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    )
    metadata_tibble <- get_open_targets_GWAS_metadata_tibble(
      GWAS_inputs_tibble = GWAS_input_tibble |>
        dplyr::mutate(Category = NA_character_, .before = 1),
      open_targets_study_dataset_path = open_targets_study_dataset_path,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ) |>
      dplyr::transmute(
        GWAS_ID,
        sourceId = GWAS_config_tibble$sourceId[[1]],
        sourceType = GWAS_config_tibble$sourceType[[1]],
        studyId,
        finemappingMethod,
        variant_weighting_mode,
        confidence,
        credibleSetProbability = GWAS_input_tibble$credibleSetProbability[[1]],
        sourceRelease = open_targets_release,
        sample_size,
        n_credible_set_loci,
        dplyr::across(dplyr::matches("^ancestry_(EUR|EAS|AFR|AMR|SAS|OTH)$"))
      )
  } else if (identical(GWAS_config_tibble$sourceType[[1]], "local_file")) {
    source_tibble <- arrow::read_parquet(local_source_file)
    validate_local_finemapped_GWAS_tibble(source_tibble)
    variants_tibble <- get_local_finemapped_GWAS_variants_tibble(
      GWAS_config_tibble = GWAS_config_tibble,
      source_tibble = source_tibble
    )
    metadata_tibble <- get_local_finemapped_GWAS_metadata_tibble(
      variants_tibble = source_tibble,
      GWAS_config_tibble = GWAS_config_tibble
    )
    GWAS_input_tibble <- metadata_tibble
  } else {
    stop("Unsupported GWAS sourceType: ", GWAS_config_tibble$sourceType[[1]])
  }

  list(
    GWAS_ID = GWAS_config_tibble$GWAS_ID[[1]],
    sourceId = GWAS_config_tibble$sourceId[[1]],
    sourceType = GWAS_config_tibble$sourceType[[1]],
    studyId = GWAS_input_tibble$studyId[[1]],
    finemappingMethod = GWAS_input_tibble$finemappingMethod[[1]],
    variant_weighting_mode = GWAS_config_tibble$variant_weighting_mode[[1]],
    metadata_tibble = metadata_tibble,
    credible_set_GRanges = get_open_targets_credible_set_GRanges(variants_tibble)
  )
}

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
      sourceRelease = .env$release,
      open_targets_release = .env$release,
      credibleSetProbability = 0.95
    )
}

#' Get open targets credible set variants tibble
#'
#' Extract configured Open Targets credible-set variants for GWAS-gene chromVAR.
#'
#' @param GWAS_inputs_tibble GWAS configuration tibble with labels, Open Targets study IDs and finemapping methods.
#' @param open_targets_credible_set_dataset_path Local directory or dataset path
#'   for the Open Targets credible-set Parquet dataset.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_open_targets_credible_set_variants_tibble <- function(GWAS_inputs_tibble, open_targets_credible_set_dataset_path) {
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
          credibleSetProbability
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

  duplicated_variants <- variants_tibble |>
    dplyr::count(GWAS_ID, studyLocusId, variantId) |>
    dplyr::filter(n > 1)
  if (nrow(duplicated_variants) > 0) {
    stop("Open Targets method-only filtering duplicated GWAS_ID/studyLocusId/variantId rows. Review confidence values before processing.")
  }

  # Arrow returns dataset rows in a nondeterministic order.
  dplyr::arrange(variants_tibble, .data$studyLocusId, .data$variantId)
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
    )
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
      GWAS_inputs_tibble = GWAS_input_tibble,
      open_targets_study_dataset_path = open_targets_study_dataset_path,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ) |>
      dplyr::transmute(
        GWAS_ID,
        sourceId = GWAS_config_tibble$sourceId[[1]],
        sourceType = GWAS_config_tibble$sourceType[[1]],
        studyId,
        finemappingMethod,
        confidence,
        credibleSetProbability = GWAS_input_tibble$credibleSetProbability[[1]],
        sourceRelease = open_targets_release,
        sample_size,
        n_credible_set_loci,
        dplyr::across(dplyr::matches("^ancestry_(EUR|EAS|AFR|AMR|SAS|OTH)$"))
      )
  } else {
    source_tibble <- arrow::read_parquet(local_source_file)
    variants_tibble <- get_local_finemapped_GWAS_variants_tibble(
      GWAS_config_tibble = GWAS_config_tibble,
      source_tibble = source_tibble
    )
    metadata_tibble <- get_local_finemapped_GWAS_metadata_tibble(
      variants_tibble = source_tibble,
      GWAS_config_tibble = GWAS_config_tibble
    )
    GWAS_input_tibble <- metadata_tibble
  }

  list(
    GWAS_ID = GWAS_config_tibble$GWAS_ID[[1]],
    sourceId = GWAS_config_tibble$sourceId[[1]],
    sourceType = GWAS_config_tibble$sourceType[[1]],
    studyId = GWAS_input_tibble$studyId[[1]],
    finemappingMethod = GWAS_input_tibble$finemappingMethod[[1]],
    metadata_tibble = metadata_tibble,
    credible_set_GRanges = get_open_targets_credible_set_GRanges(variants_tibble)
  )
}

#' Get GWAS credible-set similarity matrix
#'
#' Sum each GWAS's posterior probabilities per variant, normalize them to unit
#' mass and compare GWAS pairs by their shared mass, `sum(pmin(PIP_1, PIP_2))`,
#' which is one minus the total variation distance. Only exact variant matches
#' count.
#'
#' @param GWAS_input_records List of GWAS input records.
#' @return Symmetric GWAS-by-GWAS matrix with values from 0 to 1.
#' @keywords internal

get_GWAS_credible_set_similarity_matrix <- function(GWAS_input_records) {
  GWAS_IDs <- purrr::map_chr(GWAS_input_records, "GWAS_ID")
  PIP_tibble <- GWAS_input_records |>
    purrr::map_dfr(\(record) tibble::tibble(
      GWAS_ID = record$GWAS_ID,
      variantId = record$credible_set_GRanges$variantId,
      PIP = record$credible_set_GRanges$posteriorProbability
    )) |>
    dplyr::summarise(PIP = sum(.data$PIP), .by = c(GWAS_ID, variantId)) |>
    dplyr::mutate(PIP = .data$PIP / sum(.data$PIP), .by = GWAS_ID)
  shared_tibble <- PIP_tibble |>
    dplyr::inner_join(PIP_tibble, by = "variantId", relationship = "many-to-many") |>
    dplyr::summarise(similarity = sum(pmin(.data$PIP.x, .data$PIP.y)), .by = c(GWAS_ID.x, GWAS_ID.y))

  similarity_matrix <- matrix(0, length(GWAS_IDs), length(GWAS_IDs), dimnames = list(GWAS_IDs, GWAS_IDs))
  similarity_matrix[cbind(shared_tibble$GWAS_ID.x, shared_tibble$GWAS_ID.y)] <- shared_tibble$similarity
  similarity_matrix
}

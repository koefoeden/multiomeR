#' Download open targets dataset
#'
#' Download all Parquet shards for one Open Targets dataset release.
#'
#' @param dataset_url Open Targets dataset URL ending in release and dataset path components.
#' @return Normalized local directory path under the targets store containing
#'   the downloaded Parquet files.
#' @keywords internal

download_open_targets_dataset <- function(dataset_url) {
  dataset_match <- stringr::str_match(dataset_url, "/platform/([^/]+)/output/([^/]+)/?$")
  if (any(is.na(dataset_match))) {
    stop("Open Targets dataset URL does not match the expected platform release/output/dataset shape: ", dataset_url)
  }

  open_targets_release <- dataset_match[, 2]
  dataset_name <- dataset_match[, 3]
  output_path <- file.path(targets::tar_config_get("store"), "files", "OpenTargets", open_targets_release, dataset_name)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  dataset_index <- paste(readLines(dataset_url, warn = FALSE), collapse = "\n")
  parquet_files <- stringr::str_match_all(dataset_index, 'href="([^"]+\\.parquet)"')[[1]][, 2]
  if (length(parquet_files) == 0) {
    stop("No Parquet files found at Open Targets dataset URL: ", dataset_url)
  }

  purrr::walk(parquet_files, \(file_name) {
    destfile <- file.path(output_path, file_name)
    if (file.exists(destfile) && file.info(destfile)$size > 0) {
      return(invisible(NULL))
    }

    tempfile <- paste0(destfile, ".tmp")
    if (file.exists(tempfile)) {
      unlink(tempfile)
    }

    utils::download.file(paste0(dataset_url, file_name), tempfile, mode = "wb")
    if (!file.exists(tempfile) || file.info(tempfile)$size == 0) {
      stop("Open Targets download failed or produced an empty file: ", file_name)
    }
    if (!file.rename(tempfile, destfile)) {
      stop("Could not move downloaded Open Targets file into place: ", destfile)
    }
  })

  normalizePath(output_path, mustWork = TRUE)
}

#' Get open targets GWAS metadata tibble
#'
#' Join configured GWAS inputs with Open Targets study and credible-set summaries.
#'
#' @param GWAS_inputs_tibble GWAS configuration tibble with labels, Open Targets study IDs, finemapping methods, and weighting modes.
#' @param open_targets_study_dataset_path Local directory or dataset path for the
#'   Open Targets `study` Parquet dataset.
#' @param open_targets_credible_set_dataset_path Local directory or dataset path
#'   for the Open Targets credible-set Parquet dataset.
#' @return A GWAS metadata tibble with sample-size and credible-set counts for
#'   configured study/finemapping-method combinations. Missing studies error.
#' @keywords internal

get_open_targets_GWAS_metadata_tibble <- function(GWAS_inputs_tibble, open_targets_study_dataset_path, open_targets_credible_set_dataset_path) {
  input_studies <- unique(GWAS_inputs_tibble$studyId)

  study_tibble <- arrow::open_dataset(open_targets_study_dataset_path) |>
    dplyr::filter(studyType == "gwas", studyId %in% input_studies) |>
    dplyr::select(
      studyId,
      nCases,
      nControls,
      nSamples,
      discoverySamples
    ) |>
    dplyr::collect()

  missing_studies <- setdiff(input_studies, study_tibble$studyId)
  if (length(missing_studies) > 0) {
    stop("Configured Open Targets studyId(s) missing from study: ", paste(missing_studies, collapse = ", "))
  }

  credible_set_summary_tibble <- arrow::open_dataset(open_targets_credible_set_dataset_path) |>
    dplyr::filter(studyType == "gwas", studyId %in% input_studies) |>
    dplyr::select(studyId, finemappingMethod, confidence, studyLocusId, sampleSize) |>
    dplyr::collect() |>
    dplyr::semi_join(
      GWAS_inputs_tibble |> dplyr::select(studyId, finemappingMethod),
      by = c("studyId", "finemappingMethod")
    ) |>
    dplyr::group_by(studyId, finemappingMethod) |>
    dplyr::summarise(
      confidence = paste(sort(unique(confidence)), collapse = "; "),
      n_credible_set_loci = dplyr::n_distinct(studyLocusId),
      max_locus_sample_size = suppressWarnings(max(sampleSize, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(max_locus_sample_size = dplyr::if_else(is.infinite(max_locus_sample_size), NA_real_, max_locus_sample_size))

  ancestry_alias_tibble <- tibble::tribble(
    ~ancestry, ~ancestry_group,
    "European", "EUR",
    "Finnish", "EUR",
    "East Asian", "EAS",
    "African American or Afro-Caribbean", "AFR",
    "African unspecified", "AFR",
    "Sub-Saharan African", "AFR",
    "Hispanic or Latin American", "AMR",
    "South Asian", "SAS",
    "Central Asian", "OTH",
    "Greater Middle Eastern", "OTH",
    "Other", "OTH",
    "NR", "OTH"
  )

  ancestry_wide_tibble <- study_tibble |>
    dplyr::select(studyId, discoverySamples) |>
    tidyr::unnest(discoverySamples, keep_empty = TRUE) |>
    dplyr::mutate(sampleSize = as.numeric(sampleSize), ancestry = dplyr::coalesce(ancestry, "NR")) |>
    dplyr::left_join(ancestry_alias_tibble, by = "ancestry") |>
    dplyr::mutate(ancestry_group = dplyr::coalesce(ancestry_group, "OTH")) |>
    dplyr::summarise(ancestry_sample_size = sum(sampleSize, na.rm = TRUE), .by = c(studyId, ancestry_group)) |>
    tidyr::complete(studyId, ancestry_group = c("EUR", "EAS", "AFR", "AMR", "SAS", "OTH"), fill = list(ancestry_sample_size = 0)) |>
    dplyr::mutate(
      total_ancestry_sample_size = sum(ancestry_sample_size),
      ancestry_fraction = dplyr::if_else(total_ancestry_sample_size > 0, ancestry_sample_size / total_ancestry_sample_size, 0),
      .by = studyId
    ) |>
    dplyr::select(studyId, ancestry_group, ancestry_fraction) |>
    tidyr::pivot_wider(names_from = ancestry_group, values_from = ancestry_fraction, names_prefix = "ancestry_", values_fill = 0)

  GWAS_inputs_tibble |>
    dplyr::left_join(study_tibble, by = "studyId") |>
    dplyr::left_join(credible_set_summary_tibble, by = c("studyId", "finemappingMethod")) |>
    dplyr::left_join(ancestry_wide_tibble, by = "studyId") |>
    dplyr::mutate(
      sample_size = dplyr::coalesce(as.numeric(nSamples), as.numeric(nCases) + as.numeric(nControls), max_locus_sample_size)
    ) |>
    dplyr::select(
      Category,
      GWAS_ID,
      studyId,
      finemappingMethod,
      variant_weighting_mode,
      confidence,
      open_targets_release,
      sample_size,
      n_credible_set_loci,
      dplyr::matches("^ancestry_(EUR|EAS|AFR|AMR|SAS|OTH)$")
    )
}

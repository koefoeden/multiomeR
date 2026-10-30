#' Download open targets dataset
#'
#' Download all Parquet shards for one Open Targets dataset release.
#'
#' @param dataset_url Open Targets dataset URL ending in release and dataset path components.
#' @return Directory path under the targets store containing the downloaded
#'   Parquet files.
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
  parquet_files <- parse_open_targets_listing(dataset_index)

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

  output_path
}

#' Share GWAS samples among ancestry groups
#'
#' The GWAS Catalog gives a sample group's size, not its split among the
#' ancestry categories it lists. A group whose categories fall in several
#' ancestry groups is therefore counted as mixed (`MIX`). "NR" (not reported)
#' is ignored, and groups that report no ancestry are left out.
#'
#' @param sample_tibble One row per sample group and ancestry category, with
#'   `sample_group`, `ancestry` and the group's `sample_size`.
#' @return One row of `ancestry_*` shares, all missing when no group reports an
#'   ancestry.
#' @keywords internal

get_GWAS_ancestry_share_tibble <- function(sample_tibble) {
  ancestry_groups <- c("EUR", "EAS", "AFR", "AMR", "SAS", "OTH", "MIX")
  group_tibble <- sample_tibble |>
    dplyr::filter(ancestry != "NR") |>
    dplyr::mutate(ancestry_group = dplyr::recode_values(
      ancestry,
      c("European", "Finnish") ~ "EUR",
      "East Asian" ~ "EAS",
      c("African American or Afro-Caribbean", "African unspecified", "Sub-Saharan African") ~ "AFR",
      "Hispanic or Latin American" ~ "AMR",
      "South Asian" ~ "SAS",
      default = "OTH"
    )) |>
    dplyr::summarise(
      ancestry_group = if (dplyr::n_distinct(ancestry_group) == 1L) dplyr::first(ancestry_group) else "MIX",
      sample_size = dplyr::first(sample_size),
      .by = sample_group
    )
  sizes <- tapply(group_tibble$sample_size, factor(group_tibble$ancestry_group, ancestry_groups), sum, default = 0)
  shares <- if (nrow(group_tibble) > 0L) sizes / sum(sizes) else sizes * NA_real_
  tibble::as_tibble(as.list(rlang::set_names(as.numeric(shares), paste0("ancestry_", ancestry_groups))))
}

#' Get open targets GWAS metadata tibble
#'
#' Join configured GWAS inputs with Open Targets study and credible-set summaries.
#' Open Targets splits a GWAS Catalog sample group equally among its listed
#' ancestries, so rows of equal size within a study are taken as one group.
#'
#' @param GWAS_inputs_tibble GWAS configuration tibble with labels, Open Targets study IDs, and finemapping methods.
#' @param open_targets_study_dataset_path Local directory or dataset path for the
#'   Open Targets `study` Parquet dataset.
#' @param open_targets_credible_set_dataset_path Local directory or dataset path
#'   for the Open Targets credible-set Parquet dataset.
#' @return A GWAS metadata tibble with sample sizes, ancestry shares and credible-set confidence for
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
      max_locus_sample_size = suppressWarnings(max(sampleSize, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(max_locus_sample_size = dplyr::if_else(is.infinite(max_locus_sample_size), NA_real_, max_locus_sample_size))

  ancestry_wide_tibble <- study_tibble |>
    dplyr::select(studyId, discoverySamples) |>
    tidyr::unnest(discoverySamples, keep_empty = TRUE) |>
    dplyr::mutate(sampleSize = as.numeric(sampleSize), ancestry = dplyr::coalesce(ancestry, "NR")) |>
    dplyr::reframe(
      get_GWAS_ancestry_share_tibble(tibble::tibble(
        sample_group = sampleSize,
        ancestry = ancestry,
        sample_size = stats::ave(sampleSize, sampleSize, FUN = sum)
      )),
      .by = studyId
    )

  GWAS_inputs_tibble |>
    dplyr::left_join(study_tibble, by = "studyId") |>
    dplyr::left_join(credible_set_summary_tibble, by = c("studyId", "finemappingMethod")) |>
    dplyr::left_join(ancestry_wide_tibble, by = "studyId") |>
    dplyr::mutate(
      sample_size = dplyr::coalesce(as.numeric(nSamples), as.numeric(nCases) + as.numeric(nControls), max_locus_sample_size)
    ) |>
    dplyr::select(
      GWAS_ID,
      studyId,
      finemappingMethod,
      confidence,
      open_targets_release,
      sample_size,
      dplyr::starts_with("ancestry_")
    )
}

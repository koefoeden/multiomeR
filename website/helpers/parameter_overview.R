library(dplyr)
library(glue)
library(jsonlite)
library(readr)
library(stringr)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

find_parameter_manifest_root <- function(start_dir = ".") {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current_dir, "cfg_pipeline_parameters.tsv"))) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Could not find cfg_pipeline_parameters.tsv in current or parent directories.", call. = FALSE)
    }
    current_dir <- parent_dir
  }
}

read_parameter_manifest <- function(manifest_file) {
  manifest <- readr::read_tsv(
    manifest_file,
    col_types = readr::cols(
      .default = readr::col_character(),
      allow_missing_after_inheritance = readr::col_logical()
    ),
    na = character(),
    show_col_types = FALSE
  )

  names(manifest) <- stringr::str_trim(names(manifest))
  dplyr::mutate(manifest, dplyr::across(where(is.character), stringr::str_trim))
}

parameter_topic_order <- function(scope) {
  switch(
    scope,
    aggregation = c(
      "required",
      "GEX processing",
      "ATAC processing",
      "Multimodal processing",
      "subgroup analysis",
      "plotting/UMAP",
      "miscellaneous"
    ),
    NULL
  )
}

parameter_overview_data <- function(manifest, scope, topic_order = NULL) {
  parameter_tibble <- manifest |>
    dplyr::mutate(manifest_row = dplyr::row_number()) |>
    dplyr::filter(.data$scope == .env$scope)

  if (nrow(parameter_tibble) == 0) {
    stop("No parameters found for manifest scope: ", scope, call. = FALSE)
  }

  topic_order <- topic_order %||% parameter_topic_order(scope) %||% unique(parameter_tibble$topic)

  parameter_tibble <- parameter_tibble |>
    dplyr::mutate(
      default_is_null = .data$default_value == "NULL",
      must_specify = !.data$allow_missing_after_inheritance & .data$default_is_null,
      resolved_value_required = !.data$allow_missing_after_inheritance,
      status = dplyr::case_when(
        .data$must_specify ~ "Must specify",
        .data$resolved_value_required ~ "Defaulted",
        TRUE ~ "Optional"
      ),
      sort_status = dplyr::case_when(
        .data$must_specify ~ 1L,
        .data$resolved_value_required ~ 2L,
        TRUE ~ 3L
      ),
      topic_sort = match(.data$topic, topic_order),
      topic_sort = dplyr::coalesce(.data$topic_sort, length(topic_order) + 1L)
    ) |>
    dplyr::arrange(.data$topic_sort, .data$sort_status, .data$topic, .data$param_name) |>
    dplyr::select(
      param_name,
      data_type,
      cardinality,
      default_value,
      allowed_values,
      examples,
      topic,
      part_of,
      description,
      status
    )

  list(parameters = parameter_tibble)
}

render_parameter_overview_fragment <- function(
  overview_data,
  scope,
  search_placeholder = NULL
) {
  search_placeholder <- search_placeholder %||% "Search names, descriptions or values"
  data_json <- jsonlite::toJSON(overview_data, auto_unbox = TRUE, null = "null")
  # Prevent a manifest value from closing the inert JSON script element.
  data_json <- gsub("<", "\\u003c", data_json, fixed = TRUE)
  script <- paste(readLines(file.path(find_parameter_manifest_root(), "website", "helpers", "parameter_overview.js")), collapse = "\n")
  search_placeholder <- htmltools::htmlEscape(search_placeholder, attribute = TRUE)

  glue::glue(
'<div class="parameter-overview" data-parameter-overview>
<div class="parameter-overview-toolbar">
<label class="parameter-search-control"><span>Find a parameter</span>
<input class="parameter-overview-search" type="search" placeholder="{search_placeholder}"></label>
<label class="parameter-topic-control"><span>Topic</span>
<select><option value="">All topics</option></select></label>
</div>
<div class="parameter-filter-bar">
<div class="parameter-filter-group" role="group" aria-label="Filter by requirement">
<button type="button" class="parameter-filter is-active" data-status="all" aria-pressed="true">All <span></span></button>
<button type="button" class="parameter-filter" data-status="Must specify" aria-pressed="false">Required <span></span></button>
<button type="button" class="parameter-filter" data-status="Defaulted" aria-pressed="false">Defaulted <span></span></button>
<button type="button" class="parameter-filter" data-status="Optional" aria-pressed="false">Optional <span></span></button>
</div>
<button type="button" class="parameter-reset">Reset filters</button>
</div>
<p class="parameter-overview-help">Required: supply a value directly or through inheritance. Defaulted: a value is provided. Optional: may remain NULL. Open a parameter for its type, allowed values and example.</p>
<p class="parameter-overview-summary" aria-live="polite" aria-atomic="true"></p>
<div class="parameter-topic-list"></div>
<p class="parameter-empty-state" hidden>No parameters match. Try another search or reset the filters.</p>
<noscript><p>Enable JavaScript to browse this reference, or read <a href="https://github.com/koefoeden/multiomeR/blob/main/cfg_pipeline_parameters.tsv">the parameter manifest</a>.</p></noscript>
<script type="application/json" class="parameter-overview-data">{data_json}</script>
</div>
<script>{script}</script>'
  )
}

emit_parameter_overview <- function(
  scope,
  manifest_file = file.path(find_parameter_manifest_root(), "cfg_pipeline_parameters.tsv"),
  search_placeholder = NULL
) {
  manifest <- read_parameter_manifest(manifest_file)
  overview_data <- parameter_overview_data(manifest, scope)
  cat("```{=html}\n", render_parameter_overview_fragment(overview_data, scope, search_placeholder), "\n```\n", sep = "")
}

render_parameter_overview_document <- function(overview_data, scope, search_placeholder = NULL) {
  scope_label <- stringr::str_to_sentence(stringr::str_replace_all(scope, "_", " "))
  fragment <- render_parameter_overview_fragment(overview_data, scope, search_placeholder)

  glue::glue(
'<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{scope_label} YAML parameters | multiomeR proof of concept</title>
  <link rel="stylesheet" href="styles.css">
</head>
<body class="parameter-overview-standalone">
<main>
{fragment}
</main>
</body>
</html>'
  )
}

render_parameter_overview_file <- function(
  scope = "aggregation",
  manifest_file = file.path(find_parameter_manifest_root(), "cfg_pipeline_parameters.tsv"),
  output_file = NULL,
  search_placeholder = NULL
) {
  repo_root <- find_parameter_manifest_root()
  if (is.null(output_file)) {
    output_name <- if (identical(scope, "aggregation")) {
      "aggregation_parameter_overview_poc.html"
    } else {
      stringr::str_glue("{scope}_parameter_overview_poc.html")
    }
    output_file <- file.path(repo_root, "website", output_name)
  }

  manifest <- read_parameter_manifest(manifest_file)
  overview_data <- parameter_overview_data(manifest, scope)
  html <- render_parameter_overview_document(overview_data, scope, search_placeholder)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_file, useBytes = TRUE)
  message("Wrote ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
  invisible(output_file)
}

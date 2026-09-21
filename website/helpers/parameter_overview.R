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

parameter_overview_data <- function(manifest, scope = NULL) {
  parameter_tibble <- manifest |>
    dplyr::filter(is.null(.env$scope) | .data$scope %in% .env$scope) |>
    dplyr::mutate(
      status = dplyr::case_when(
        .data$default_value != "NULL" ~ "Defaulted",
        !.data$allow_missing_after_inheritance ~ "Required",
        TRUE ~ "Optional"
      )
    ) |>
    dplyr::select(
      scope, param_name, short_name, data_type, cardinality, default_value,
      allowed_values, examples, part_of, description, status
    )

  if (nrow(parameter_tibble) == 0) {
    stop("No parameters found for manifest scope: ", scope, call. = FALSE)
  }
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
<div class="parameter-workflows" role="tablist" aria-label="Workflow"></div>
<div class="parameter-overview-toolbar">
<label class="parameter-search-control"><span>Search this workflow</span>
<input class="parameter-overview-search" type="search" placeholder="{search_placeholder}"></label>
<button type="button" class="parameter-reset">Clear search</button>
</div>
<p class="parameter-overview-help">Required: supply a value directly or through inheritance. Defaulted: a non-null default is provided. Optional: may remain unset. Within each group, parameters are sorted by short name.</p>
<p class="parameter-overview-summary" aria-live="polite" aria-atomic="true"></p>
<div class="parameter-panel" id="parameter-panel" role="tabpanel">
<div class="parameter-section-list"></div>
<p class="parameter-empty-state" hidden>No parameters match in this workflow. Clear the search or choose another tab.</p>
</div>
<noscript><p>Enable JavaScript to browse this reference, or read <a href="https://github.com/koefoeden/multiomeR/blob/main/cfg_pipeline_parameters.tsv">the parameter manifest</a>.</p></noscript>
<script type="application/json" class="parameter-overview-data">{data_json}</script>
</div>
<script>{script}</script>'
  )
}

emit_parameter_overview <- function(
  scope,
  manifest_file = file.path(find_parameter_manifest_root(), "website", "data", "public_defaults", "cfg_pipeline_parameters.tsv"),
  search_placeholder = NULL
) {
  manifest <- read_parameter_manifest(manifest_file)
  overview_data <- parameter_overview_data(manifest, scope)
  cat("```{=html}\n", render_parameter_overview_fragment(overview_data, scope, search_placeholder), "\n```\n", sep = "")
}

render_parameter_overview_document <- function(overview_data, scope, search_placeholder = NULL) {
  scope_label <- if (is.null(scope)) "All" else stringr::str_to_sentence(stringr::str_replace_all(scope, "_", " "))
  styles <- paste(readLines(file.path(find_parameter_manifest_root(), "website", "styles.css")), collapse = "\n")
  fragment <- render_parameter_overview_fragment(overview_data, scope, search_placeholder)

  glue::glue(
'<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{scope_label} parameters | multiomeR</title>
  <style>{styles}</style>
</head>
<body class="parameter-overview-standalone">
<main>
<h1>{scope_label} parameters</h1>
<p>Search settings across the pipeline and optional modules. Open a parameter for details or use its permalink to share it.</p>
{fragment}
</main>
</body>
</html>'
  )
}

render_parameter_overview_file <- function(
  scope = NULL,
  manifest_file = file.path(find_parameter_manifest_root(), "website", "data", "public_defaults", "cfg_pipeline_parameters.tsv"),
  output_file = NULL,
  search_placeholder = NULL
) {
  repo_root <- find_parameter_manifest_root()
  if (is.null(output_file)) {
    output_file <- file.path(repo_root, "website", "parameters.html")
  }

  manifest <- read_parameter_manifest(manifest_file)
  overview_data <- parameter_overview_data(manifest, scope)
  html <- render_parameter_overview_document(overview_data, scope, search_placeholder)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_file, useBytes = TRUE)
  message("Wrote ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
  invisible(output_file)
}

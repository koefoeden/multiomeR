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
  scope_label <- stringr::str_to_sentence(stringr::str_replace_all(scope, "_", " "))
  search_placeholder <- search_placeholder %||% stringr::str_glue("Search {scope_label} parameters")
  data_json <- jsonlite::toJSON(overview_data, auto_unbox = TRUE, pretty = TRUE, null = "null")

  glue::glue(
'<div class="parameter-overview" data-parameter-overview>
<div class="parameter-overview-toolbar">
<label class="parameter-search-control">
<span class="visually-hidden">{search_placeholder}</span>
<input class="parameter-overview-search" type="search" placeholder="{search_placeholder}">
</label>
<div class="parameter-filter-group" role="group" aria-label="Filter parameters by requirement status">
<button type="button" class="parameter-filter is-active" data-status="all" aria-pressed="true">
All <span class="parameter-filter-count"></span>
</button>
<button type="button" class="parameter-filter" data-status="Must specify" aria-pressed="false">
Required <span class="parameter-filter-count"></span>
</button>
<button type="button" class="parameter-filter" data-status="Defaulted" aria-pressed="false">
Defaulted <span class="parameter-filter-count"></span>
</button>
<button type="button" class="parameter-filter" data-status="Optional" aria-pressed="false">
Optional <span class="parameter-filter-count"></span>
</button>
</div>
</div>
<p class="parameter-overview-summary" aria-live="polite"></p>
<div class="parameter-topic-list"></div>
<div class="parameter-empty-state" role="status">No {scope_label} parameters match the current filters.</div>
<noscript><p class="parameter-noscript">Enable JavaScript to browse the parameter reference.</p></noscript>
<script type="application/json" class="parameter-overview-data">
{data_json}
</script>
</div>
<script>
(() => {{
  const widget = document.currentScript.previousElementSibling;
  const data = JSON.parse(widget.querySelector(".parameter-overview-data").textContent);
  const parameters = data.parameters;
  const searchInput = widget.querySelector(".parameter-overview-search");
  const filterButtons = Array.from(widget.querySelectorAll(".parameter-filter"));
  const summary = widget.querySelector(".parameter-overview-summary");
  const topicList = widget.querySelector(".parameter-topic-list");
  const emptyState = widget.querySelector(".parameter-empty-state");
  let activeStatus = "all";

  const text = (value) => value === null || value === undefined ? "" : String(value);
  const hasValue = (value) => text(value).trim() !== "";
  const escapeHtml = (value) => text(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll(\'"\', "&quot;")
    .replaceAll("\'", "&#039;");

  function parameterSearchText(parameter) {{
    return [
      parameter.param_name,
      parameter.description,
      parameter.topic,
      parameter.part_of,
      parameter.data_type,
      parameter.cardinality,
      parameter.default_value,
      parameter.allowed_values,
      parameter.examples,
      parameter.status
    ].map(text).join(" ").toLowerCase();
  }}

  function badgeClass(status) {{
    if (status === "Must specify") return "must";
    if (status === "Defaulted") return "defaulted";
    return "optional";
  }}

  function statusLabel(status) {{
    return status === "Must specify" ? "Required" : status;
  }}

  function renderValue(value, fallback) {{
    return hasValue(value)
      ? `<span class="parameter-detail-value">${{escapeHtml(value)}}</span>`
      : `<span class="parameter-detail-value parameter-empty-value">${{escapeHtml(fallback)}}</span>`;
  }}

  function renderParameterRow(parameter, searchIsActive) {{
    const statusClass = badgeClass(parameter.status);
    const open = searchIsActive ? " open" : "";

    return `
      <details class="parameter-row ${{statusClass}}" data-param="${{escapeHtml(parameter.param_name)}}"${{open}}>
        <summary class="parameter-head">
          <span class="parameter-name-wrap">
            <code class="parameter-name">${{escapeHtml(parameter.param_name)}}</code>
            <span class="parameter-status ${{statusClass}}">${{escapeHtml(statusLabel(parameter.status))}}</span>
          </span>
          <span class="parameter-description">${{escapeHtml(parameter.description)}}</span>
        </summary>
        <div class="parameter-details">
          <div class="parameter-detail">
            <span class="parameter-detail-label">Default</span>
            ${{renderValue(parameter.default_value, "missing allowed")}}
          </div>
          <div class="parameter-detail">
            <span class="parameter-detail-label">Type</span>
            <span class="parameter-detail-value">${{escapeHtml(parameter.data_type)}} / ${{escapeHtml(parameter.cardinality)}}</span>
          </div>
          <div class="parameter-detail">
            <span class="parameter-detail-label">Allowed values</span>
            ${{renderValue(parameter.allowed_values, "any value matching the type")}}
          </div>
          <div class="parameter-detail">
            <span class="parameter-detail-label">Example</span>
            ${{renderValue(parameter.examples, "no example yet")}}
          </div>
          <div class="parameter-detail">
            <span class="parameter-detail-label">Used by</span>
            ${{renderValue(parameter.part_of, "not assigned")}}
          </div>
        </div>
      </details>
    `;
  }}

  function renderTopicGroup(topic, topicParameters, searchIsActive) {{
    const requiredCount = topicParameters.filter((parameter) => parameter.status === "Must specify").length;
    const shouldOpen = searchIsActive || activeStatus !== "all" || requiredCount > 0;
    const open = shouldOpen ? " open" : "";
    const requiredText = requiredCount > 0 ? ` · ${{requiredCount}} required` : "";

    return `
      <details class="parameter-topic"${{open}}>
        <summary>
          <span class="parameter-topic-title">${{escapeHtml(topic)}}</span>
          <span class="parameter-topic-count">${{topicParameters.length}}${{requiredText}}</span>
        </summary>
        <div class="parameter-list">
          ${{topicParameters.map((parameter) => renderParameterRow(parameter, searchIsActive)).join("")}}
        </div>
      </details>
    `;
  }}

  function updateFilterCounts() {{
    filterButtons.forEach((button) => {{
      const status = button.dataset.status;
      const count = status === "all"
        ? parameters.length
        : parameters.filter((parameter) => parameter.status === status).length;
      button.querySelector(".parameter-filter-count").textContent = count;
    }});
  }}

  function render() {{
    const query = searchInput.value.trim().toLowerCase();
    const visible = parameters.filter((parameter) => {{
      const matchesSearch = !query || parameterSearchText(parameter).includes(query);
      const matchesStatus = activeStatus === "all" || parameter.status === activeStatus;
      return matchesSearch && matchesStatus;
    }});
    const topicGroups = new Map();

    visible.forEach((parameter) => {{
      if (!topicGroups.has(parameter.topic)) topicGroups.set(parameter.topic, []);
      topicGroups.get(parameter.topic).push(parameter);
    }});

    topicList.innerHTML = Array.from(topicGroups)
      .map(([topic, topicParameters]) => renderTopicGroup(topic, topicParameters, query.length > 0))
      .join("");
    emptyState.style.display = visible.length === 0 ? "block" : "none";
    topicList.style.display = visible.length === 0 ? "none" : "grid";

    const required = parameters.filter((parameter) => parameter.status === "Must specify").length;
    summary.textContent = query || activeStatus !== "all"
      ? `Showing ${{visible.length}} of ${{parameters.length}} parameters`
      : `${{parameters.length}} parameters · ${{required}} required`;
  }}

  searchInput.addEventListener("input", render);
  filterButtons.forEach((button) => {{
    button.addEventListener("click", () => {{
      activeStatus = button.dataset.status;
      filterButtons.forEach((item) => {{
        const active = item === button;
        item.classList.toggle("is-active", active);
        item.setAttribute("aria-pressed", String(active));
      }});
      render();
    }});
  }});
  updateFilterCounts();
  render();
}})();
</script>'
  )
}

emit_parameter_overview <- function(
  scope,
  manifest_file = file.path(find_parameter_manifest_root(), "cfg_pipeline_parameters.tsv"),
  search_placeholder = NULL
) {
  manifest <- read_parameter_manifest(manifest_file)
  overview_data <- parameter_overview_data(manifest, scope)
  cat(render_parameter_overview_fragment(overview_data, scope, search_placeholder))
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

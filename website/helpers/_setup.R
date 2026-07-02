knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%",
  eval = FALSE
)

repo_root <- if (file.exists("DESCRIPTION")) {
  normalizePath(".", winslash = "/", mustWork = TRUE)
} else {
  normalizePath("..", winslash = "/", mustWork = TRUE)
}
knitr::opts_knit$set(root.dir = repo_root)

emit_mermaid <- function(path, theme_path = "website/figures/common_theme.mmd") {
  mermaid_lines <- readLines(path, warn = FALSE)

  if (!is.null(theme_path)) {
    mermaid_lines <- c(
      readLines(theme_path, warn = FALSE),
      mermaid_lines
    )
  }

  cat(
    "```{mermaid}\n",
    paste(mermaid_lines, collapse = "\n"),
    "\n```\n",
    sep = ""
  )
}

emit_yaml_template_entry <- function(path, key) {
  yaml_lines <- readLines(path, warn = FALSE)
  start_line <- which(startsWith(yaml_lines, paste0(key, ":")))

  if (length(start_line) != 1) {
    stop("Expected exactly one top-level YAML entry named '", key, "' in ", path, call. = FALSE)
  }

  top_level_lines <- grep("^[[:alnum:]_.-]+:", yaml_lines)
  next_top_level_line <- top_level_lines[top_level_lines > start_line][1]
  end_line <- if (is.na(next_top_level_line)) length(yaml_lines) else next_top_level_line - 1
  entry_lines <- yaml_lines[start_line:end_line]

  while (length(entry_lines) > 0 && identical(utils::tail(entry_lines, 1), "")) {
    entry_lines <- utils::head(entry_lines, -1)
  }

  cat(
    "```{.yaml filename=\"YAML\"}\n",
    paste(entry_lines, collapse = "\n"),
    "\n```\n",
    sep = ""
  )
}

github_repo <- "https://github.com/koefoeden/multiomeR/tree/main"
force_recreate_graph <- identical(
  tolower(Sys.getenv("FORCE_RECREATE_GRAPH", "false")),
  "true"
)

pipeline_github_file <- file.path(github_repo, "_targets.R")
datasets_config_file <- "cfg_datasets.yaml"
aggregations_config_file <- "cfg_aggregations.yaml"
pipeline_parameters_file <- "cfg_pipeline_parameters.tsv"
datasets_github_file <- file.path(github_repo, "cfg_datasets.yaml")
reactions_github_file <- file.path(github_repo, "cfg_reactions.tsv")
aggregations_github_file <- file.path(github_repo, "cfg_aggregations.yaml")
pipeline_parameters_github_file <- file.path(github_repo, "cfg_pipeline_parameters.tsv")

module_cfg_template_file <- switch(
  pipeline_name,
  "differential_analyses" = "module_differential_analyses/cfg_template.yaml",
  "genetic_enrichment" = "module_genetic_enrichment/cfg_template.yaml",
  NULL
)

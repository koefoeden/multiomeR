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

github_repo <- "https://github.com/koefoeden/multiomeR/tree/main"
force_recreate_graph <- identical(
  tolower(Sys.getenv("FORCE_RECREATE_GRAPH", "false")),
  "true"
)

pipeline_github_file <- file.path(github_repo, "_targets.R")
datasets_template_file <- "cfg_datasets_template.yaml"
datasets_github_file <- file.path(github_repo, datasets_template_file)
reactions_github_file <- file.path(github_repo, "cfg_reactions_template.tsv")

module_cfg_template_file <- switch(
  pipeline_name,
  "differential_analyses" = "module_differential_analyses/cfg_template.yaml",
  "genetic_enrichment" = "module_genetic_enrichment/cfg_template.yaml",
  NULL
)

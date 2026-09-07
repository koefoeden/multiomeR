multiomeR_project_root <- rprojroot::find_root(
  rprojroot::has_file("pixi.toml"),
  path = getwd()
)

source_project_file <- function(path, envir = parent.frame()) {
  sys.source(file.path(multiomeR_project_root, path), envir = envir)
}

require_reference_version <- function(package, expected) {
  observed <- as.character(utils::packageVersion(package))
  if (!identical(observed, expected)) {
    testthat::fail(paste0(
      "Reference tests require ", package, " ", expected,
      "; found ", observed, "."
    ))
  }
  invisible(observed)
}

load_project_test_runtime <- function() {
  if (!exists("load_project_runtime", mode = "function", inherits = TRUE)) {
    source_project_file("R/bootstrap_helpers.R", envir = globalenv())
  }
  load_project_runtime()
}

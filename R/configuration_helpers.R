#' Locate a file in this checkout's selected configuration directory
#'
#' The optional, untracked configuration.local contains one directory path.
#' Relative directory paths are anchored at the repository root. Paths inside
#' configuration files retain their existing interpretation.
#' @param filename Configuration basename, including its extension.
#' @param project_root Repository root; explicit roots are useful for validation.
#' @param must_exist Whether the requested file must exist. Disabled modules
#'   may omit their configuration file.
#' @return Absolute path to the selected configuration file.
configuration_path <- function(filename, project_root = get_project_root(), must_exist = TRUE) {
  selection_file <- file.path(project_root, "configuration.local")
  directory <- "configuration"
  if (file.exists(selection_file)) {
    directory <- trimws(readLines(selection_file, warn = FALSE))
    if (length(directory) != 1L || !nzchar(directory)) {
      stop(selection_file, " must contain exactly one non-empty directory path.", call. = FALSE)
    }
  }
  directory <- path.expand(directory)
  if (!startsWith(directory, "/")) directory <- file.path(project_root, directory)
  if (!dir.exists(directory)) {
    stop("Selected configuration directory does not exist: ", directory, call. = FALSE)
  }
  path <- file.path(normalizePath(directory, winslash = "/", mustWork = TRUE), filename)
  if (must_exist && !file.exists(path)) {
    stop("Missing file in selected configuration directory: ", path, call. = FALSE)
  }
  path
}

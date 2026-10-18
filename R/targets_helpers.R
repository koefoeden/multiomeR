# Add packages to the default target packages.
w_def <- function(packages) {
  targets::tar_option_get("packages") %>%
    c(packages) %>%
    unique()
}

#' Open a BPCells directory with a signature of its files
#'
#' An opened BPCells object records only its path, dimensions and names, so
#' targets hashes it identically after the directory's content changes and
#' skips its dependents. The attached signature of file sizes and modification
#' times, which targets itself uses to decide whether to rehash a file, changes
#' whenever the directory is rewritten. The opening target reruns only after
#' targets detects new content, so dependents rerun exactly when it changes.
#' BPCells records the resolved absolute path; the object keeps `path` as
#' given instead, normally relative to the project root, so objects built from
#' it still open after the store moves.
#'
#' @param path BPCells matrix or fragments directory.
#' @param open Function that opens `path`, such as a BPCells opener or
#'   `open_BPCells_counts_as_double()`.
#' @return The opened object with a `file_signature` attribute.
#' @keywords internal
open_BPCells_dir <- function(path, open = BPCells::open_matrix_dir) {
  files <- sort(list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE))
  info <- file.info(file.path(path, files))
  object <- open(path)
  # A converting opener, such as open_BPCells_counts_as_double(), wraps the directory object.
  if (methods::.hasSlot(object, "dir")) object@dir <- path else object@matrix@dir <- path
  attr(object, "file_signature") <- digest::digest(
    paste(files, info$size, format(as.numeric(info$mtime), digits = 15), sep = ":", collapse = ";"),
    algo = "xxhash64", serialize = FALSE
  )
  object
}

#' Open BPCells counts as doubles
#'
#' GEX counts are stored as `uint32_t`; reading them as doubles gives every
#' consumer the values it got when they were stored as doubles.
#'
#' @param path BPCells matrix directory.
#' @return The lazily converted matrix.
#' @keywords internal
open_BPCells_counts_as_double <- function(path) {
  BPCells::convert_matrix_type(BPCells::open_matrix_dir(path), "double")
}

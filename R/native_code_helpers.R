#' Compile and load a standalone native source file once per R session
#'
#' The library name includes the source checksum, so an unchanged source is
#' reused from the loaded DLLs and an edited one is rebuilt. Each process builds
#' in its own temporary directory to avoid concurrent writes.
#' @param native_source_file C or C++ source compiled with `R CMD SHLIB`.
#' @param library_prefix Prefix of the loaded library name.
#' @return Name of the loaded DLL, for use as the `PACKAGE` of `.Call()`.
load_native_library <- function(native_source_file, library_prefix) {
  library_name <- paste0(library_prefix, "_", substr(unname(tools::md5sum(native_source_file)), 1, 12))
  if (library_name %in% names(getLoadedDLLs())) {
    return(library_name)
  }
  build_dir <- tempfile(paste0(library_prefix, "_"))
  dir.create(build_dir)
  build_source_file <- file.path(build_dir, basename(native_source_file))
  stopifnot(file.copy(native_source_file, build_source_file))
  shared_library_file <- file.path(build_dir, paste0(library_name, .Platform$dynlib.ext))
  build_result <- processx::run(
    command = file.path(R.home("bin"), "R"),
    args = c("CMD", "SHLIB", "-o", shared_library_file, build_source_file),
    wd = build_dir,
    error_on_status = FALSE
  )
  if (build_result$status != 0L) {
    stop("Could not compile ", native_source_file, ":\n", build_result$stdout, build_result$stderr, call. = FALSE)
  }
  dyn.load(shared_library_file)[["name"]]
}

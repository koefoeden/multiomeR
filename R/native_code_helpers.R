#' Compile and load a standalone native source file once per user
#'
#' The library name hashes the source, the R build and the build variables, so
#' an unchanged library is reused from the loaded DLLs or from the user's build
#' cache and an edited one is rebuilt. Each process builds in its own temporary
#' directory and installs the library by an atomic rename, so concurrent workers
#' never load a partial build and later workers skip the compiler.
#' @param native_source_file C or C++ source compiled with `R CMD SHLIB`.
#' @param library_prefix Prefix of the loaded library name.
#' @param env Named character vector of extra build variables, such as `PKG_CXXFLAGS`.
#' @return Name of the loaded DLL, for use as the `PACKAGE` of `.Call()`.
load_native_library <- function(native_source_file, library_prefix, env = character()) {
  library_name <- paste0(library_prefix, "_", substr(rlang::hash(list(
    unname(tools::md5sum(native_source_file)), R.version.string, R.home(), env
  )), 1, 12))
  if (library_name %in% names(getLoadedDLLs())) {
    return(library_name)
  }
  cache_dir <- file.path(tools::R_user_dir("multiomeR", "cache"), "native_libraries")
  shared_library_file <- file.path(cache_dir, paste0(library_name, .Platform$dynlib.ext))
  if (!file.exists(shared_library_file)) {
    build_dir <- tempfile(paste0(library_prefix, "_"))
    dir.create(build_dir)
    build_source_file <- file.path(build_dir, basename(native_source_file))
    stopifnot(file.copy(native_source_file, build_source_file))
    built_library_file <- file.path(build_dir, basename(shared_library_file))
    build_result <- processx::run(
      command = file.path(R.home("bin"), "R"),
      args = c("CMD", "SHLIB", "-o", built_library_file, build_source_file),
      wd = build_dir,
      env = c("current", env),
      error_on_status = FALSE
    )
    if (build_result$status != 0L) {
      stop("Could not compile ", native_source_file, ":\n", build_result$stdout, build_result$stderr, call. = FALSE)
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    staged_library_file <- tempfile(library_name, tmpdir = cache_dir)
    stopifnot(file.copy(built_library_file, staged_library_file), file.rename(staged_library_file, shared_library_file))
  }
  dyn.load(shared_library_file)[["name"]]
}

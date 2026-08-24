#' With annotation hub cache lock
#'
#' Evaluate code while holding a file lock around the AnnotationHub cache.
#'
#' @param expr Expression evaluated in the caller environment, for example while holding a file lock or suppressing matching warnings.
#' @param timeout_ms Milliseconds to wait for the AnnotationHub cache lock before failing.
#' @return The value of `expr`; concurrent callers wait for the cache lock up to
#'   `timeout_ms`.
#' @keywords internal

with_annotation_hub_cache_lock <- function(expr, timeout_ms = 60 * 60 * 1000) {
  cache_dir <- file.path(targets::tar_config_get("store"), "files", "AnnotationHub")
  fs::dir_create(cache_dir)

  lock_path <- file.path(cache_dir, "annotationhub-cache.lock")
  lock_unsupported <- FALSE
  lock <- tryCatch(
    filelock::lock(lock_path, exclusive = TRUE, timeout = timeout_ms),
    error = \(err) {
      if (!grepl("No locks available", conditionMessage(err), fixed = TRUE)) {
        stop(err)
      }
      lock_unsupported <<- TRUE
      NULL
    }
  )
  if (is.null(lock) && lock_unsupported) {
    lock_dir <- paste0(lock_path, ".dir")
    timeout_at <- Sys.time() + timeout_ms / 1000
    while (!dir.create(lock_dir, showWarnings = FALSE)) {
      lock_age_secs <- difftime(Sys.time(), file.info(lock_dir)$mtime, units = "secs")
      if (!is.na(lock_age_secs) && lock_age_secs > timeout_ms / 1000) {
        unlink(lock_dir, recursive = TRUE)
        next
      }
      if (Sys.time() > timeout_at) {
        stop("Timed out waiting for AnnotationHub cache lock: ", lock_path, call. = FALSE)
      }
      Sys.sleep(1)
    }
    on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
    return(force(expr))
  }
  if (is.null(lock)) {
    stop("Timed out waiting for AnnotationHub cache lock: ", lock_path, call. = FALSE)
  }

  on.exit(filelock::unlock(lock), add = TRUE)
  force(expr)
}

#' With annotation hub cache lock
#'
#' Evaluate code while holding a lock around the AnnotationHub cache. The lock
#' is a directory, because creating one is atomic on NFS too: file locks there
#' can fail for one caller while succeeding for another, so callers that fell
#' back to a different lock ran together and read partly downloaded files.
#'
#' @param expr Expression evaluated in the caller environment, for example while holding a file lock or suppressing matching warnings.
#' @param timeout_ms Milliseconds to wait for the AnnotationHub cache lock before failing; an
#'   older lock is taken to be left by a crashed caller and removed.
#' @return The value of `expr`; concurrent callers wait for the cache lock up to
#'   `timeout_ms`.
#' @keywords internal

with_annotation_hub_cache_lock <- function(expr, timeout_ms = 60 * 60 * 1000) {
  cache_dir <- file.path(targets::tar_config_get("store"), "files", "AnnotationHub")
  fs::dir_create(cache_dir)

  lock_dir <- file.path(cache_dir, "annotationhub-cache.lock.dir")
  timeout_at <- Sys.time() + timeout_ms / 1000
  while (!dir.create(lock_dir, showWarnings = FALSE)) {
    lock_age_secs <- difftime(Sys.time(), file.info(lock_dir)$mtime, units = "secs")
    if (!is.na(lock_age_secs) && lock_age_secs > timeout_ms / 1000) {
      unlink(lock_dir, recursive = TRUE)
      next
    }
    if (Sys.time() > timeout_at) {
      stop("Timed out waiting for AnnotationHub cache lock: ", lock_dir, call. = FALSE)
    }
    Sys.sleep(1)
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  force(expr)
}

#' Load the peak-gene-specific batched Kenward-Roger kernel
#'
#' Each worker compiles at most once per source revision into its own temporary
#' cache, avoiding concurrent writes to shared build artifacts.
load_peak_gene_KR_kernel <- function(native_source_file) {
  environment <- new.env(parent = baseenv())
  Rcpp::sourceCpp(
    native_source_file, env = environment,
    cacheDir = file.path(tempdir(), "peak_gene_KR"), showOutput = FALSE
  )
  environment$peak_gene_KR_batch_cpp
}

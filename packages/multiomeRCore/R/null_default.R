#' Default for NULL values
#'
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

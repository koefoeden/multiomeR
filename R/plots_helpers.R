ignore_aes_in_color_legend <- function() {
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
}


symmetric_limits <- function(x) {
  max <- max(abs(x))
  c(-max, max)
}


get_num_facet_rows <- function(in_plot) {
  if (!inherits(in_plot, "ggplot")) {
    return(1)
  }
  n_facet_rows <- ggplot2::ggplot_build(in_plot)$layout$layout$ROW %>% unique() %>% length()
  return(n_facet_rows)
}

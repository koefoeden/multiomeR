CONST_UMAP_ggplot2_theme <- ggplot2::theme(
  axis.ticks = ggplot2::element_blank(),
  axis.text = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank(),
  panel.grid = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(hjust = 0.5),
  axis.line.x = ggplot2::element_line(color = "black"),
  axis.line.y = ggplot2::element_line(color = "black"),
  legend.position = "none"
)


#' Get n scaled alpha size aesthetics
#'
#' Map the number of plotted observations to alpha and point-size values for dense scatter plots.
#'
#' @param data Data frame or tibble whose row count is used as the plot density.
#'   Zero-row inputs are rejected.
#' @param n_small Row-count lower anchor. Inputs at or below this count use
#'   `alpha_max` and `size_max`.
#' @param n_large Row-count upper anchor. Inputs at or above this count use
#'   `alpha_min` and `size_min`.
#' @param alpha_min,alpha_max Minimum and maximum point alpha values returned
#'   after log10 scaling by row count.
#' @param size_min,size_max Minimum and maximum point sizes returned after
#'   log10 scaling by row count.
#' @return List with numeric `alpha` and `size` entries for ggplot point layers.
#' @keywords internal

get_n_scaled_alpha_size_aesthetics <- function(data, n_small = 5000, n_large = 5e5, alpha_min = 0.5, alpha_max = 1, size_min = 0.05, size_max = 0.2) {
  n <- nrow(data)
  if (n == 0) {
    stop("Data has zero rows.")
  }

  # Compute logarithms for scaling
  lns <- log10(n_small)
  lnl <- log10(n_large)
  lnn <- log10(n)

  # Ensure lnn is within the defined bounds
  lnn <- pmin(pmax(lnn, lns), lnl)

  # Calculate alpha and size based on the number of observations
  alpha <- alpha_max - (alpha_max - alpha_min) * (lnn - lns) / (lnl - lns)
  size <- size_max - (size_max - size_min) * (lnn - lns) / (lnl - lns)

  return(list(alpha = alpha, size = size))
}

get_BPCells_plot_embedding_aesthetics <- function(data, rasterize = TRUE) {
  get_n_scaled_alpha_size_aesthetics(
    data = data,
    size_min = if (rasterize) 1 else 0.2,
    size_max = if (rasterize) 2 else 0.6
  )
}

apply_alpha_to_plot_embedding <- function(plot, alpha) {
  if (inherits(plot, "ggplot") && length(plot$layers) > 0) {
    plot$layers[[1]]$aes_params$alpha <- alpha
  }
  plot
}

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

get_num_facet_cols <- function(in_plot) {
  if (!inherits(in_plot, "ggplot")) {
    return(1)
  } else {
    n_facet_cols <- ggplot2::ggplot_build(in_plot)$layout$layout$COL %>% unique() %>% length()
    return(n_facet_cols)
  }
}


# Currently saves plot(s) to a single file.

# reduces the data footprint of a ggplot object through converting to a grob and then back to a ggplot
# prevents quosures from being captured before serialisation.

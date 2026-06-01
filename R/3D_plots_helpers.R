#' Get 3D UMAP widget
#'
#' Build an interactive plotly 3D scatter widget from coordinate columns in a cell-embedding tibble.
#'
#' @param cell_embeddings_tibble Tibble with one row per cell and coordinate columns named from `dim_name` plus a grouping column for hover labels/colors.
#' @param dim_name Coordinate-name prefix; the widget reads `<dim_name>_1`, `<dim_name>_2`, and `<dim_name>_3`.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @return A plotly htmlwidget for interactive inspection; it is not written to disk by this helper.
#' @keywords internal

get_3D_UMAP_widget <- function(cell_embeddings_tibble, dim_name = "PCAharmonyUMAPlocal3D", cluster_col = "PCA_harmony_SNN_cluster") {
  dim_names <- stringr::str_c(dim_name, c(1, 2, 3), sep = "_")

  dim_args <- list(
    title = "",
    autorange = TRUE,
    showspikes = FALSE,
    showgrid = TRUE,
    zeroline = FALSE,
    showline = FALSE,
    autotick = TRUE,
    ticks = "",
    showticklabels = FALSE
  )

  html_widget <- plotly::plot_ly(
    data = cell_embeddings_tibble,
    x = ~ .data[[dim_names[1]]],
    y = ~ .data[[dim_names[2]]],
    z = ~ .data[[dim_names[3]]],
    color = ~ .data[[cluster_col]],
    alpha = 1,
    colors = viridis::turbo(length(unique(cell_embeddings_tibble[[cluster_col]]))),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 1),
    text = ~ .data[[cluster_col]],
    hoverinfo = "text"
  ) %>%
    plotly::layout(showlegend = FALSE, scene = list(xaxis = dim_args, yaxis = dim_args, zaxis = dim_args))

  return(html_widget)
}

#' Run 3D UMAP from embedding matrix
#'
#' Run a three-component cosine UMAP on selected embedding dimensions and return coordinates as a cell tibble.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param n_neighbors UMAP neighbor count; clipped below the number of input cells where the helper does that internally.
#' @param min_dist UMAP minimum-distance parameter controlling how tightly local neighborhoods are packed.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param col_prefix Prefix assigned to generated coordinate columns, for example `LSI_UMAP` gives `LSI_UMAP_1` and `LSI_UMAP_2`.
#' @return A tibble keyed by cell barcode with three UMAP coordinate columns.
#' @keywords internal

run_3D_UMAP_from_embedding_matrix <- function(
  embedding_matrix,
  dims,
  n_neighbors,
  min_dist,
  seed = 1,
  dim_prefix = "PCA_",
  col_prefix = "PCAharmonyUMAPlocal3D"
) {
  run_UMAP_from_embedding_matrix(
    embedding_matrix = embedding_matrix,
    dims = dims,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    seed = seed,
    dim_prefix = dim_prefix,
    col_prefix = col_prefix,
    n_components = 3
  ) |>
    embedding_matrix_to_tibble()
}

generate_static_images_from_HTML_widget <- function(angles_deg_vec, html_widget) {
  out_dir <- get_structured_file_path(override_suffix = "")

  output_files <- angles_deg_vec %>%
    purrr::set_names() %>%
    purrr::imap(
      ~ {
        angle_rad <- .x * pi / 180 # Convert degrees to radians
        camera <- list(eye = list(x = 2 * cos(angle_rad), y = 2 * sin(angle_rad), z = 1.25)) # Adjust Z to change the elevation angle
        p_camera <- plotly::layout(html_widget, scene = list(camera = camera))
        out_file <- stringr::str_glue("{out_dir}/frame_{.y}.png")
        plotly::save_image(p = p_camera, file = out_file, width = 1600, height = 1600)
        return(out_file)
      }
    ) %>%
    unlist()

  return(output_files)
}

generate_GIF <- function(image_paths_vec) {
  out_file <- get_structured_file_path(filetype = "gif")

  gifski::gifski(
    png_files = image_paths_vec,
    gif_file = out_file,
    width = 1600,
    height = 1600,
    delay = 12 / length(image_paths_vec)
  )

  return(out_file)
}

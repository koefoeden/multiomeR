repo_root <- rprojroot::find_root(rprojroot::has_file("_targets.R"))
setwd(repo_root)

source("R/bootstrap_helpers.R")
load_project_runtime()

library(patchwork)

figure_text_family <- "Nimbus Sans"
ggplot2::theme_update(text = ggplot2::element_text(family = figure_text_family))

minimal_UMAP_theme <- ggplot2::theme(
  axis.title = ggplot2::element_blank(),
  axis.text = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  plot.title = ggplot2::element_blank()
)

figure_text_theme <- ggplot2::theme(
  text = ggplot2::element_text(family = figure_text_family)
)

panel_title_theme <- ggplot2::theme(
  text = ggplot2::element_text(family = figure_text_family),
  plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
  plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
  plot.title.position = "plot"
)

set_text_layer_family <- function(plot) {
  if (inherits(plot, "ggplot")) {
    for (layer_index in seq_along(plot$layers)) {
      if (grepl("Text|Label", class(plot$layers[[layer_index]]$geom)[[1]])) {
        plot$layers[[layer_index]]$aes_params$family <- figure_text_family
      }
    }
  }
  if (inherits(plot, "patchwork")) {
    for (plot_index in seq_along(plot$patches$plots)) {
      plot$patches$plots[[plot_index]] <- set_text_layer_family(plot$patches$plots[[plot_index]])
    }
  }
  plot
}

read_minimal_UMAP <- function(plot_paths, pick_regex, title, show_legend = FALSE) {
  image_path <- pick_plot_path(plot_paths, pick_regex = pick_regex)
  plot <- readRDS(plot_object_path(image_path)) +
    minimal_UMAP_theme +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8, face = "bold"),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )
  if (show_legend) {
    plot + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, title = NULL))
  } else {
    plot + ggplot2::guides(colour = "none")
  }
}

figure_2_UMAP_panel <- function() {
  GEX_UMAP <- read_minimal_UMAP(
    targets::tar_read(harmony.categorical.UMAPs.GEX.mixed_human_7x),
    "PCA_harmony_SNN_cluster_cell_type",
    "GEX"
  )
  ATAC_UMAP <- read_minimal_UMAP(
    targets::tar_read(categorical.UMAPs.ATAC.mixed_human_7x),
    "LSI_harmony_SNN_cluster_cell_type",
    "ATAC"
  )
  WNN_UMAP <- read_minimal_UMAP(
    targets::tar_read(categorical.UMAPs.WNN.mixed_human_7x),
    "WNN_harmony_SNN_cluster_cell_type",
    "WNN",
    show_legend = TRUE
  )

  UMAPs <- (GEX_UMAP | patchwork::plot_spacer() | ATAC_UMAP | patchwork::plot_spacer() | WNN_UMAP) +
    patchwork::plot_layout(widths = c(1, 0.08, 1, 0.08, 1))

  (patchwork::guide_area() / UMAPs) +
    patchwork::plot_layout(heights = c(0.12, 1), guides = "collect") &
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = ggplot2::element_text(size = 7),
      legend.key.size = grid::unit(3, "mm")
    ) &
    figure_text_theme
}

figure_2_TF_activity_panel <- function(plot) {
  plot +
    ggplot2::labs(
      title = "Transcription factor activity",
      x = NULL,
      y = "Cell type",
      fill = "Mean Z-score"
    ) +
    ggplot2::scale_x_discrete(labels = \(labels) sub("_+MA.*$", "", labels)) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title.position = "top",
        barheight = grid::unit(18, "mm"),
        barwidth = grid::unit(3, "mm")
      )
    ) +
    ggplot2::theme(
      legend.key.size = grid::unit(2.5, "mm"),
      legend.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 6),
      legend.spacing.y = grid::unit(0, "mm"),
      legend.box.spacing = grid::unit(0, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    figure_text_theme +
    panel_title_theme
}

figure_2_peak_gene_link_panel <- function(plot, group = "B") {
  plot$data <- plot$data[plot$data$group == group, , drop = FALSE]
  plot$data$group <- droplevels(plot$data$group)
  plot$patches$layout$heights <- grid::unit(c(1, 1, 1), "null")
  plot$patches$plots[[1]]$layers[[2]]$aes_params$size <- 1.6
  plot <- set_text_layer_family(plot)

  plot +
    ggplot2::labs(x = NULL, y = "Insertions") +
    patchwork::plot_annotation(
      title = "Peak-gene links",
      theme = panel_title_theme
    ) &
    ggplot2::guides(colour = "none") &
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) &
    figure_text_theme
}

figure_2_single_cell_GWAS_panel <- function(plot) {
  set_text_layer_family(plot) +
    minimal_UMAP_theme +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.18))) +
    ggplot2::labs(
      title = "Single-nucleus GWAS enrichment",
      subtitle = "Monocyte count"
    ) +
    figure_text_theme +
    panel_title_theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8.5, margin = ggplot2::margin(b = 3)),
      plot.margin = ggplot2::margin(2, 5, 2, 5),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7)
    )
}

strip_GWAS_author_year_suffix <- function(GWAS_ID) {
  sub("_[^_]+[0-9]{4}$", "", as.character(GWAS_ID)) |>
    gsub("([a-z])([A-Z])", "\\1 \\2", x = _) |>
    sub("^BCell", "B-cell", x = _)
}

figure_2_GWAS_heatmap_panel <- function(plot) {
  score_data <- plot$patches$plots[[2]]$data
  required_columns <- c("cluster", "GWAS_ID", "relative_deviation", "support_label")
  missing_columns <- setdiff(required_columns, colnames(score_data))
  if (!is.data.frame(score_data) || length(missing_columns)) {
    stop(
      "Unexpected chromVAR deviation heatmap plot-object structure; missing score columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  score_plot <- ggplot2::ggplot(
    score_data,
    ggplot2::aes(x = .data$cluster, y = .data$GWAS_ID, fill = .data$relative_deviation)
  ) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.15) +
    ggplot2::geom_text(ggplot2::aes(label = .data$support_label), size = 2.8, na.rm = TRUE) +
    ggplot2::scale_y_discrete(
      drop = FALSE,
      labels = strip_GWAS_author_year_suffix
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      name = "Relative deviation",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(18, "mm"),
        barheight = grid::unit(3, "mm")
      )
    ) +
    ggplot2::coord_fixed(ratio = 0.62) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7, base_family = figure_text_family) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, vjust = 1, size = 6.5),
      axis.text.y = ggplot2::element_text(size = 6.5),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "top",
      plot.margin = ggplot2::margin(2, 3, 1, 2)
    )

  patchwork::wrap_plots(score_plot) +
    patchwork::plot_annotation(
      title = "Cell-type GWAS enrichment",
      theme = panel_title_theme
    ) &
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 6),
      legend.key.size = grid::unit(2.5, "mm"),
      legend.box.spacing = grid::unit(0, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    ) &
    figure_text_theme
}

output_dir <- "manuscript_figures/outputs/combined_targets_output_figures_7x"
benchmark_plot_rds <- file.path(
  "manuscript_figures",
  "outputs",
  "benchmark",
  "multimodal_seurat_walltime_plot.rds"
)
if (!file.exists(benchmark_plot_rds)) {
  stop(
    "Missing benchmark plot object: ", benchmark_plot_rds, "\n",
    "Run manuscript_figures/render_benchmark_walltime_plot.R first.",
    call. = FALSE
  )
}

panel_specs <- tibble::tribble(
  ~figure_number , ~tag , ~plot_call                                                                                                                                                                       , ~pick_regex                  , ~title        , ~caption , ~plot_modifier                                                   ,
  # "S1"           , "A"  , quote(targets::tar_read(categorical.UMAPs.WNN.mixed_human_7x))                                                      , "WNN_harmony_SNN_cluster_cell_type" , NA_character_ , "UMAP embedding derived from native RNA+ATAC weighted-nearest-neighbor graph. Nuclei are colored by cluster-level cell type labels." , NA                                          ,
  # "S1"           , "B"  , quote(targets::tar_read(markers_dot_plot.GEX.mixed_human_7x))                                                       , NA_character_                       , NA_character_ , "Expression of canonical marker genes across cell types."                                                                            , quote(plot + ggplot2::labs(y = "celltype")) ,
  # "S1"           , "C"  , quote(targets::tar_read(coverage_tracks_plots.ATAC.mixed_human_7x))                                                 , "MS4A1"                             , NA_character_ , "Chromatin accessibility track at the MS4A1 locus."                                                                                  , NA                                          ,
  "2"            , "A"  , quote(figure_2_UMAP_panel())                                                                                                                                                     , NA_character_                , NA_character_ , NA       , NA                                                               ,
  "2"            , "B"  , quote(targets::tar_read(motif_family_accessibility_heatmap.ATAC.mixed_human_7x))                                                                                                , NA_character_                , NA_character_ , NA       , quote(figure_2_TF_activity_panel(plot))                          ,
  "2"            , "C"  , quote(targets::tar_read(peak_gene_correlation_top_link_ATAC_tracks_plots.peak_gene_correlation.WNN.mixed_human_7x))                                                              , "Brain_rank001_"             , NA_character_ , NA       , quote(figure_2_peak_gene_link_panel(plot, group = "Brain"))      ,
  "2"            , "D"  , quote(targets::tar_read(TRS_UMAPs.WNN_harmony_SNN.SCAVENGE.single_nucleus.genetic_enrichment.mixed_human_7x))                                                                    , "MonocyteCount_Vuckovic2020" , NA_character_ , NA       , quote(figure_2_single_cell_GWAS_panel(plot))                     ,
  "2"            , "E"  , quote(targets::tar_read(chromVAR_deviation_heatmap.cell_type_pseudobulk.genetic_enrichment.mixed_human_7x))                                                                      , NA_character_                , NA_character_ , NA       , quote(figure_2_GWAS_heatmap_panel(plot))                         ,
  "3"            , "A"  , quote(readRDS(benchmark_plot_rds))                                                                                                                                             , NA_character_                , NA_character_ , NA       , NA                                                               ,
)

figure_specs <- tibble::tribble(
  ~figure_number , ~height , ~layout                      , ~caption      ,
  # "S1"           ,       9 , quote((A / B / C))           , "Visualization of the observed cell types and expected markers in the PBMC dataset." ,
  "2"            , 7.25    , quote(A / (B | C) / (D | E)) , NA_character_ ,
  "3"            , 2.5     , quote(A)                     , NA_character_ ,
)

expected_pngs <- file.path(output_dir, paste0(figure_specs$figure_number, ".png"))
stale_pngs <- setdiff(list.files(output_dir, pattern = "[.]png$", full.names = TRUE), expected_pngs)
unlink(stale_pngs)

render_figures(
  panel_specs = panel_specs,
  figure_specs = figure_specs,
  output_dir = output_dir,
  figure_theme = figure_text_theme
)

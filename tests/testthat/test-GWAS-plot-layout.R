source_project_file("R/GWAS_plot_helpers.R")

make_GWAS_layout_fixture <- function() {
  metadata <- tibble::tibble(
    GWAS_ID = factor(c("Trait_A", "Trait_B", "Trait_C"), levels = c("Trait_C", "Trait_B", "Trait_A")),
    Category = c("A", "B", "B"),
    finemappingMethod = c("PICS", "SuSiE-inf", "PICS"),
    n_credible_set_loci = c(3, 100, 1000),
    sample_size = c(1000, 100000, 3000000)
  )
  scores <- tidyr::expand_grid(GWAS_ID = as.character(metadata$GWAS_ID), cluster = c("Small", "Large")) |>
    dplyr::left_join(metadata |> dplyr::mutate(GWAS_ID = as.character(GWAS_ID)), by = "GWAS_ID") |>
    dplyr::mutate(
      median_score = seq_len(dplyr::n()),
      n_cells = rep(c(71, 30000), 3),
      significance = factor(NA_character_, levels = c("P <= 0.05", "P < 0.01"))
    )
  list(metadata = metadata, scores = scores)
}

testthat::test_that("tiles, annotation bars and empty dotplots share exact row extents", {
  fixture <- make_GWAS_layout_fixture()
  tracks <- plot_GWAS_metadata_tracks(fixture$metadata)
  for (point_size_col in list(NULL, "significance")) {
    plot <- plot_GWAS_by_cluster_heatmap(fixture$scores, tracks, point_size_col = point_size_col)
    for (i in seq_len(5)) {
      panels <- ggplot2::ggplot_build(plot[[i]])$layout$panel_params
      for (panel in panels) testthat::expect_equal(panel$y.range, c(0.5, 3.5))
    }
    grob <- patchwork::patchworkGrob(plot)
    body <- grob$layout[grepl("^panel.*-[1-5]$", grob$layout$name), ]
    support <- grob$layout[grob$layout$name == "panel-6", ]
    heatmap <- grob$layout[grob$layout$name == "panel-5", ]
    testthat::expect_equal(nrow(body), 5)
    testthat::expect_length(unique(body$t), 1)
    testthat::expect_length(unique(body$b), 1)
    testthat::expect_equal(support$l, heatmap$l)
    testthat::expect_equal(support$r, heatmap$r)
    testthat::expect_equal(sum(grob$layout$name == "guide-box"), 1)
  }
})

testthat::test_that("nuclei labels clear bars and support can be omitted", {
  fixture <- make_GWAS_layout_fixture()
  tracks <- plot_GWAS_metadata_tracks(fixture$metadata)
  plot <- plot_GWAS_by_cluster_heatmap(fixture$scores, tracks)
  support <- ggplot2::ggplot_build(plot[[6]])
  testthat::expect_true(all(support$data[[2]]$hjust < 0))
  testthat::expect_equal(support$data[[2]]$y, support$data[[1]]$y)
  testthat::expect_gt(support$layout$panel_params[[1]]$y.range[[2]], max(fixture$scores$n_cells))
  without_support <- plot_GWAS_by_cluster_heatmap(fixture$scores, tracks, show_feature_support = FALSE)
  testthat::expect_s3_class(patchwork::patchworkGrob(without_support), "gtable")
  testthat::expect_s3_class(plot_GWAS_by_cluster_heatmap(fixture$scores[0, ], tracks), "empty_plot_list")
})

source_project_file("R/processing_helpers.R")
source_project_file("R/QC_helpers.R")

testthat::test_that("QC violin thresholds retain their GEM well association", {
  thresholds <- get_GEM_well_QC_exclude_threshold_tibble(
    GEM_well_QC_exclude_list = list(
      well_1 = "nCount_RNA < 5",
      well_2 = "nCount_RNA > 35",
      well_without_cutoff = NULL
    ),
    feature_names = "nCount_RNA"
  )

  testthat::expect_identical(thresholds$GEM_well_ID, c("well_1", "well_2"))
  testthat::expect_identical(thresholds$threshold, c(5, 35))
  testthat::expect_identical(thresholds$ymin, c(-Inf, 35))
  testthat::expect_identical(thresholds$ymax, c(5, Inf))
})

testthat::test_that("QC violin cutoff shading is confined to each GEM well", {
  plots <- plot_per_dataset_QC_violins(
    metadata_tibble = tibble::tibble(
      GEM_well_ID = rep(c("well_1", "well_2"), each = 20),
      dataset = rep(c("dataset_1", "dataset_2"), each = 20),
      nCount_RNA = seq_len(40)
    ),
    feature_names = "nCount_RNA",
    GEM_well_QC_exclude_list = list(
      well_1 = "nCount_RNA < 5",
      well_2 = "nCount_RNA > 35"
    ),
    show_dataset_legend = TRUE
  )

  threshold_rects <- plots$nCount_RNA$layers[[1]]$data
  testthat::expect_equal(threshold_rects$xmin, c(0.55, 1.55))
  testthat::expect_equal(threshold_rects$xmax, c(1.45, 2.45))
  testthat::expect_identical(
    plots$nCount_RNA$theme$legend.position,
    "bottom"
  )
  testthat::expect_no_error(ggplot2::ggplot_build(plots$nCount_RNA))
})

testthat::test_that("GEM-well QC comparisons use manifest selection and labels", {
  plots <- plot_GEM_well_QC_comparisons(
    metadata_tibble = tibble::tibble(
      GEM_well_ID = rep(c("well_1", "well_2"), each = 20),
      dataset = rep(c("dataset_1", "dataset_2"), each = 20),
      nCount_RNA = seq_len(40),
      log10_nCount_RNA = log10(seq_len(40)),
      scDblFinder.score_GEX = seq_len(40) / 40
    ),
    QC_metric_manifest_tibble = tibble::tribble(
      ~metric_id, ~display_name, ~description, ~available_from_stage, ~plot_min_q, ~plot_max_q, ~do_plot,
      "nCount_RNA", "RNA UMI count", "Total RNA UMI count.", "GEM_well", NA, 0.99, TRUE,
      "log10_nCount_RNA", "Log10 RNA UMI count", "Log RNA UMI count.", "GEM_well", NA, NA, FALSE,
      "scDblFinder.score_GEX", "GEX doublet score", "Doublet score.", "GEX", NA, NA, TRUE
    ),
    QC_exclude_per_GEM_well_list = list(
      well_1 = "nCount_RNA < 5",
      well_2 = "nCount_RNA > 35"
    )
  )

  testthat::expect_named(plots, "nCount_RNA")
  testthat::expect_identical(plots$nCount_RNA$labels$title, "RNA UMI count")
  testthat::expect_identical(
    plots$nCount_RNA$labels$subtitle,
    "nCount_RNA: Total RNA UMI count."
  )
  testthat::expect_null(plots$nCount_RNA$labels$caption)
  testthat::expect_identical(nrow(plots$nCount_RNA$data), 38L)
  testthat::expect_equal(range(plots$nCount_RNA$data$value), c(1, 39))
  testthat::expect_equal(
    plots$nCount_RNA$layers[[1]]$data$xmin,
    c(0.55, 1.55)
  )
  testthat::expect_no_error(ggplot2::ggplot_build(plots$nCount_RNA))
})

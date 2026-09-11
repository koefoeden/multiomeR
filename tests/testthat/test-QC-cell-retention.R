source_project_file("R/QC_helpers.R")

testthat::test_that("retention reports preserve completely lost and initially empty wells", {
  before <- tibble::tibble(GEM_well_ID = c("kept", "kept", "lost"), barcode_w_prefix = c("a", "b", "c"))
  after <- before[1, ]
  counts <- summarize_QC_cell_retention(before, after, c("kept", "lost", "empty"),
    "checkpoint:1_pre-aggregation-QC", list("First filters" = c("b", "c")))
  testthat::expect_identical(counts$GEM_well_ID, c("kept", "lost", "empty"))
  testthat::expect_identical(counts$input_cells, c(2L, 1L, 0L))
  testthat::expect_identical(counts$retained_cells, c(1L, 0L, 0L))
  testthat::expect_identical(counts$excluded_cells, c(1L, 1L, 0L))
  testthat::expect_equal(counts$retained_fraction, c(0.5, 0, NA_real_))
  cumulative <- summarize_QC_cell_retention(after, after[FALSE, ],
    c("kept", "lost", "empty"), "checkpoint:3_GEX-QC", list("Second filters" = "a"), counts)
  testthat::expect_identical(cumulative[1:3, ], counts)
  testthat::expect_identical(cumulative$stage,
    rep(c("checkpoint:1_pre-aggregation-QC", "checkpoint:3_GEX-QC"), each = 3))
  testthat::expect_identical(cumulative$retained_cells[4:6], c(0L, 0L, 0L))
  testthat::expect_error(summarize_QC_cell_retention(before, after,
    c("kept", "lost", "empty"), "Broken QC", list("First filters" = c("b", "c")), counts), "do not connect")
  testthat::expect_s3_class(ggplot2::ggplot_build(plot_QC_cell_retention(cumulative)), "ggplot_built")
})

testthat::test_that("discard branches count overlapping actions once and conserve nuclei", {
  before <- tibble::tibble(GEM_well_ID = "well", barcode_w_prefix = letters[1:4])
  counts <- summarize_QC_cell_retention(before, before[1, ], "well", "checkpoint:3_GEX-QC",
    list("Called doublets" = c("b", "c"), "High-doublet clusters" = c("c", "d")))
  testthat::expect_identical(counts$input_cells, c(4L, 2L))
  testthat::expect_identical(counts$excluded_cells, c(2L, 1L))
  testthat::expect_identical(counts$retained_cells, c(2L, 1L))
  plot <- plot_QC_cell_retention(counts)
  polygons <- plot$layers[[1]]$data
  discarded <- split(polygons[polygons$outcome == "Discarded", ],
    polygons$flow[polygons$outcome == "Discarded"])
  testthat::expect_length(discarded, 2L)
  testthat::expect_equal(vapply(discarded, function(branch) tail(branch$y, 1) - branch$y[1], numeric(1)),
    c("1 excluded" = 2, "2 excluded" = 1))
  testthat::expect_error(summarize_QC_cell_retention(before, before[1, ], "well", "Broken QC",
    list("Incomplete reasons" = "b")), "do not reproduce")
  empty <- summarize_QC_cell_retention(before[FALSE, ], before[FALSE, ], "well", "Empty QC", list())
  testthat::expect_s3_class(ggplot2::ggplot_build(plot_QC_cell_retention(empty)), "ggplot_built")
})

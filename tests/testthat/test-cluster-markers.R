source_project_file("R/processing_multimodal_helpers.R")
source_project_file("R/processing_GEX_helpers.R")
source_project_file("R/QC_helpers.R")

testthat::test_that("size-ranked cluster IDs preserve membership and barcode order", {
  load_project_test_runtime()
  original <- stats::setNames(factor(c("200", "2", "10", "200", "2", "10", "200", "99"),
    levels = c("99", "200", "10", "2", "unused")), paste0("cell", 1:8))
  ranked <- number_clusters_by_size(original)
  testthat::expect_identical(names(ranked), names(original))
  testthat::expect_identical(as.character(ranked), c("1", "2", "3", "1", "2", "3", "1", "4"))
  testthat::expect_identical(outer(as.character(original), as.character(original), `==`),
    outer(as.character(ranked), as.character(ranked), `==`))
  testthat::expect_identical(number_clusters_by_size(ranked), ranked)
  testthat::expect_identical(filter_clusters_by_min_barcodes(ranked, 2), droplevels(ranked[1:7]))
})

make_cluster_marker_fixture <- function() {
  cluster <- rep(c("2", "10", "4", "5", "6", "99"), each = 6)
  cell_type <- rep(c("Pair", "Pair", "Multiple", "Multiple", "Multiple", "Singleton"), each = 6)
  metadata <- tibble::tibble(
    barcode_w_prefix = paste0("cell", seq_along(cluster)),
    PCA_harmony_SNN_cluster = factor(cluster, levels = c("99", "10", "6", "5", "4", "2", "unused")),
    PCA_harmony_SNN_cluster_cell_type = cell_type
  )
  set.seed(1)
  counts <- matrix(stats::rpois(36 * 30, 4), nrow = 30)
  counts[1, cluster == "2"] <- 40
  counts[2, cluster == "10"] <- 40
  counts[3, cluster == "4"] <- 40
  counts[4, cluster == "5"] <- 40
  counts[5, cluster == "6"] <- 40
  dimnames(counts) <- list(paste0("gene", seq_len(nrow(counts))), metadata$barcode_w_prefix)
  list(metadata = metadata, counts = counts)
}

testthat::test_that("branch records omit singletons and unused factor levels", {
  fixture <- make_cluster_marker_fixture()
  groups <- make_cluster_marker_groups(fixture$metadata)
  testthat::expect_identical(groups$cell_type, c("Multiple", "Pair"))
  testthat::expect_identical(groups$clusters[[2]], c("2", "10"))
  testthat::expect_type(groups$metadata[[2]]$cluster, "character")
  testthat::expect_equal(nrow(make_cluster_marker_groups(
    fixture$metadata[fixture$metadata$PCA_harmony_SNN_cluster_cell_type == "Singleton", ]
  )), 0L)
})

testthat::test_that("two clusters give one directed comparison with a restricted background", {
  fixture <- make_cluster_marker_fixture()
  group <- make_cluster_marker_groups(fixture$metadata)[2, ]
  result <- get_cluster_markers_from_matrix(fixture$counts, group)
  testthat::expect_identical(unique(as.character(result$markers$cluster)), "2")
  testthat::expect_gt(result$markers$avg_log2FC[result$markers$gene == "gene1"], 0)
  testthat::expect_lt(result$markers$avg_log2FC[result$markers$gene == "gene2"], 0)
  expected_background <- rowMeans(fixture$counts[, fixture$metadata$PCA_harmony_SNN_cluster == "10"])
  testthat::expect_equal(result$markers$background_mean, unname(expected_background[result$markers$gene]))
  fixture$counts[, fixture$metadata$PCA_harmony_SNN_cluster_cell_type != "Pair"] <- 100000
  testthat::expect_equal(result, get_cluster_markers_from_matrix(fixture$counts, group))
  testthat::expect_equal(result$markers$p_val_adj, p.adjust(result$markers$p_val_raw, "BH"))
  plot <- plot_cluster_marker_volcano(result)
  testthat::expect_s3_class(plot$facet, "FacetNull")
  testthat::expect_match(plot$labels$subtitle, "Cluster 2 - cluster 10", fixed = TRUE)
  testthat::expect_match(plot$labels$subtitle, "BH-adjusted across genes separately for each contrast", fixed = TRUE)
})

testthat::test_that("multicluster results use pooled within-type backgrounds and per-contrast BH", {
  fixture <- make_cluster_marker_fixture()
  group <- make_cluster_marker_groups(fixture$metadata)[1, ]
  result <- get_cluster_markers_from_matrix(fixture$counts, group)
  testthat::expect_setequal(as.character(result$markers$cluster), c("4", "5", "6"))
  for (cluster in result$clusters) {
    markers <- result$markers[as.character(result$markers$cluster) == cluster, ]
    background <- fixture$metadata$PCA_harmony_SNN_cluster_cell_type == "Multiple" &
      fixture$metadata$PCA_harmony_SNN_cluster != cluster
    means <- rowMeans(fixture$counts[, background])
    testthat::expect_equal(markers$background_mean, unname(means[markers$gene]))
    testthat::expect_equal(markers$p_val_adj, p.adjust(markers$p_val_raw, "BH"))
  }
  plot <- plot_cluster_marker_volcano(result)
  testthat::expect_s3_class(plot$facet, "FacetWrap")
  testthat::expect_equal(nrow(ggplot2::ggplot_build(plot)$layout$layout), 3L)
  testthat::expect_match(plot$labels$subtitle, "pooled remaining clusters within this cell type", fixed = TRUE)
})

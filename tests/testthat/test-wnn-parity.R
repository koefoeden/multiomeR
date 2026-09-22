source_project_file("tests/testthat/fixtures/algorithm-deviations.R")

# Thresholds sit about 0.01 below the values observed with Seurat 5.5.0 at production
# settings, so a real loss of agreement fails rather than passing a loose floor.
testthat::test_that("integration: BPCells WNN meets the stress-fixture similarity thresholds", {
  load_algorithm_test_runtime()
  require_reference_version("Seurat", "5.5.0")
  metrics <- get_wnn_similarity_metrics(
    embeddings = make_wnn_fixture(
      n_cells = 160L,
      n_dimensions = 8L,
      n_clusters = 4L,
      ATAC_cluster_order = c(2L, 1L, 4L, 3L),
      seed = 847
    ),
    k = 15L,
    candidate_k = 50L
  )

  validate_wnn_fixture(
    label = "stress",
    metrics = metrics,
    weight_threshold = 0.92,
    mean_overlap_threshold = 0.98,
    q25_overlap_threshold = 0.95
  )
})

testthat::test_that("integration: BPCells WNN meets production similarity thresholds", {
  load_algorithm_test_runtime()
  require_reference_version("Seurat", "5.5.0")
  embeddings <- make_wnn_fixture(
    n_cells = 400L,
    n_dimensions = 12L,
    n_clusters = 5L,
    ATAC_cluster_order = c(2L, 1L, 4L, 5L, 3L),
    seed = 848
  )
  # k = 20 is the manifest default; k = 50 is the most common configured value.
  thresholds <- list(`20` = c(0.97, 0.98), `50` = c(0.98, 0.99))
  for (k in c(20L, 50L)) {
    validate_wnn_fixture(
      label = paste0("production k = ", k),
      metrics = get_wnn_similarity_metrics(embeddings, k = k, candidate_k = 200L),
      weight_threshold = thresholds[[as.character(k)]][[1]],
      mean_overlap_threshold = thresholds[[as.character(k)]][[2]],
      q25_overlap_threshold = 0.95
    )
  }
})

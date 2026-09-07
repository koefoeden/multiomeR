source_project_file("tests/testthat/fixtures/algorithm-deviations.R")

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
    weight_threshold = 0.80,
    mean_overlap_threshold = 0.85,
    q25_overlap_threshold = 0.75
  )
})

testthat::test_that("integration: BPCells WNN meets production-like similarity thresholds", {
  load_algorithm_test_runtime()
  require_reference_version("Seurat", "5.5.0")
  metrics <- get_wnn_similarity_metrics(
    embeddings = make_wnn_fixture(
      n_cells = 400L,
      n_dimensions = 12L,
      n_clusters = 5L,
      ATAC_cluster_order = c(2L, 1L, 4L, 5L, 3L),
      seed = 848
    ),
    k = 30L,
    candidate_k = 200L
  )

  validate_wnn_fixture(
    label = "production-like",
    metrics = metrics,
    weight_threshold = 0.90,
    mean_overlap_threshold = 0.94,
    q25_overlap_threshold = 0.90
  )
})

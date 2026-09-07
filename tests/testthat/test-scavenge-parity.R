source_project_file("tests/testthat/fixtures/algorithm-deviations.R")

testthat::test_that("integration: SCAVENGE propagation matches independent references", {
  load_algorithm_test_runtime()
  context <- make_SCAVENGE_test_context()
  observed_precise <- run_sparse_random_walk_with_restart(
    NN_graph = context$graph,
    seed_cells = context$seed_cells,
    restart_prob = context$restart_prob,
    stationary_cutoff = 1e-12
  )
  observed_precomputed <- run_sparse_random_walk_with_restart(
    NN_graph = context$graph,
    seed_cells = context$seed_cells,
    restart_prob = context$restart_prob,
    stationary_cutoff = 1e-12,
    transition_matrix = context$transition
  )
  expect_at_most(
    max(abs(observed_precise - observed_precomputed)),
    0,
    "SCAVENGE precomputed-transition delta"
  )

  restart <- setNames(numeric(nrow(context$graph)), rownames(context$graph))
  restart[context$seed_cells] <- 1 / length(context$seed_cells)
  closed_form <- solve(
    diag(nrow(context$graph)) -
      (1 - context$restart_prob) * as.matrix(context$transition),
    context$restart_prob * restart
  )
  expect_at_most(
    max(abs(observed_precise - closed_form)),
    1e-10,
    "SCAVENGE random-walk closed-form delta"
  )

  reference_precise <- reference_SCAVENGE_random_walk(
    context$weighted_graph,
    context$seed_cells,
    context$restart_prob,
    stationary_cutoff = 1e-12
  )
  expect_at_most(
    max(abs(observed_precise - reference_precise)),
    1e-12,
    "SCAVENGE pinned-reference propagation delta"
  )
})

testthat::test_that("integration: SCAVENGE degree-matched sampling matches the reference", {
  load_algorithm_test_runtime()
  context <- make_SCAVENGE_test_context()
  testthat::expect_identical(
    context$seed_idx,
    get_SCAVENGE_seed_index(context$z_score, seed_percent = 0.05)
  )
  testthat::expect_identical(context$reference_samples, context$native_samples)

  singleton_graph <- Matrix::Matrix(
    matrix(
      c(0, 1, 0, 1, 0, 1, 0, 1, 0),
      nrow = 3,
      dimnames = list(paste0("singleton", 1:3), paste0("singleton", 1:3))
    ),
    sparse = TRUE
  )
  singleton_samples <- sample_SCAVENGE_degree_matched_seed_indices(
    singleton_graph,
    setNames(c(FALSE, TRUE, FALSE), rownames(singleton_graph)),
    permutation_times = 10L
  )
  testthat::expect_true(
    all(vapply(singleton_samples, identical, logical(1), 2L))
  )
})

testthat::test_that("integration: SCAVENGE permutation statistics match across cores", {
  load_algorithm_test_runtime()
  context <- make_SCAVENGE_test_context()
  native_source_file <- file.path(
    multiomeR_project_root,
    "src",
    "scavenge_random_walk.cpp"
  )
  native_statistics_1_core <- run_SCAVENGE_permutation_statistics(
    context$transition,
    context$native_samples,
    context$observed_score,
    context$cluster_index_record,
    cores = 1L,
    restart_prob = context$restart_prob,
    native_source_file = native_source_file
  )
  native_statistics_2_cores <- run_SCAVENGE_permutation_statistics(
    context$transition,
    context$native_samples,
    context$observed_score,
    context$cluster_index_record,
    cores = 2L,
    restart_prob = context$restart_prob,
    native_source_file = native_source_file
  )

  testthat::expect_identical(
    native_statistics_1_core,
    native_statistics_2_cores
  )
  testthat::expect_equal(
    unname(native_statistics_1_core$cell_exceedance_counts),
    unname(get_SCAVENGE_reference_exceedance_counts(context))
  )
})

testthat::test_that("integration: SCAVENGE end-to-end outputs match the pinned reference", {
  load_algorithm_test_runtime()
  context <- make_SCAVENGE_test_context()
  reference_exceedance_counts <- get_SCAVENGE_reference_exceedance_counts(context)
  reference_scores <- get_SCAVENGE_reference_scores(context)

  set.seed(431)
  observed_result <- get_SCAVENGE_result_from_chromVAR_z_score_record(
    chromVAR_z_score_record = list(
      GWAS_ID = "fixture",
      z_score_vec = context$z_score
    ),
    NN_graph = context$weighted_graph,
    metadata_tibble = context$metadata_tibble,
    graph_name = "PCA_harmony_SNN",
    cores = 2,
    permutation_times = context$permutation_times,
    restart_prob = context$restart_prob,
    seed_percent = 0.05,
    scale_percent = 0.1
  )
  observed <- observed_result$TRS_tibble
  cluster_summary <- observed_result$TRS_summary_tibble

  testthat::expect_equal(nrow(cluster_summary), 12L)
  testthat::expect_true(
    all(cluster_summary$permutation_times == context$permutation_times)
  )
  testthat::expect_true(
    all(
      cluster_summary$permutation_p_value >=
        1 / (context$permutation_times + 1)
    )
  )
  testthat::expect_true(all(cluster_summary$permutation_p_value <= 1))
  expect_at_most(
    max(abs(observed$score - reference_scores[observed$barcode_w_prefix])),
    1e-12,
    "SCAVENGE trait-relevance-score delta"
  )
  testthat::expect_equal(
    unname(round(observed$p_val * context$permutation_times)),
    unname(reference_exceedance_counts)
  )
  reference_significant <- reference_exceedance_counts <=
    0.05 * context$permutation_times
  testthat::expect_identical(
    unname(observed$score_is_sig),
    unname(reference_significant)
  )
})

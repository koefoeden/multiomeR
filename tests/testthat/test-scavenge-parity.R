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

testthat::test_that("integration: SCAVENGE trait-relevance scores match the pinned reference", {
  load_algorithm_test_runtime()
  context <- make_SCAVENGE_test_context()
  testthat::expect_identical(
    context$seed_idx,
    get_SCAVENGE_seed_index(context$z_score, seed_percent = 0.05)
  )

  reference_scores <- get_SCAVENGE_reference_scores(context)
  observed <- get_SCAVENGE_TRS_tibble(
    chromVAR_z_score_record = list(
      GWAS_ID = "fixture",
      z_score_vec = context$z_score
    ),
    NN_graph = context$weighted_graph,
    restart_prob = context$restart_prob,
    seed_percent = 0.05,
    scale_percent = 0.1
  )
  expect_at_most(
    max(abs(observed$score - reference_scores[observed$barcode_w_prefix])),
    1e-12,
    "SCAVENGE trait-relevance-score delta"
  )
})

#!/usr/bin/env Rscript

if (!exists("weighted_nearest_neighbors_BPCells", mode = "function")) {
  source("R/bootstrap_helpers.R")
  load_project_runtime()
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

expect_at_least <- function(observed, threshold, label) {
  if (!is.finite(observed) || observed < threshold) {
    fail(
      label,
      ": observed ", format(observed, digits = 6),
      ", required at least ", format(threshold, digits = 6)
    )
  }
}

expect_at_most <- function(observed, threshold, label) {
  if (!is.finite(observed) || observed > threshold) {
    fail(
      label,
      ": observed ", format(observed, digits = 6),
      ", required at most ", format(threshold, digits = 6)
    )
  }
}

expect_reference_version <- function(package, expected) {
  observed <- as.character(utils::packageVersion(package))
  if (!identical(observed, expected)) {
    fail(package, " reference version changed from ", expected, " to ", observed, "; review validation tolerances")
  }
  observed
}

make_wnn_fixture <- function(n_cells, n_dimensions, n_clusters, ATAC_cluster_order, seed) {
  stopifnot(n_cells %% n_clusters == 0L, length(ATAC_cluster_order) == n_clusters)
  set.seed(seed)
  cells <- sprintf("cell%03d", seq_len(n_cells))
  cluster <- rep(seq_len(n_clusters), each = n_cells / n_clusters)
  centers <- matrix(stats::rnorm(n_clusters * n_dimensions, sd = 2), nrow = n_clusters)
  quality <- seq(0.15, 1.25, length.out = n_cells)

  RNA <- centers[cluster, , drop = FALSE] + matrix(
    stats::rnorm(n_cells * n_dimensions, sd = rep(quality, each = n_dimensions)),
    nrow = n_cells,
    byrow = TRUE
  )
  ATAC <- centers[ATAC_cluster_order[cluster], , drop = FALSE] + matrix(
    stats::rnorm(n_cells * n_dimensions, sd = rep(rev(quality), each = n_dimensions)),
    nrow = n_cells,
    byrow = TRUE
  )

  rownames(RNA) <- rownames(ATAC) <- cells
  colnames(RNA) <- paste0("RNA_", seq_len(n_dimensions))
  colnames(ATAC) <- paste0("ATAC_", seq_len(n_dimensions))
  list(RNA = RNA, ATAC = ATAC)
}

run_seurat_wnn_reference <- function(embeddings, k, candidate_k) {
  cells <- rownames(embeddings$RNA)
  counts <- Matrix::sparseMatrix(
    i = rep(1L, length(cells)),
    j = seq_along(cells),
    x = 1,
    dims = c(1L, length(cells)),
    dimnames = list("dummy", cells)
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["rna_fixture"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings$RNA,
    key = "RNA_",
    assay = "RNA"
  )
  object[["atac_fixture"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings$ATAC,
    key = "ATAC_",
    assay = "RNA"
  )

  set.seed(847)
  object <- Seurat::FindMultiModalNeighbors(
    object = object,
    reduction.list = list("rna_fixture", "atac_fixture"),
    dims.list = rep(list(seq_len(ncol(embeddings$RNA))), 2L),
    k.nn = k,
    knn.range = candidate_k,
    l2.norm = TRUE,
    modality.weight.name = c("RNA.weight", "ATAC.weight"),
    verbose = FALSE
  )

  list(
    weights = object[[]][, c("RNA.weight", "ATAC.weight"), drop = FALSE],
    nn_idx = object@neighbors$weighted.nn@nn.idx
  )
}

get_wnn_similarity_metrics <- function(embeddings, k, candidate_k) {
  observed <- weighted_nearest_neighbors_BPCells(
    embeddings_list = embeddings,
    k = k,
    candidate_k = candidate_k,
    threads = 1,
    ef = 1000
  )
  reference <- run_seurat_wnn_reference(embeddings, k = k, candidate_k = candidate_k)

  weight_correlations <- c(
    RNA = stats::cor(observed$modality_weights$RNA, reference$weights$RNA.weight, method = "spearman"),
    ATAC = stats::cor(observed$modality_weights$ATAC, reference$weights$ATAC.weight, method = "spearman")
  )
  neighbor_overlap <- vapply(seq_len(nrow(reference$nn_idx)), function(cell_idx) {
    length(intersect(reference$nn_idx[cell_idx, ], observed$nn_idx[cell_idx, ])) / k
  }, numeric(1))

  c(
    RNA_weight_spearman = weight_correlations[["RNA"]],
    ATAC_weight_spearman = weight_correlations[["ATAC"]],
    neighbor_overlap_mean = mean(neighbor_overlap),
    neighbor_overlap_q25 = as.numeric(stats::quantile(neighbor_overlap, 0.25))
  )
}

validate_wnn_fixture <- function(label, metrics, weight_threshold, mean_overlap_threshold, q25_overlap_threshold) {
  expect_at_least(
    min(metrics[c("RNA_weight_spearman", "ATAC_weight_spearman")]),
    weight_threshold,
    paste0("WNN ", label, " modality-weight Spearman correlation")
  )
  expect_at_least(
    metrics[["neighbor_overlap_mean"]],
    mean_overlap_threshold,
    paste0("WNN ", label, " mean neighbor-set overlap")
  )
  expect_at_least(
    metrics[["neighbor_overlap_q25"]],
    q25_overlap_threshold,
    paste0("WNN ", label, " first-quartile neighbor-set overlap")
  )

  cat(
    "WNN ", label, " validation ok",
    ": weight Spearman RNA=", format(metrics[["RNA_weight_spearman"]], digits = 4),
    ", ATAC=", format(metrics[["ATAC_weight_spearman"]], digits = 4),
    "; neighbor overlap mean=", format(metrics[["neighbor_overlap_mean"]], digits = 4),
    ", q25=", format(metrics[["neighbor_overlap_q25"]], digits = 4),
    "\n",
    sep = ""
  )
}

validate_wnn <- function() {
  seurat_version <- expect_reference_version("Seurat", "5.5.0")

  stress_metrics <- get_wnn_similarity_metrics(
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
    metrics = stress_metrics,
    weight_threshold = 0.80,
    mean_overlap_threshold = 0.85,
    q25_overlap_threshold = 0.75
  )

  production_metrics <- get_wnn_similarity_metrics(
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
    label = paste0("production-like against Seurat ", seurat_version),
    metrics = production_metrics,
    weight_threshold = 0.90,
    mean_overlap_threshold = 0.94,
    q25_overlap_threshold = 0.90
  )
}

# The compact reference functions below reproduce SCAVENGE 1.0.2 at commit
# 8ee8b173d965009a696b2a590d5b17b28b7cf851. Keeping the small fixture local
# makes CI deterministic without installing SCAVENGE's historical dependency
# stack or fetching network resources.
reference_SCAVENGE_random_walk <- function(graph, seed_cells, restart_prob = 0.05, stationary_cutoff = 1e-5) {
  graph <- methods::as(graph != 0, "dMatrix")
  graph <- t(t(graph) / Matrix::colSums(graph))
  restart <- setNames(numeric(nrow(graph)), rownames(graph))
  restart[seed_cells] <- 1
  restart <- restart / sum(restart)
  # randomWalk_sparse() transposes before calling its iterator, which
  # transposes once more. Preserve both operations in this compact reference.
  transition <- Matrix::t(Matrix::t(graph))
  score <- restart
  delta <- 1
  while (delta > stationary_cutoff) {
    next_score <- as.numeric(
      ((1 - restart_prob) * transition) %*% score
    ) + restart_prob * restart
    delta <- sum(abs(next_score - score))
    score <- next_score
  }
  names(score) <- rownames(graph)
  score
}

reference_SCAVENGE_seed_index <- function(z_score, seed_percent = 0.05) {
  p_value <- stats::pnorm(z_score, lower.tail = FALSE)
  if (sum(p_value <= 0.05) / length(p_value) > seed_percent) {
    rank(-z_score) <= floor(seed_percent * length(z_score))
  } else {
    p_value <= 0.05
  }
}

reference_SCAVENGE_scores <- function(propagation_score, z_score, scale_percent = 0.1) {
  ceiling <- stats::quantile(propagation_score, 0.95, names = FALSE)
  capped <- propagation_score
  capped[capped > ceiling] <- ceiling
  scaled <- (capped - min(capped)) / (max(capped) - min(capped))
  scale_idx <- rank(-z_score) <= floor(scale_percent * length(z_score))
  scaled * mean(z_score[scale_idx])
}

reference_SCAVENGE_sample_seed_indices <- function(
  graph,
  seed_idx,
  permutation_times
) {
  cell_table <- data.frame(cell = seq_len(nrow(graph)), degree = Matrix::colSums(graph))
  seed_table <- data.frame(
    seed = which(seed_idx),
    degree = Matrix::colSums(graph[, seed_idx, drop = FALSE])
  ) |>
    with(data.frame(table(degree)))
  cells_by_degree <- tapply(cell_table$cell, cell_table$degree, list)
  cells_by_degree <- cells_by_degree[names(cells_by_degree) %in% seed_table$degree]

  lapply(seq_len(permutation_times), function(permutation) {
    cells_by_degree |>
      mapply(FUN = sample, seed_table$Freq) |>
      unlist(use.names = FALSE) |>
      sort()
  })
}

make_SCAVENGE_fixture <- function() {
  n_cells <- 60L
  cells <- sprintf("cell%02d", seq_len(n_cells))
  adjacency <- matrix(0, n_cells, n_cells, dimnames = list(cells, cells))
  template_edges <- rbind(
    c(1L, 2L), c(1L, 3L), c(1L, 4L), c(2L, 3L),
    c(3L, 4L), c(4L, 5L), c(5L, 6L), c(5L, 7L),
    c(6L, 7L), c(7L, 8L), c(8L, 9L), c(8L, 10L),
    c(9L, 10L)
  )
  edge_index <- 0L
  add_edge <- function(from, to) {
    edge_index <<- edge_index + 1L
    weight <- c(0.25, 0.5, 1)[[(edge_index - 1L) %% 3L + 1L]]
    adjacency[from, to] <<- weight
    adjacency[to, from] <<- weight
  }
  for (block in 0:5) {
    block_start <- block * 10L
    for (edge in seq_len(nrow(template_edges))) {
      add_edge(
        block_start + template_edges[edge, 1],
        block_start + template_edges[edge, 2]
      )
    }
  }
  for (block in 0:5) {
    add_edge(block * 10L + 10L, ((block + 1L) %% 6L) * 10L + 10L)
  }
  z_score <- setNames(seq(0, 0.01, length.out = n_cells), cells)
  z_score[cells[1:3]] <- c(5, 4.5, 4)
  list(
    graph = Matrix::Matrix(adjacency, sparse = TRUE),
    z_score = z_score,
    block = rep(sprintf("block_%d", 1:6), each = 10L)
  )
}

validate_SCAVENGE <- function() {
  fixture <- make_SCAVENGE_fixture()
  weighted_graph <- fixture$graph
  graph <- get_SCAVENGE_adjacency_matrix(weighted_graph)
  z_score <- fixture$z_score
  seed_cells <- names(z_score)[seq_len(3L)]
  restart_prob <- 0.05

  observed_precise <- run_sparse_random_walk_with_restart(
    NN_graph = graph,
    seed_cells = seed_cells,
    restart_prob = restart_prob,
    stationary_cutoff = 1e-12
  )
  transition <- get_SCAVENGE_transition_matrix(graph)
  observed_precomputed <- run_sparse_random_walk_with_restart(
    NN_graph = graph,
    seed_cells = seed_cells,
    restart_prob = restart_prob,
    stationary_cutoff = 1e-12,
    transition_matrix = transition
  )
  transition_reuse_delta <- max(abs(observed_precise - observed_precomputed))
  expect_at_most(transition_reuse_delta, 0, "SCAVENGE precomputed-transition delta")

  restart <- setNames(numeric(nrow(graph)), rownames(graph))
  restart[seed_cells] <- 1 / length(seed_cells)
  closed_form <- solve(
    diag(nrow(graph)) - (1 - restart_prob) * as.matrix(transition),
    restart_prob * restart
  )
  closed_form_delta <- max(abs(observed_precise - closed_form))
  expect_at_most(closed_form_delta, 1e-10, "SCAVENGE random-walk closed-form delta")

  reference_precise <- reference_SCAVENGE_random_walk(
    weighted_graph,
    seed_cells,
    restart_prob,
    stationary_cutoff = 1e-12
  )
  reference_propagation_delta <- max(abs(observed_precise - reference_precise))
  expect_at_most(
    reference_propagation_delta,
    1e-12,
    "SCAVENGE pinned-reference propagation delta"
  )
  observed_score <- run_sparse_random_walk_with_restart(
    graph,
    seed_cells,
    restart_prob
  )
  reference_propagation <- reference_SCAVENGE_random_walk(
    weighted_graph,
    seed_cells,
    restart_prob
  )

  seed_idx <- reference_SCAVENGE_seed_index(z_score, seed_percent = 0.05)
  if (!identical(seed_idx, get_SCAVENGE_seed_index(z_score, seed_percent = 0.05))) {
    fail("SCAVENGE seed selection differs from the reference fixture")
  }
  reference_cells <- names(reference_propagation)[reference_propagation != 0]
  reference_scores <- reference_SCAVENGE_scores(
    propagation_score = reference_propagation[reference_cells],
    z_score = z_score,
    scale_percent = 0.1
  )

  permutation_times <- 199L
  set.seed(431)
  reference_samples <- reference_SCAVENGE_sample_seed_indices(
    graph,
    seed_idx,
    permutation_times
  )
  set.seed(431)
  native_samples <- sample_SCAVENGE_degree_matched_seed_indices(
    graph,
    seed_idx,
    permutation_times
  )
  if (!identical(reference_samples, native_samples)) {
    fail("SCAVENGE degree-matched seed samples differ from the pinned reference")
  }
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
  if (!all(vapply(singleton_samples, identical, logical(1), 2L))) {
    fail("SCAVENGE singleton degree strata do not retain their only candidate")
  }

  metadata_tibble <- tibble::tibble(
    barcode_w_prefix = rownames(graph),
    PCA_harmony_SNN_cluster_named = fixture$block,
    PCA_harmony_SNN_cluster_cell_type = fixture$block
  )
  cluster_index_record <- get_SCAVENGE_cluster_index_record(
    metadata_tibble,
    rownames(graph),
    "PCA_harmony_SNN"
  )
  native_statistics_1_core <- run_SCAVENGE_permutation_statistics(
    transition,
    native_samples,
    observed_score,
    cluster_index_record,
    cores = 1L,
    restart_prob = restart_prob,
    native_source_file = "src/scavenge_random_walk.cpp"
  )
  native_statistics_2_cores <- run_SCAVENGE_permutation_statistics(
    transition,
    native_samples,
    observed_score,
    cluster_index_record,
    cores = 2L,
    restart_prob = restart_prob,
    native_source_file = "src/scavenge_random_walk.cpp"
  )
  if (!identical(native_statistics_1_core, native_statistics_2_cores)) {
    fail("SCAVENGE permutation statistics differ across core counts")
  }

  reference_permutation_scores <- vapply(
    reference_samples,
    function(sampled_indices) {
      reference_SCAVENGE_random_walk(
        graph,
        rownames(graph)[sampled_indices],
        restart_prob
      )
    },
    numeric(nrow(graph))
  )
  reference_exceedance_counts <- rowSums(
    reference_permutation_scores > observed_score
  )
  if (!all(
    native_statistics_1_core$cell_exceedance_counts ==
      reference_exceedance_counts
  )) {
    fail("SCAVENGE streaming exceedance counts differ from the pinned reference")
  }

  set.seed(431)
  observed_result <- get_SCAVENGE_result_from_chromVAR_z_score_record(
    chromVAR_z_score_record = list(GWAS_ID = "fixture", z_score_vec = z_score),
    NN_graph = weighted_graph,
    metadata_tibble = metadata_tibble,
    graph_name = "PCA_harmony_SNN",
    cores = 2,
    permutation_times = permutation_times,
    restart_prob = restart_prob,
    seed_percent = 0.05,
    scale_percent = 0.1
  )
  observed <- observed_result$TRS_tibble
  cluster_summary <- observed_result$TRS_summary_tibble
  if (
    nrow(cluster_summary) != 12L ||
      any(cluster_summary$permutation_times != permutation_times) ||
      any(cluster_summary$permutation_p_value < 1 / (permutation_times + 1)) ||
      any(cluster_summary$permutation_p_value > 1)
  ) {
    fail("SCAVENGE cluster-permutation summary is malformed")
  }
  score_delta <- max(abs(
    observed$score - reference_scores[observed$barcode_w_prefix]
  ))
  expect_at_most(
    score_delta,
    1e-12,
    "SCAVENGE trait-relevance-score delta"
  )
  observed_exceedance_counts <- round(
    observed$p_val * permutation_times
  )
  if (!all(observed_exceedance_counts == reference_exceedance_counts)) {
    fail("SCAVENGE integrated empirical P-values do not preserve exceedance counts")
  }

  reference_significant <- reference_exceedance_counts <=
    0.05 * permutation_times
  if (!identical(unname(observed$score_is_sig), unname(reference_significant))) {
    fail("SCAVENGE significant-cell calls differ from the pinned reference")
  }

  cat(
    "SCAVENGE validation ok: reference 1.0.2@8ee8b173d965",
    "; closed-form max delta=", format(closed_form_delta, scientific = TRUE, digits = 3),
    "; precomputed-transition max delta=", format(transition_reuse_delta, scientific = TRUE, digits = 3),
    "; pinned-reference max delta=", format(reference_propagation_delta, scientific = TRUE, digits = 3),
    "; score max delta=", format(score_delta, scientific = TRUE, digits = 3),
    "; exact streamed exceedance counts",
    "; exact significant-cell calls",
    "\n",
    sep = ""
  )
}

validate_wnn()
validate_SCAVENGE()
cat("algorithm deviation validation ok\n")

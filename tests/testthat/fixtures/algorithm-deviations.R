load_algorithm_test_runtime <- function() {
  if (!exists("weighted_nearest_neighbors_BPCells", mode = "function")) {
    load_project_test_runtime()
  }
}

expect_at_least <- function(observed, threshold, label) {
  testthat::expect_true(
    is.finite(observed) && observed >= threshold,
    info = paste0(label, ": observed ", observed, ", required >= ", threshold)
  )
}

expect_at_most <- function(observed, threshold, label) {
  testthat::expect_true(
    is.finite(observed) && observed <= threshold,
    info = paste0(label, ": observed ", observed, ", required <= ", threshold)
  )
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

make_SCAVENGE_test_context <- function(permutation_times = 199L) {
  fixture <- make_SCAVENGE_fixture()
  weighted_graph <- fixture$graph
  graph <- get_SCAVENGE_adjacency_matrix(weighted_graph)
  z_score <- fixture$z_score
  seed_cells <- names(z_score)[seq_len(3L)]
  restart_prob <- 0.05
  transition <- get_SCAVENGE_transition_matrix(graph)
  observed_score <- run_sparse_random_walk_with_restart(
    graph,
    seed_cells,
    restart_prob
  )
  seed_idx <- reference_SCAVENGE_seed_index(z_score, seed_percent = 0.05)

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

  list(
    fixture = fixture,
    weighted_graph = weighted_graph,
    graph = graph,
    z_score = z_score,
    seed_cells = seed_cells,
    seed_idx = seed_idx,
    restart_prob = restart_prob,
    transition = transition,
    observed_score = observed_score,
    reference_samples = reference_samples,
    native_samples = native_samples,
    metadata_tibble = metadata_tibble,
    cluster_index_record = cluster_index_record,
    permutation_times = permutation_times
  )
}

get_SCAVENGE_reference_exceedance_counts <- function(context) {
  reference_permutation_scores <- vapply(
    context$reference_samples,
    function(sampled_indices) {
      reference_SCAVENGE_random_walk(
        context$graph,
        rownames(context$graph)[sampled_indices],
        context$restart_prob
      )
    },
    numeric(nrow(context$graph))
  )
  rowSums(reference_permutation_scores > context$observed_score)
}

get_SCAVENGE_reference_scores <- function(context) {
  reference_propagation <- reference_SCAVENGE_random_walk(
    context$weighted_graph,
    context$seed_cells,
    context$restart_prob
  )
  reference_cells <- names(reference_propagation)[reference_propagation != 0]
  reference_SCAVENGE_scores(
    propagation_score = reference_propagation[reference_cells],
    z_score = context$z_score,
    scale_percent = 0.1
  )
}

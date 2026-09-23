get_SCAVENGE_seed_index <- function(z_score_vec, seed_percent) {
  if (seed_percent <= 0 || seed_percent >= 1) {
    stop("seed_percent must be between 0 and 1.")
  }

  seed_idx <- stats::pnorm(z_score_vec, lower.tail = FALSE) <= 0.05
  max_seed_count <- max(1L, floor(seed_percent * length(z_score_vec)))
  if (sum(seed_idx) > max_seed_count) {
    seed_idx <- rank(-z_score_vec) <= max_seed_count
  }
  seed_idx
}

cap_values_by_quantile <- function(x, q_ceiling = 0.95) {
  pmin(x, stats::quantile(x, q_ceiling, names = FALSE, na.rm = TRUE))
}

min_max_scale_vec <- function(x) {
  x_range <- range(x, na.rm = TRUE)
  if (x_range[[1]] == x_range[[2]]) {
    return(x * 0)
  }
  (x - x_range[[1]]) / (x_range[[2]] - x_range[[1]])
}

#' Prepare a SCAVENGE adjacency matrix
#'
#' Convert a sparse weighted neighbor graph to the binary adjacency matrix
#' expected by the reference SCAVENGE implementation.
#'
#' @param NN_graph Sparse cell-by-cell neighbor graph.
#' @return Sparse binary adjacency matrix with the input dimnames.
#' @keywords internal

get_SCAVENGE_adjacency_matrix <- function(NN_graph) {
  adjacency_matrix <- Matrix::drop0(methods::as(
    methods::as(NN_graph, "dMatrix"),
    "generalMatrix"
  ))
  adjacency_matrix@x[] <- 1
  adjacency_matrix
}

#' Build a sparse SCAVENGE transition matrix
#'
#' Column-normalize a sparse nearest-neighbor graph for random walk with restart.
#'
#' @param NN_graph Sparse cell-by-cell adjacency matrix with identical row and
#'   column names and no degree-zero cells.
#' @return Column-normalized sparse transition matrix.
#' @keywords internal

get_SCAVENGE_transition_matrix <- function(NN_graph) {
  col_sums <- Matrix::colSums(NN_graph)
  if (any(col_sums == 0)) {
    stop("NN_graph contains degree-zero cells.")
  }
  transition_matrix <- methods::as(
    methods::as(NN_graph, "dMatrix"),
    "generalMatrix"
  )
  transition_matrix@x <- transition_matrix@x *
    rep.int(1 / col_sums, diff(transition_matrix@p))
  colnames(transition_matrix) <- NULL
  transition_matrix
}

#' Run sparse random walk with restart
#'
#' Propagate seed-cell signal over a sparse nearest-neighbor graph.
#'
#' @param NN_graph Sparse cell-by-cell adjacency matrix with identical row and
#'   column names and no degree-zero cells.
#' @param seed_cells Character vector of cell names used as restart seeds; all
#'   values must be row names of `NN_graph`.
#' @param restart_prob Probability of restarting at seed cells on each iteration;
#'   must be between 0 and 1.
#' @param stationary_cutoff L1-change threshold used to stop iterations once the
#'   score vector is stationary.
#' @param max_iter Maximum random-walk iterations before returning the latest score.
#' @param transition_matrix Precomputed column-normalized transition matrix.
#' @return Named numeric propagation score vector aligned to `NN_graph` row names.
#' @keywords internal

run_sparse_random_walk_with_restart <- function(
  NN_graph,
  seed_cells,
  restart_prob,
  stationary_cutoff = 1e-5,
  max_iter = 10000,
  transition_matrix = get_SCAVENGE_transition_matrix(NN_graph)
) {
  if (restart_prob <= 0 || restart_prob >= 1) {
    stop("restart_prob must be between 0 and 1.")
  }
  if (!all(seed_cells %in% rownames(NN_graph))) {
    stop("seed_cells contains cells not found in NN_graph.")
  }

  restart_vec <- numeric(nrow(NN_graph))
  names(restart_vec) <- rownames(NN_graph)
  restart_vec[seed_cells] <- 1 / length(seed_cells)

  score_vec <- restart_vec
  for (iteration in seq_len(max_iter)) {
    next_score_vec <- as.numeric((1 - restart_prob) * (transition_matrix %*% score_vec) + restart_prob * restart_vec)
    delta <- sum(abs(next_score_vec - score_vec))
    score_vec <- next_score_vec
    if (delta <= stationary_cutoff) {
      break
    }
  }

  names(score_vec) <- rownames(NN_graph)
  score_vec
}

drop_SCAVENGE_degree_zero_cells <- function(NN_graph) {
  repeat {
    keep_cells <- Matrix::colSums(NN_graph) != 0
    if (all(keep_cells)) {
      return(NN_graph)
    }
    NN_graph <- NN_graph[keep_cells, keep_cells, drop = FALSE]
  }
}

#' Propagate SCAVENGE seed signal from a chromVAR z-score record
#'
#' Restrict the binary graph to cells with finite z scores and neighbors, select
#' seed cells, and propagate their signal. Cells without propagated signal are
#' then removed together with any cells this leaves without neighbors.
#'
#' @param chromVAR_z_score_record List containing `GWAS_ID` and named cell-level
#'   `z_score_vec`.
#' @param NN_graph Sparse cell-by-cell neighbor graph; cells are intersected with
#'   the z-score vector before scoring.
#' @param restart_prob Restart probability for random-walk propagation.
#' @param seed_percent Fraction of highest z-score cells used as seed cells.
#' @return List with `GWAS_ID`, the filtered `z_score_vec` used for seed
#'   selection, and for the retained cells a binary adjacency `graph`, named
#'   logical `seed_idx`, and named `propagation_score_vec`. Scores are zero when
#'   no cell qualifies as a seed.
#' @keywords internal

get_SCAVENGE_propagation_record <- function(
  chromVAR_z_score_record,
  NN_graph,
  restart_prob,
  seed_percent
) {
  z_score_vec <- chromVAR_z_score_record$z_score_vec
  shared_cells <- intersect(names(z_score_vec), rownames(NN_graph))
  z_score_vec <- z_score_vec[shared_cells]
  graph <- get_SCAVENGE_adjacency_matrix(NN_graph[shared_cells, shared_cells])

  finite_z_score_cell_idx <- which(is.finite(z_score_vec) & z_score_vec <= 1000)
  graph <- drop_SCAVENGE_degree_zero_cells(graph[finite_z_score_cell_idx, finite_z_score_cell_idx])
  z_score_vec <- z_score_vec[rownames(graph)]
  seed_idx <- if (length(z_score_vec) == 0) logical() else get_SCAVENGE_seed_index(z_score_vec, seed_percent = seed_percent)
  record <- list(
    GWAS_ID = chromVAR_z_score_record$GWAS_ID,
    z_score_vec = z_score_vec,
    graph = graph,
    seed_idx = seed_idx,
    propagation_score_vec = z_score_vec * 0
  )
  if (!any(seed_idx)) {
    return(record)
  }

  propagation_score_vec <- run_sparse_random_walk_with_restart(
    NN_graph = graph,
    seed_cells = names(z_score_vec)[seed_idx],
    restart_prob = restart_prob
  )
  graph <- drop_SCAVENGE_degree_zero_cells(
    graph[propagation_score_vec != 0, propagation_score_vec != 0, drop = FALSE]
  )
  record$graph <- graph
  record$seed_idx <- seed_idx[rownames(graph)]
  record$propagation_score_vec <- propagation_score_vec[rownames(graph)]
  record
}

get_empty_TRS_tibble <- function() {
  tibble::tibble(
    barcode_w_prefix = character(),
    score = numeric(),
    GWAS_ID = character(),
    seed_idx = logical()
  )
}

#' Get SCAVENGE trait-relevance scores from a chromVAR z-score record
#'
#' Cap the propagation scores at their 0.95 quantile, min-max scale them, and
#' multiply by the mean z score of the top 1% of cells (at least one).
#'
#' @inheritParams get_SCAVENGE_propagation_record
#' @return Cell-level tibble with `barcode_w_prefix`, `score`, `GWAS_ID`, and
#'   `seed_idx`.
#' @keywords internal

get_SCAVENGE_TRS_tibble <- function(
  chromVAR_z_score_record,
  NN_graph,
  restart_prob,
  seed_percent
) {
  propagation_record <- get_SCAVENGE_propagation_record(
    chromVAR_z_score_record = chromVAR_z_score_record,
    NN_graph = NN_graph,
    restart_prob = restart_prob,
    seed_percent = seed_percent
  )
  propagation_score_vec <- propagation_record$propagation_score_vec
  if (length(propagation_score_vec) == 0) {
    return(get_empty_TRS_tibble())
  }

  z_score_vec <- sort(propagation_record$z_score_vec, decreasing = TRUE)
  scale_factor <- mean(z_score_vec[seq_len(max(1L, floor(0.01 * length(z_score_vec))))])
  tibble::tibble(
    barcode_w_prefix = names(propagation_score_vec),
    score = unname(min_max_scale_vec(cap_values_by_quantile(propagation_score_vec, q_ceiling = 0.95)) * scale_factor),
    GWAS_ID = propagation_record$GWAS_ID,
    seed_idx = unname(propagation_record$seed_idx)
  )
}

get_empty_TRS_summary_tibble <- function() {
  tibble::tibble(
    GWAS_ID = character(),
    grouping_col = character(),
    cluster = character(),
    n_cells = integer(),
    median_score = numeric(),
    mean_score = numeric(),
    q25_score = numeric(),
    q75_score = numeric(),
    min_score = numeric(),
    max_score = numeric()
  )
}

#' Summarize SCAVENGE TRS by groups
#'
#' Summarize SCAVENGE TRS distributions by graph-specific metadata groups.
#'
#' @param TRS_tibble Cell-level SCAVENGE TRS tibble keyed by `barcode_w_prefix`.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param graph_name Name of the neighbor graph or modality used in SCAVENGE group summaries.
#' @param group_by_cols Suffixes appended to `graph_name` to find metadata
#'   grouping columns, for example `cluster_named`.
#' @return Summary tibble with cell counts and score quantiles for each
#'   GWAS/group combination.
#' @keywords internal

summarize_SCAVENGE_TRS_by_groups <- function(
  TRS_tibble,
  metadata_tibble,
  graph_name,
  group_by_cols = c("cluster_named", "cluster_cell_type")
) {
  if (nrow(TRS_tibble) == 0) {
    return(get_empty_TRS_summary_tibble())
  }

  group_cols <- stringr::str_c(graph_name, "_", group_by_cols)
  cell_group_tibble <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::select(barcode_w_prefix, dplyr::all_of(group_cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(group_cols),
      names_to = "grouping_col",
      values_to = "cluster"
    ) |>
    dplyr::mutate(cluster = as.character(cluster)) |>
    dplyr::filter(!is.na(cluster), cluster != "")

  if (nrow(cell_group_tibble) == 0) {
    return(get_empty_TRS_summary_tibble())
  }

  TRS_tibble |>
    dplyr::inner_join(cell_group_tibble, by = "barcode_w_prefix", relationship = "many-to-many") |>
    dplyr::group_by(GWAS_ID, grouping_col, cluster) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_score = stats::median(score, na.rm = TRUE),
      mean_score = mean(score, na.rm = TRUE),
      q25_score = as.numeric(stats::quantile(score, 0.25, na.rm = TRUE, names = FALSE)),
      q75_score = as.numeric(stats::quantile(score, 0.75, na.rm = TRUE, names = FALSE)),
      min_score = min(score, na.rm = TRUE),
      max_score = max(score, na.rm = TRUE),
      .groups = "drop"
    )
}

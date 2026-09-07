get_SCAVENGE_seed_index <- function(z_score_vec, seed_percent = 0.05, p_value_cutoff = 0.05) {
  if (seed_percent <= 0 || seed_percent >= 1) {
    stop("seed_percent must be between 0 and 1.")
  }

  seed_idx <- stats::pnorm(z_score_vec, lower.tail = FALSE) <= p_value_cutoff
  max_seed_count <- max(1L, floor(seed_percent * length(z_score_vec)))
  if (sum(seed_idx) > max_seed_count) {
    seed_idx <- rank(-z_score_vec) <= max_seed_count
  }
  seed_idx
}

get_SCAVENGE_scale_factor <- function(z_score_vec, scale_percent = 0.01) {
  if (scale_percent <= 0 || scale_percent >= 1) {
    stop("scale_percent must be between 0 and 1.")
  }

  top_count <- max(1L, floor(scale_percent * length(z_score_vec)))
  mean(sort(z_score_vec, decreasing = TRUE)[seq_len(top_count)])
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

if (!exists("SCAVENGE_native_state_env", inherits = FALSE)) {
  SCAVENGE_native_state_env <- new.env(parent = emptyenv())
  SCAVENGE_native_state_env$dll_name <- NULL
}

load_SCAVENGE_native_library <- function(native_source_file) {
  if (!is.null(SCAVENGE_native_state_env$dll_name)) {
    return(SCAVENGE_native_state_env$dll_name)
  }

  build_dir <- tempfile("multiomeR_scavenge_")
  dir.create(build_dir)
  build_source_file <- file.path(build_dir, basename(native_source_file))
  if (!file.copy(native_source_file, build_source_file)) {
    stop(
      "Could not copy the SCAVENGE native source into its temporary build directory.",
      call. = FALSE
    )
  }
  shared_library_file <- file.path(
    build_dir,
    paste0("multiomeR_scavenge", .Platform$dynlib.ext)
  )
  compile_library <- function(env = character()) {
    system2(
      command = file.path(R.home("bin"), "R"),
      args = c("CMD", "SHLIB", "-o", shared_library_file, build_source_file),
      stdout = TRUE,
      stderr = TRUE,
      env = env
    )
  }
  build_output <- compile_library(
    c("PKG_CXXFLAGS=-fopenmp", "PKG_LIBS=-fopenmp")
  )
  build_status <- attr(build_output, "status")
  if (!is.null(build_status) && build_status != 0L) {
    openmp_output <- build_output
    unlink(c(
      shared_library_file,
      sub("[.]cpp$", ".o", build_source_file)
    ))
    build_output <- compile_library()
    build_status <- attr(build_output, "status")
    if (!is.null(build_status) && build_status != 0L) {
      stop(
        "Could not compile the SCAVENGE native random-walk helper with ",
        "OpenMP or its serial fallback:\n",
        paste(c(openmp_output, build_output), collapse = "\n"),
        call. = FALSE
      )
    }
  }

  loaded_library <- dyn.load(shared_library_file)
  SCAVENGE_native_state_env$dll_name <- loaded_library[["name"]]
  SCAVENGE_native_state_env$dll_name
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
  restart_prob = 0.05,
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

#' Sample SCAVENGE seeds within exact degree strata
#'
#' Reproduce the reference implementation's sequential base-R sampling while
#' storing only the sampled seed indices needed by the streamed native walks.
#'
#' @param NN_graph Binary sparse adjacency matrix without degree-zero cells.
#' @param seed_idx Named logical vector identifying observed seed cells.
#' @param permutation_times Number of degree-matched seed samples.
#' @return List of sorted one-based integer cell indices, one per permutation.
#' @keywords internal

sample_SCAVENGE_degree_matched_seed_indices <- function(
  NN_graph,
  seed_idx,
  permutation_times
) {
  seed_idx <- seed_idx[rownames(NN_graph)]
  degree_vec <- Matrix::colSums(NN_graph)
  seed_counts_by_degree <- table(degree_vec[seed_idx])
  cells_by_degree <- split(
    seq_along(degree_vec),
    degree_vec
  )[names(seed_counts_by_degree)]
  if (any(lengths(cells_by_degree) == 0L)) {
    stop("Could not resolve every SCAVENGE seed-degree group.")
  }
  sample_sizes <- as.integer(seed_counts_by_degree)
  lapply(seq_len(permutation_times), function(permutation) {
    sampled_indices <- Map(
      function(cell_indices, sample_size) {
        if (length(cell_indices) == 1L) {
          cell_indices
        } else {
          sample(cell_indices, sample_size)
        }
      },
      cells_by_degree,
      sample_sizes
    )
    sort(as.integer(unlist(sampled_indices, use.names = FALSE)))
  })
}

get_SCAVENGE_cluster_index_record <- function(
  metadata_tibble,
  cell_names,
  graph_name
) {
  grouping_cols <- paste0(
    graph_name,
    c("_cluster_named", "_cluster_cell_type")
  )
  aligned_metadata <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE)
  aligned_metadata <- aligned_metadata[
    match(cell_names, aligned_metadata$barcode_w_prefix),
    grouping_cols,
    drop = FALSE
  ]
  membership_tibble <- aligned_metadata |>
    dplyr::mutate(cell_index = seq_along(cell_names) - 1L) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(grouping_cols),
      names_to = "grouping_col",
      values_to = "cluster"
    ) |>
    dplyr::mutate(cluster = as.character(cluster)) |>
    dplyr::filter(!is.na(cluster), cluster != "")
  groups <- membership_tibble |>
    dplyr::distinct(grouping_col, cluster) |>
    dplyr::arrange(grouping_col, cluster) |>
    dplyr::mutate(group_id = dplyr::row_number())
  indexed_membership <- membership_tibble |>
    dplyr::inner_join(groups, by = c("grouping_col", "cluster")) |>
    dplyr::arrange(group_id, cell_index)
  cluster_cells <- split(
    indexed_membership$cell_index,
    factor(indexed_membership$group_id, levels = groups$group_id)
  )

  groups$n_cells <- lengths(cluster_cells)
  list(
    groups = groups,
    cluster_offsets = as.integer(c(0L, cumsum(lengths(cluster_cells)))),
    cluster_cell_indices = as.integer(unlist(cluster_cells, use.names = FALSE))
  )
}

run_SCAVENGE_permutation_statistics <- function(
  transition_matrix,
  sampled_cell_idx_list,
  observed_score_vec,
  cluster_index_record,
  cores,
  restart_prob,
  native_source_file
) {
  permutation_times <- length(sampled_cell_idx_list)
  core_count <- min(cores, permutation_times)
  chunk_count <- min(permutation_times, 4L * core_count)
  chunk_sizes <- rep(permutation_times %/% chunk_count, chunk_count)
  chunk_sizes[seq_len(permutation_times %% chunk_count)] <-
    chunk_sizes[seq_len(permutation_times %% chunk_count)] + 1L
  dll_name <- load_SCAVENGE_native_library(native_source_file)

  statistics <- .Call(
    "multiomeR_scavenge_permutation_statistics",
    transition_matrix@p,
    transition_matrix@i,
    transition_matrix@x,
    as.integer(c(0L, cumsum(lengths(sampled_cell_idx_list)))),
    as.integer(unlist(sampled_cell_idx_list, use.names = FALSE) - 1L),
    as.integer(c(0L, cumsum(chunk_sizes))),
    unname(observed_score_vec),
    cluster_index_record$cluster_offsets,
    cluster_index_record$cluster_cell_indices,
    as.double(restart_prob),
    as.double(1e-5),
    as.integer(10000L),
    as.integer(core_count),
    PACKAGE = dll_name
  )
  statistics$cell_exceedance_counts <-
    rowSums(statistics$cell_exceedance_counts)
  statistics
}

summarize_SCAVENGE_cluster_permutations <- function(
  observed_score_vec,
  cluster_index_record,
  cluster_statistics,
  permutation_times,
  GWAS_ID,
  graph_name
) {
  cluster_indices <- purrr::map(
    seq_len(nrow(cluster_index_record$groups)),
    \(i) seq.int(
      cluster_index_record$cluster_offsets[[i]] + 1L,
      cluster_index_record$cluster_offsets[[i + 1L]]
    )
  )
  cell_indices <- purrr::map(
    cluster_indices,
    \(idx) cluster_index_record$cluster_cell_indices[idx] + 1L
  )
  observed_median_raw_score <- purrr::map_dbl(
    cell_indices,
    \(idx) stats::median(observed_score_vec[idx])
  )

  cluster_index_record$groups |>
    dplyr::mutate(
      GWAS_ID = GWAS_ID,
      graph = graph_name,
      observed_median_raw_score = observed_median_raw_score,
      null_median_raw_score = apply(
        cluster_statistics,
        1,
        stats::median
      ),
      null_q95_raw_score = apply(
        cluster_statistics,
        1,
        stats::quantile,
        probs = 0.95,
        names = FALSE
      ),
      exceedance_count = rowSums(
        cluster_statistics >= observed_median_raw_score
      ),
      permutation_times = permutation_times,
      permutation_p_value = (exceedance_count + 1) / (permutation_times + 1)
    ) |>
    dplyr::mutate(
      permutation_p_adj_BH = stats::p.adjust(
        permutation_p_value,
        method = "BH"
      ),
      .by = grouping_col
    ) |>
    dplyr::select(
      GWAS_ID,
      graph,
      grouping_col,
      cluster,
      observed_median_raw_score,
      null_median_raw_score,
      null_q95_raw_score,
      exceedance_count,
      permutation_times,
      permutation_p_value,
      permutation_p_adj_BH
    )
}

get_empty_TRS_tibble <- function() {
  tibble::tibble(
    barcode_w_prefix = character(),
    score = numeric(),
    GWAS_ID = character(),
    seed_idx = logical(),
    p_val = numeric(),
    log10_p_val = numeric(),
    score_is_sig = logical()
  )
}

#' Get SCAVENGE results from a chromVAR z-score record
#'
#' Compute cell-level SCAVENGE TRS and empirical cluster-level significance from
#' the same degree-matched permutation random walks.
#'
#' @param chromVAR_z_score_record List containing `GWAS_ID` and named cell-level
#'   `z_score_vec`.
#' @param NN_graph Sparse cell-by-cell neighbor graph; cells are intersected with
#'   the z-score vector before scoring.
#' @param metadata_tibble Cell metadata containing the graph-specific named and
#'   cell-type cluster columns.
#' @param graph_name Graph-name prefix used to resolve cluster columns.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param max_z_score Upper z-score cap for finite-cell filtering before seed selection.
#' @param permutation_times Number of degree-matched permutations used for cell-
#'   and cluster-level empirical P-values.
#' @param restart_prob Restart probability for random-walk propagation.
#' @param seed_percent Fraction of highest z-score cells used as seed cells.
#' @param scale_percent Upper quantile used to derive the TRS scale factor from
#'   filtered z scores.
#' @param native_source_file Tracked C++ source for the shared-memory random-walk
#'   kernel.
#' @return List containing a cell-level `TRS_tibble` and a group-level
#'   `TRS_summary_tibble` with empirical cluster P-values.
#' @keywords internal

get_SCAVENGE_result_from_chromVAR_z_score_record <- function(
  chromVAR_z_score_record,
  NN_graph,
  metadata_tibble,
  graph_name,
  cores,
  max_z_score = 1000,
  permutation_times = 1000,
  restart_prob = 0.05,
  seed_percent = 0.05,
  scale_percent = 0.01,
  native_source_file = file.path(
    get_project_root(),
    "src",
    "scavenge_random_walk.cpp"
  )
) {
  z_score_vec <- chromVAR_z_score_record$z_score_vec
  GWAS_ID <- chromVAR_z_score_record$GWAS_ID
  empty_result <- function() {
    list(
      TRS_tibble = get_empty_TRS_tibble(),
      TRS_summary_tibble = get_empty_TRS_summary_tibble()
    )
  }
  if (permutation_times < 1) {
    stop("permutation_times must be at least 1.")
  }

  shared_cells <- intersect(names(z_score_vec), rownames(NN_graph))
  z_score_vec <- z_score_vec[shared_cells]
  NN_graph <- get_SCAVENGE_adjacency_matrix(
    NN_graph[shared_cells, shared_cells]
  )

  finite_z_score_cell_idx <- which(is.finite(z_score_vec) & z_score_vec <= max_z_score)
  z_score_vec_filtered <- z_score_vec[finite_z_score_cell_idx]
  z_score_filtered_graph <- NN_graph[finite_z_score_cell_idx, finite_z_score_cell_idx]

  deg0_filtered_graph <- drop_SCAVENGE_degree_zero_cells(z_score_filtered_graph)
  deg0_z_score_vec_filtered <- z_score_vec_filtered[rownames(deg0_filtered_graph)]
  if (length(deg0_z_score_vec_filtered) == 0) {
    return(empty_result())
  }

  is_seed_bool_vec <- get_SCAVENGE_seed_index(deg0_z_score_vec_filtered, seed_percent = seed_percent)
  if (!any(is_seed_bool_vec)) {
    TRS_tibble <- tibble::tibble(
      barcode_w_prefix = names(deg0_z_score_vec_filtered),
      score = 0,
      GWAS_ID = GWAS_ID,
      seed_idx = FALSE,
      p_val = 1,
      log10_p_val = 0,
      score_is_sig = FALSE
    )
    TRS_summary_tibble <- summarize_SCAVENGE_TRS_by_groups(
      TRS_tibble = TRS_tibble,
      metadata_tibble = metadata_tibble,
      graph_name = graph_name
    ) |>
      dplyr::mutate(
        graph = graph_name,
        observed_median_raw_score = 0,
        null_median_raw_score = NA_real_,
        null_q95_raw_score = NA_real_,
        exceedance_count = NA_integer_,
        permutation_times = 0L,
        permutation_p_value = 1,
        permutation_p_adj_BH = 1,
        .after = GWAS_ID
      )
    return(list(
      TRS_tibble = TRS_tibble,
      TRS_summary_tibble = TRS_summary_tibble
    ))
  }

  # Calculate network propagation score and filter out cells with no score
  net_prop_score_named_vec <- run_sparse_random_walk_with_restart(
    NN_graph = deg0_filtered_graph,
    seed_cells = rownames(deg0_filtered_graph)[is_seed_bool_vec],
    restart_prob = restart_prob
  )
  zero_net_prop_score <- net_prop_score_named_vec == 0
  net_prop_score_named_vec_filtered <- net_prop_score_named_vec[!zero_net_prop_score]
  kept_cells <- names(net_prop_score_named_vec_filtered)
  if (length(kept_cells) == 0) {
    return(empty_result())
  }

  triple_filtered_graph <- drop_SCAVENGE_degree_zero_cells(deg0_filtered_graph[kept_cells, kept_cells, drop = FALSE])
  kept_cells <- rownames(triple_filtered_graph)
  net_prop_score_named_vec_filtered <- net_prop_score_named_vec_filtered[kept_cells]
  if (length(kept_cells) == 0) {
    return(empty_result())
  }

  # Cap, scale, and multiply by scale factor
  scale_factor <- get_SCAVENGE_scale_factor(deg0_z_score_vec_filtered, scale_percent = scale_percent)

  cell_named_TRS_vec <- net_prop_score_named_vec_filtered %>%
    cap_values_by_quantile(q_ceiling = 0.95) %>%
    min_max_scale_vec() %>%
    magrittr::multiply_by(scale_factor)

  seed_idx <- is_seed_bool_vec[kept_cells]
  cluster_index_record <- get_SCAVENGE_cluster_index_record(
    metadata_tibble = metadata_tibble,
    cell_names = kept_cells,
    graph_name = graph_name
  )
  sampled_cell_idx_list <- sample_SCAVENGE_degree_matched_seed_indices(
    NN_graph = triple_filtered_graph,
    seed_idx = seed_idx,
    permutation_times = permutation_times
  )
  permutation_statistics <- run_SCAVENGE_permutation_statistics(
    transition_matrix = get_SCAVENGE_transition_matrix(triple_filtered_graph),
    sampled_cell_idx_list = sampled_cell_idx_list,
    observed_score_vec = net_prop_score_named_vec_filtered,
    cluster_index_record = cluster_index_record,
    cores = cores,
    restart_prob = restart_prob,
    native_source_file = native_source_file
  )

  TRS_tibble <- cell_named_TRS_vec |>
    tibble::enframe(name = "barcode_w_prefix", value = "score") |>
    dplyr::mutate(
      GWAS_ID = GWAS_ID,
      seed_idx = unname(seed_idx),
      p_val = permutation_statistics$cell_exceedance_counts /
        permutation_times,
      log10_p_val = -log10(p_val),
      score_is_sig = p_val <= 0.05
    )
  cluster_permutation_tibble <- summarize_SCAVENGE_cluster_permutations(
    observed_score_vec = net_prop_score_named_vec_filtered,
    cluster_index_record = cluster_index_record,
    cluster_statistics = permutation_statistics$cluster_statistics,
    permutation_times = permutation_times,
    GWAS_ID = GWAS_ID,
    graph_name = graph_name
  )
  TRS_summary_tibble <- summarize_SCAVENGE_TRS_by_groups(
    TRS_tibble = TRS_tibble,
    metadata_tibble = metadata_tibble,
    graph_name = graph_name
  ) |>
    dplyr::left_join(
      cluster_permutation_tibble,
      by = c("GWAS_ID", "grouping_col", "cluster")
    )

  list(
    TRS_tibble = TRS_tibble,
    TRS_summary_tibble = TRS_summary_tibble
  )
}

get_empty_TRS_summary_tibble <- function() {
  tibble::tibble(
    GWAS_ID = character(),
    grouping_col = character(),
    cluster = character(),
    n_cells = integer(),
    n_sig = integer(),
    prop_sig = numeric(),
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
#' @return Summary tibble with cell counts, significant-cell fractions, and score
#'   quantiles for each GWAS/group combination.
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

  full_group_cols <- stringr::str_c(graph_name, "_", group_by_cols)
  available_group_cols <- intersect(full_group_cols, colnames(metadata_tibble))
  if (length(available_group_cols) == 0) {
    return(get_empty_TRS_summary_tibble())
  }

  cell_group_tibble <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::select(barcode_w_prefix, dplyr::all_of(available_group_cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(available_group_cols),
      names_to = "grouping_col",
      values_to = "cluster"
    ) |>
    dplyr::mutate(cluster = as.character(cluster)) |>
    dplyr::filter(!is.na(cluster), cluster != "")

  if (nrow(cell_group_tibble) == 0) {
    return(get_empty_TRS_summary_tibble())
  }

  TRS_tibble |>
    dplyr::inner_join(cell_group_tibble, by = "barcode_w_prefix") |>
    dplyr::group_by(GWAS_ID, grouping_col, cluster) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_sig = sum(score_is_sig, na.rm = TRUE),
      prop_sig = n_sig / n_cells,
      median_score = stats::median(score, na.rm = TRUE),
      mean_score = mean(score, na.rm = TRUE),
      q25_score = as.numeric(stats::quantile(score, 0.25, na.rm = TRUE, names = FALSE)),
      q75_score = as.numeric(stats::quantile(score, 0.75, na.rm = TRUE, names = FALSE)),
      min_score = min(score, na.rm = TRUE),
      max_score = max(score, na.rm = TRUE),
      .groups = "drop"
    )
}

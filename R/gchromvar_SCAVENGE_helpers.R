#' Filter credible set variants
#'
#' Filter and optionally reweight credible-set variants before peak overlap.
#'
#' @param credible_set_GRanges GRanges object containing credible set GRanges coordinates and metadata.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param ... Additional arguments forwarded to the variant weighting function.
#' @return A GRanges object containing retained variants, with any weighting
#'   function side effects applied to metadata columns.
#' @keywords internal

filter_credible_set_variants <- function(credible_set_GRanges, posterior_probability_cutoff = NULL, posterior_probability_weighting_function = NULL, ...) {
  if (!is.null(posterior_probability_cutoff)) {
    credible_set_GRanges <- S4Vectors::subset(credible_set_GRanges, posteriorProbability > posterior_probability_cutoff)
  }

  if (!is.null(posterior_probability_weighting_function)) {
    credible_set_GRanges <- posterior_probability_weighting_function(credible_set_GRanges, ...)
  }

  credible_set_GRanges
}

#' Sum variant weights in peaks
#'
#' Sum overlapping variant weights for each peak range.
#'
#' @param variant_GRanges GRanges object containing variant GRanges coordinates and metadata.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param weight_col Column in `variant_tibble` to draw as vertical PIP/weight values.
#' @return Named numeric vector with one value per peak; peaks without variants
#'   receive zero.
#' @keywords internal

sum_variant_weights_in_peaks <- function(variant_GRanges, peak_ranges, weight_col = "posteriorProbability") {
  if (!weight_col %in% names(GenomicRanges::mcols(variant_GRanges))) {
    stop("variant_GRanges is missing weight column: ", weight_col)
  }

  peak_variant_overlaps_hits <- GenomicRanges::findOverlaps(
    query = peak_ranges,
    subject = variant_GRanges
  )

  out <- numeric(length(peak_ranges))
  names(out) <- get_peak_names_from_GRanges(peak_ranges)
  if (length(peak_variant_overlaps_hits) == 0) {
    return(out)
  }

  peak_idx <- S4Vectors::queryHits(peak_variant_overlaps_hits)
  variant_idx <- S4Vectors::subjectHits(peak_variant_overlaps_hits)
  variant_weights <- GenomicRanges::mcols(variant_GRanges)[[weight_col]][variant_idx]
  variant_weights[is.na(variant_weights)] <- 0

  out <- Matrix::sparseMatrix(
    i = peak_idx,
    j = rep.int(1L, length(peak_idx)),
    x = variant_weights,
    dims = c(length(peak_ranges), 1)
  ) |>
    Matrix::rowSums()
  names(out) <- get_peak_names_from_GRanges(peak_ranges)
  out
}

#' Get summed posterior probabilities per peak
#'
#' Convert credible-set variant PIPs into peak-level GWAS weights.
#'
#' @param credible_set_GRanges GRanges object containing credible set GRanges coordinates and metadata.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param GWAS_ID Configured GWAS label used in target names, plots, and Open Targets joins.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param weight_transform Weight post-processing mode. `cap_1` caps summed peak
#'   weights at 1; `sum` leaves summed weights unchanged.
#' @param ... Additional arguments forwarded to `filter_credible_set_variants()`.
#' @return Named numeric vector of peak weights, with names matching peak ranges.
#'   Errors if no credible-set variants overlap peaks for the GWAS.
#' @keywords internal

get_summed_posterior_probabilities_per_peak <- function(
  credible_set_GRanges,
  peak_ranges,
  GWAS_ID,
  posterior_probability_cutoff = NULL,
  posterior_probability_weighting_function = NULL,
  weight_transform = "cap_1",
  ...
) {
  credible_set_GRanges <- credible_set_GRanges |>
    filter_credible_set_variants(
      posterior_probability_cutoff = posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function,
      ...
    )

  posterior_probability_sums_per_peak <- sum_variant_weights_in_peaks(
    variant_GRanges = credible_set_GRanges,
    peak_ranges = peak_ranges,
    weight_col = "posteriorProbability"
  )

  if (sum(posterior_probability_sums_per_peak > 0) == 0) {
    stop(sprintf("No peaks overlap with credible-set variants for GWAS trait: %s - skipped.", GWAS_ID))
  }

  if (identical(weight_transform, "cap_1")) {
    n_capped <- sum(posterior_probability_sums_per_peak > 1)
    if (n_capped > 0) {
      message(sprintf(
        "Capped %s GWAS_chromVAR peak weight(s) above 1 for %s; max uncapped weight was %.3f.",
        n_capped,
        GWAS_ID,
        max(posterior_probability_sums_per_peak)
      ))
    }
    posterior_probability_sums_per_peak <- pmin(posterior_probability_sums_per_peak, 1)
  } else if (!identical(weight_transform, "sum")) {
    stop("Unsupported weight_transform: ", weight_transform)
  }

  posterior_probability_sums_per_peak
}

#' Get GWAS chromVAR peak weight record
#'
#' Build one peak-weight record for GWAS chromVAR scoring.
#'
#' @param GWAS_input_record Single GWAS branch record containing the study, finemapping method, and weighting mode.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param posterior_probability_cutoff Minimum posterior probability/PIP retained before assigning variants to peaks.
#' @param posterior_probability_weighting_function Function applied to credible-set variants before peak weights are summed; receives the variant tibble.
#' @param weight_transform Optional function or scalar transform applied to variant weights before aggregation.
#' @param ... Additional arguments forwarded to peak-weight construction helpers.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_GWAS_chromVAR_peak_weight_record <- function(
  GWAS_input_record,
  peak_ranges,
  posterior_probability_cutoff = NULL,
  posterior_probability_weighting_function = NULL,
  weight_transform = "cap_1",
  ...
) {
  GWAS_ID <- GWAS_input_record$GWAS_ID
  list(
    GWAS_ID = GWAS_ID,
    peak_weights_vec = get_summed_posterior_probabilities_per_peak(
      credible_set_GRanges = GWAS_input_record$credible_set_GRanges,
      peak_ranges = peak_ranges,
      GWAS_ID = GWAS_ID,
      posterior_probability_cutoff = posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function,
      weight_transform = weight_transform,
      ...
    )
  )
}

get_GWAS_chromVAR_peak_weight_summary_tibble <- function(peak_weight_records) {
  total_peaks <- length(peak_weight_records[[1]]$peak_weights_vec)

  peak_weight_records |>
    purrr::map_dfr(\(peak_weight_record) {
      tibble::tibble(
        GWAS_ID = peak_weight_record$GWAS_ID,
        total_capped_posteriorProbability = sum(peak_weight_record$peak_weights_vec),
        frac_overlapped_peaks = sum(peak_weight_record$peak_weights_vec > 0) / total_peaks
      )
    })
}

get_GWAS_chromVAR_peak_weight_matrix <- function(peak_weight_records, RSE_ATAC) {
  peak_weights_matrix <- peak_weight_records |>
    purrr::map("peak_weights_vec") |>
    do.call(what = cbind)
  colnames(peak_weights_matrix) <- purrr::map_chr(peak_weight_records, "GWAS_ID")

  peak_weights_matrix |>
    align_peak_weights_to_RSE(RSE_ATAC = RSE_ATAC) |>
    Matrix::Matrix(sparse = TRUE)
}

#' Plot GWAS chromVAR peak weights summary
#'
#' Plot how many peaks receive nonzero GWAS chromVAR weights per trait.
#'
#' @param peak_weight_records List of GWAS peak-weight records, each containing
#'   `GWAS_ID` and a named `peak_weights_vec`.
#' @param overlap_threshold Number of nonzero peak overlaps used as the dashed
#'   threshold in the overlap-fraction facet.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_chromVAR_peak_weights_summary <- function(peak_weight_records, overlap_threshold = 100) {
  total_peaks <- length(peak_weight_records[[1]]$peak_weights_vec)
  overlap_threshold_fraction <- overlap_threshold / total_peaks

  plot_tibble <- peak_weight_records |>
    get_GWAS_chromVAR_peak_weight_summary_tibble() %>%
    tidyr::pivot_longer(!dplyr::matches("GWAS_ID|_sheet"))

  vline_tibble <- tibble::tibble(
    name = "frac_overlapped_peaks",
    xintercept = overlap_threshold_fraction
  )

  plot_tibble %>%
    ggplot2::ggplot(ggplot2::aes(y = GWAS_ID, x = value)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      data = vline_tibble,
      ggplot2::aes(xintercept = xintercept),
      inherit.aes = FALSE,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~name, scales = "free_x") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = stringr::str_glue("Peak x GWAS posterior-probability overlap summary, from {total_peaks} total peaks"),
      caption = stringr::str_glue(
        "Dashed line in the frac_overlapped_peaks facet marks {overlap_threshold} overlapped peaks ({scales::percent(overlap_threshold_fraction)} of all peaks)."
      )
    )
}

align_peak_weights_to_RSE <- function(peak_weights_matrix, RSE_ATAC) {
  peak_names <- get_peak_names_from_GRanges(SummarizedExperiment::rowRanges(RSE_ATAC))
  if (is.null(rownames(peak_weights_matrix))) {
    if (nrow(peak_weights_matrix) != length(peak_names)) {
      stop("Unnamed peak_weights_matrix must have the same number of rows as RSE_ATAC peaks.")
    }
    rownames(peak_weights_matrix) <- peak_names
  }

  missing_peaks <- setdiff(peak_names, rownames(peak_weights_matrix))
  if (length(missing_peaks) > 0) {
    stop("Missing peak weight rows for ", length(missing_peaks), " ATAC peak(s).")
  }

  peak_weights_matrix[peak_names, , drop = FALSE]
}

#' Get GWAS chromVAR z-score chunk record
#'
#' Compute one GWAS chromVAR z-score vector for one cell chunk.
#'
#' @param peak_weight_record One GWAS record containing `GWAS_ID` and named
#'   peak weights.
#' @param chunk_context_record Cell-chunk record containing counts, background,
#'   chunk ID, and cell names.
#' @param RSE_ATAC RangedSummarizedExperiment for ATAC peaks, with row ranges aligned to peak-level matrices.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_GWAS_chromVAR_z_score_chunk_record <- function(peak_weight_record, chunk_context_record, RSE_ATAC) {
  peak_weights_matrix <- matrix(peak_weight_record$peak_weights_vec, ncol = 1)
  rownames(peak_weights_matrix) <- names(peak_weight_record$peak_weights_vec)
  colnames(peak_weights_matrix) <- peak_weight_record$GWAS_ID
  peak_weights_matrix <- align_peak_weights_to_RSE(
    peak_weights_matrix = peak_weights_matrix,
    RSE_ATAC = RSE_ATAC
  )
  peak_weights_matrix <- Matrix::Matrix(peak_weights_matrix, sparse = TRUE)

  chunk_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = chunk_context_record$counts),
    rowRanges = SummarizedExperiment::rowRanges(RSE_ATAC)
  )
  SummarizedExperiment::rowData(chunk_obj) <- SummarizedExperiment::rowData(RSE_ATAC)

  z_score_matrix <- betterChromVAR::computeDeviationsAnalytic(
    object = chunk_obj,
    background = chunk_context_record$background,
    annotations = peak_weights_matrix,
    verbose = FALSE,
    retSE = FALSE,
    compute = "z"
  )$z

  z_score_vec <- z_score_matrix[peak_weight_record$GWAS_ID, ]
  names(z_score_vec) <- chunk_context_record$cell_names
  tibble::tibble(
    GWAS_ID = peak_weight_record$GWAS_ID,
    chunk_id = chunk_context_record$chunk_id,
    z_score_vec = list(z_score_vec)
  )
}

combine_GWAS_chromVAR_z_score_chunk_records <- function(chromVAR_z_score_chunk_records) {
  chromVAR_z_score_chunk_records <- dplyr::arrange(chromVAR_z_score_chunk_records, chunk_id)
  list(
    GWAS_ID = chromVAR_z_score_chunk_records$GWAS_ID[[1]],
    z_score_vec = unlist(chromVAR_z_score_chunk_records$z_score_vec, use.names = TRUE)
  )
}

get_SCAVENGE_seed_index <- function(z_score_vec, seed_percent = 0.05, p_value_cutoff = 0.05) {
  if (seed_percent <= 0 || seed_percent >= 1) {
    stop("seed_percent must be between 0 and 1.")
  }

  seed_idx <- stats::pnorm(z_score_vec, lower.tail = FALSE) <= p_value_cutoff
  max_seed_count <- max(1L, floor(seed_percent * length(z_score_vec)))
  if (sum(seed_idx) > max_seed_count) {
    seed_idx <- rank(-z_score_vec, ties.method = "first") <= max_seed_count
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

scale_GWAS_heatmap_scores <- function(heatmap_data, group_cols = character(), score_col = "median_score") {
  if (nrow(heatmap_data) == 0) {
    return(heatmap_data)
  }

  heatmap_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("GWAS_ID", group_cols)))) |>
    dplyr::mutate(!!score_col := min_max_scale_vec(.data[[score_col]])) |>
    dplyr::ungroup()
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
#' @return Named numeric propagation score vector aligned to `NN_graph` row names.
#' @keywords internal

run_sparse_random_walk_with_restart <- function(NN_graph, seed_cells, restart_prob = 0.05, stationary_cutoff = 1e-5, max_iter = 10000) {
  if (restart_prob <= 0 || restart_prob >= 1) {
    stop("restart_prob must be between 0 and 1.")
  }
  if (!all(seed_cells %in% rownames(NN_graph))) {
    stop("seed_cells contains cells not found in NN_graph.")
  }

  col_sums <- Matrix::colSums(NN_graph)
  if (any(col_sums == 0)) {
    stop("NN_graph contains degree-zero cells.")
  }

  transition_matrix <- NN_graph %*% Matrix::Diagonal(x = 1 / col_sums)
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

#' Get degree matched permutation p values
#'
#' Estimate empirical SCAVENGE p values from degree-matched seed permutations.
#'
#' @param NN_graph Sparse cell adjacency matrix used for random-walk propagation.
#' @param seed_idx Logical vector, named by cell, identifying observed seed cells.
#' @param observed_score_vec Named observed propagation scores to compare against
#'   permutations.
#' @param permutation_times Number of degree-matched random seed sets to sample.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param restart_prob Restart probability passed to the random-walk helper.
#' @return Tibble with one row per cell and empirical permutation p value.
#' @keywords internal

get_degree_matched_permutation_p_values <- function(NN_graph, seed_idx, observed_score_vec, permutation_times = 1000, cores = 1, restart_prob = 0.05) {
  if (permutation_times < 1) {
    stop("permutation_times must be at least 1.")
  }

  observed_score_vec <- observed_score_vec[rownames(NN_graph)]
  seed_idx <- seed_idx[rownames(NN_graph)]
  degree_vec <- Matrix::colSums(NN_graph)
  cells_by_degree <- split(seq_along(degree_vec), degree_vec)
  seed_counts_by_degree <- table(degree_vec[seed_idx])

  run_permutations <- function(n_permutations) {
    exceedance_counts <- integer(length(observed_score_vec))
    for (i in seq_len(n_permutations)) {
      sampled_cell_idx <- names(seed_counts_by_degree) |>
        purrr::map(\(degree) {
          degree_cell_idx <- cells_by_degree[[degree]]
          degree_cell_idx[sample.int(length(degree_cell_idx), seed_counts_by_degree[[degree]])]
        }) |>
        unlist(use.names = FALSE) |>
        sort()

      permutation_score_vec <- run_sparse_random_walk_with_restart(
        NN_graph = NN_graph,
        seed_cells = rownames(NN_graph)[sampled_cell_idx],
        restart_prob = restart_prob
      )

      exceedance_counts <- exceedance_counts + as.integer(permutation_score_vec > observed_score_vec)
    }
    exceedance_counts
  }

  core_count <- min(cores, permutation_times)
  chunk_sizes <- rep(permutation_times %/% core_count, core_count)
  chunk_sizes[seq_len(permutation_times %% core_count)] <- chunk_sizes[seq_len(permutation_times %% core_count)] + 1L
  chunk_sizes <- chunk_sizes[chunk_sizes > 0]

  chunk_results <- if (core_count > 1) {
    parallel::mclapply(chunk_sizes, run_permutations, mc.cores = core_count)
  } else {
    lapply(chunk_sizes, run_permutations)
  }

  p_val <- (Reduce(`+`, chunk_results) + 1) / (permutation_times + 1)
  tibble::tibble(
    rowname = names(observed_score_vec),
    seed_idx = unname(seed_idx),
    p_val = p_val
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

#' Get SCAVENGE TRS from a chromVAR z-score record
#'
#' Convert one GWAS chromVAR z-score record into SCAVENGE TRS scores.
#'
#' @param chromVAR_z_score_record List containing `GWAS_ID` and named cell-level
#'   `z_score_vec`.
#' @param NN_graph Sparse cell-by-cell neighbor graph; cells are intersected with
#'   the z-score vector before scoring.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param max_z_score Upper z-score cap for finite-cell filtering before seed selection.
#' @param permutation_times Number of degree-matched permutations used for cell
#'   p values.
#' @param restart_prob Restart probability for random-walk propagation.
#' @param seed_percent Fraction of highest z-score cells used as seed cells.
#' @param scale_percent Upper quantile used to derive the TRS scale factor from
#'   filtered z scores.
#' @return Tibble with cell barcode, TRS score, GWAS ID, seed flag, permutation
#'   p value, `-log10(p)`, and significance flag.
#' @keywords internal

get_SCAVENGE_TRS_from_chromVAR_z_score_record <- function(
  chromVAR_z_score_record,
  NN_graph,
  cores,
  max_z_score = 1000,
  permutation_times = 1000,
  restart_prob = 0.05,
  seed_percent = 0.05,
  scale_percent = 0.01
) {
  z_score_vec <- chromVAR_z_score_record$z_score_vec
  GWAS_ID <- chromVAR_z_score_record$GWAS_ID

  shared_cells <- intersect(names(z_score_vec), rownames(NN_graph))
  z_score_vec <- z_score_vec[shared_cells]
  NN_graph <- NN_graph[shared_cells, shared_cells]

  finite_z_score_cell_idx <- which(is.finite(z_score_vec) & z_score_vec <= max_z_score)
  z_score_vec_filtered <- z_score_vec[finite_z_score_cell_idx]
  z_score_filtered_graph <- NN_graph[finite_z_score_cell_idx, finite_z_score_cell_idx]

  deg0_filtered_graph <- drop_SCAVENGE_degree_zero_cells(z_score_filtered_graph)
  deg0_z_score_vec_filtered <- z_score_vec_filtered[rownames(deg0_filtered_graph)]
  if (length(deg0_z_score_vec_filtered) == 0) {
    return(get_empty_TRS_tibble())
  }

  is_seed_bool_vec <- get_SCAVENGE_seed_index(deg0_z_score_vec_filtered, seed_percent = seed_percent)
  if (!any(is_seed_bool_vec)) {
    return(tibble::tibble(
      barcode_w_prefix = names(deg0_z_score_vec_filtered),
      score = 0,
      GWAS_ID = GWAS_ID,
      seed_idx = FALSE,
      p_val = 1,
      log10_p_val = 0,
      score_is_sig = FALSE
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
    return(get_empty_TRS_tibble())
  }

  triple_filtered_graph <- drop_SCAVENGE_degree_zero_cells(deg0_filtered_graph[kept_cells, kept_cells, drop = FALSE])
  kept_cells <- rownames(triple_filtered_graph)
  net_prop_score_named_vec_filtered <- net_prop_score_named_vec_filtered[kept_cells]
  if (length(kept_cells) == 0) {
    return(get_empty_TRS_tibble())
  }

  # Cap, scale, and multiply by scale factor
  scale_factor <- get_SCAVENGE_scale_factor(deg0_z_score_vec_filtered, scale_percent = scale_percent)

  cell_named_TRS_vec <- net_prop_score_named_vec_filtered %>%
    cap_values_by_quantile(q_ceiling = 0.95) %>%
    min_max_scale_vec() %>%
    magrittr::multiply_by(scale_factor)

  # Permutation test
  perm_test_mod <- get_degree_matched_permutation_p_values(
    seed_idx = is_seed_bool_vec[kept_cells],
    NN_graph = triple_filtered_graph,
    observed_score_vec = net_prop_score_named_vec_filtered,
    permutation_times = permutation_times,
    cores = cores,
    restart_prob = restart_prob
  )

  out_tibble <- cell_named_TRS_vec %>%
    tibble::enframe(name = "rowname", value = "score") %>%
    dplyr::mutate(GWAS_ID = GWAS_ID) %>%
    dplyr::left_join(perm_test_mod, by = "rowname") %>%
    dplyr::mutate(
      log10_p_val = -log10(p_val),
      score_is_sig = p_val < 0.05
    ) %>%
    dplyr::rename(barcode_w_prefix = rowname)

  return(out_tibble)
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

#' Get SCAVENGE TRS UMAP plots
#'
#' Build UMAP overlay plots from cell-level SCAVENGE TRS scores.
#'
#' @param TRS_tibble Cell-level TRS tibble with `barcode_w_prefix`, `GWAS_ID`,
#'   `score`, and significance columns.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param umap_cols Two metadata columns used as UMAP x/y coordinates.
#' @param label_col Optional metadata column used to label group centroids on
#'   the UMAP.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

get_SCAVENGE_TRS_UMAP_plots <- function(TRS_tibble, metadata_tibble, umap_cols, label_col = NULL) {
  if (nrow(TRS_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  umap_cols <- unlist(as.list(umap_cols), use.names = FALSE)
  GWAS_ID <- unique(TRS_tibble$GWAS_ID)
  score_col <- stringr::str_c("score_", GWAS_ID[[1]])
  SCAVENGE_metadata_tibble <- metadata_tibble |>
    dplyr::left_join(dplyr::select(TRS_tibble, barcode_w_prefix, score), by = "barcode_w_prefix") |>
    dplyr::rename(!!score_col := score)

  plot <- plot_UMAP_from_metadata(
    metadata_tibble = SCAVENGE_metadata_tibble,
    variable = score_col,
    umap_cols = umap_cols,
    legend_continuous = "value",
    quantile_range = c(0, 1)
  )

  if (inherits(plot, "empty_plot_list")) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot <- plot +
    ggplot2::scale_color_viridis_c(name = "TRS score") +
    ggplot2::labs(subtitle = NULL)

  if (!is.null(label_col)) {
    if (!is.character(label_col) || length(label_col) != 1 || is.na(label_col) || !nzchar(label_col)) {
      stop("`label_col` must be NULL or a non-empty length-1 character vector.", call. = FALSE)
    }
    if (!label_col %in% colnames(SCAVENGE_metadata_tibble)) {
      stop("Required label column not found: ", label_col, call. = FALSE)
    }

    label_tibble <- SCAVENGE_metadata_tibble |>
      dplyr::filter(
        !is.na(.data[[label_col]]),
        is.finite(.data[[umap_cols[[1]]]]),
        is.finite(.data[[umap_cols[[2]]]])
      ) |>
      dplyr::summarise(
        label = dplyr::first(as.character(.data[[label_col]])),
        UMAP_1 = stats::median(.data[[umap_cols[[1]]]]),
        UMAP_2 = stats::median(.data[[umap_cols[[2]]]]),
        .by = dplyr::all_of(label_col)
      ) |>
      dplyr::select(label, UMAP_1, UMAP_2)

    if (nrow(label_tibble) > 0) {
      plot <- plot +
        ggrepel::geom_label_repel(
          data = label_tibble,
          ggplot2::aes(x = UMAP_1, y = UMAP_2, label = label),
          inherit.aes = FALSE,
          size = 2.4,
          label.size = 0.15,
          label.padding = grid::unit(0.08, "lines"),
          min.segment.length = 0,
          max.overlaps = Inf,
          seed = 1
        )
    }
  }

  plot
}

#' Plot GWAS by group
#'
#' Plot GWAS score distributions across one grouping column.
#'
#' @param data_tibble Long cell- or sample-level GWAS score tibble containing
#'   `GWAS_ID`, `group_by_col`, and `plot_col`.
#' @param group_by_col Column used on the x axis and fill aesthetic.
#' @param plot_col Numeric score or `-log10(p)` column plotted on the y axis.
#' @param geom Distribution geometry: `boxplot` or `violin`. Boxplots include
#'   quasirandom points to show the summarized observations.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_by_group <- function(data_tibble, group_by_col = "cluster", plot_col = "score", geom = c("boxplot", "violin")) {
  if (nrow(data_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot_tibble <- if ("Category" %in% names(data_tibble)) {
    facet_levels <- data_tibble |>
      dplyr::distinct(Category, GWAS_ID) |>
      dplyr::arrange(Category, GWAS_ID) |>
      dplyr::mutate(facet_label = stringr::str_c(Category, "\n", GWAS_ID)) |>
      dplyr::pull(facet_label)

    data_tibble |>
      dplyr::mutate(facet_label = factor(stringr::str_c(Category, "\n", GWAS_ID), levels = facet_levels))
  } else {
    data_tibble |>
      dplyr::mutate(facet_label = factor(GWAS_ID, levels = sort(unique(GWAS_ID))))
  }

  plot <- plot_tibble %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[group_by_col]], y = .data[[plot_col]], fill = .data[[group_by_col]])) +
    ggplot2::facet_wrap(~facet_label) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (stringr::str_detect(plot_col, "p_val")) {
    plot <- plot + ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
  } else {
    plot <- plot + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey35")
  }

  if (geom[1] == "boxplot") {
    plot <- plot +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(
        width = 0.2,
        alpha = 0.55,
        size = 0.8,
        color = "grey20",
        show.legend = FALSE
      )
  } else if (geom[1] == "violin") {
    plot <- plot + ggplot2::geom_violin(position = ggplot2::position_dodge(width = 0.8))
  }

  return(plot)
}

plot_SCAVENGE_summary_sig_proportion <- function(summary_tibble) {
  if (nrow(summary_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot_tibble <- summary_tibble |>
    dplyr::mutate(
      facet_id = stringr::str_c(grouping_col, GWAS_ID, sep = "___"),
      group_reorder = tidytext::reorder_within(cluster, prop_sig, facet_id)
    )

  split(plot_tibble, plot_tibble$grouping_col) |>
    purrr::map(\(group_tibble) {
      ggplot2::ggplot(group_tibble, ggplot2::aes(x = group_reorder, y = prop_sig, fill = cluster)) +
        ggplot2::facet_wrap(~GWAS_ID, scales = "free") +
        tidytext::scale_x_reordered(labels = function(x) gsub("___.*$", "", x)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(x = "", y = "Proportion of enriched cells") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
    })
}

plot_SCAVENGE_summary_score_intervals <- function(summary_tibble) {
  if (nrow(summary_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  plot_tibble <- summary_tibble |>
    dplyr::mutate(
      facet_id = stringr::str_c(grouping_col, GWAS_ID, sep = "___"),
      group_reorder = tidytext::reorder_within(cluster, median_score, facet_id)
    )

  split(plot_tibble, plot_tibble$grouping_col) |>
    purrr::map(\(group_tibble) {
      ggplot2::ggplot(group_tibble, ggplot2::aes(x = group_reorder, fill = cluster)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = min_score, ymax = max_score), color = "grey35") +
        ggplot2::geom_crossbar(ggplot2::aes(y = median_score, ymin = q25_score, ymax = q75_score), width = 0.6) +
        ggplot2::facet_wrap(~GWAS_ID, scales = "free") +
        tidytext::scale_x_reordered(labels = function(x) gsub("___.*$", "", x)) +
        ggplot2::labs(x = "", y = "SCAVENGE TRS") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
    })
}

add_GWAS_heatmap_categories <- function(heatmap_data, GWAS_tibble) {
  heatmap_data |>
    dplyr::inner_join(GWAS_tibble |> dplyr::select(GWAS_ID, Category, dplyr::any_of("variant_weighting_mode")), by = "GWAS_ID")
}

#' Assign compartment
#'
#' Assign categorical compartments from regex patterns applied to type labels.
#'
#' @param in_tibble Input tibble containing `type_col`.
#' @param patterns Named character vector of regex patterns; names become
#'   compartment labels.
#' @param type_col Column whose values are matched after replacing `-` with `_`.
#' @param new_col Name of the output compartment column.
#' @return `in_tibble` with `new_col` added; unmatched rows receive `NA`.
#' @keywords internal

assign_compartment <- function(in_tibble, patterns, type_col, new_col = "compartment") {
  match_col <- paste0(type_col, "__compartment_match")
  patterns <- stringr::str_replace_all(patterns, "-", "_")
  conds <- purrr::map2(
    patterns,
    names(patterns),
    ~ rlang::expr(stringr::str_detect(.data[[match_col]], !!.x) ~ !!.y)
  )

  in_tibble |>
    dplyr::mutate(!!match_col := stringr::str_replace_all(as.character(.data[[type_col]]), "-", "_")) |>
    dplyr::mutate(
      !!new_col := dplyr::case_when(
        !!!conds,
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::select(-dplyr::all_of(match_col))
}

#' Get compartment metadata
#'
#' Build one-row-per-value metadata with ordered compartment labels.
#'
#' @param values Character/factor values to de-duplicate into metadata rows.
#' @param compartments_patterns Optional named regex vector used by
#'   `assign_compartment()`.
#' @param type_col Name of the value column in the returned metadata.
#' @param default_compartment Compartment label used when no pattern mapping is
#'   supplied.
#' @return Metadata tibble with `type_col` and factor `compartment`.
#' @keywords internal

get_compartment_metadata <- function(values, compartments_patterns, type_col, default_compartment) {
  metadata <- tibble::tibble(!!type_col := as.character(values)) |>
    dplyr::distinct()

  if (is.null(compartments_patterns)) {
    return(metadata |> dplyr::mutate(compartment = default_compartment))
  }

  metadata |>
    assign_compartment(compartments_patterns, type_col = type_col) |>
    dplyr::mutate(
      compartment = dplyr::coalesce(compartment, "Other"),
      compartment = factor(compartment, levels = c(names(compartments_patterns), "Other"))
    )
}

make_named_heatmap_palette <- function(values, palette = "Set3") {
  values <- sort(unique(stats::na.omit(values)))
  if (length(values) == 0) {
    return(character())
  }

  max_brewer_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  colors <- if (length(values) <= max_brewer_n) {
    RColorBrewer::brewer.pal(max(3, length(values)), palette)[seq_along(values)]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(max_brewer_n, palette))(length(values))
  }
  rlang::set_names(colors, values)
}

get_heatmap_legend_ncol <- function(values, max_row_chars = 70) {
  values <- sort(unique(stats::na.omit(as.character(values))))
  if (length(values) == 0) {
    return(1L)
  }

  label_widths <- nchar(values) + 6
  for (ncol in seq.int(length(values), 1L)) {
    nrow <- ceiling(length(values) / ncol)
    row_widths <- vapply(seq_len(nrow), \(row_idx) sum(label_widths[seq(row_idx, length(values), by = nrow)]), numeric(1))
    if (max(row_widths) <= max_row_chars) {
      return(ncol)
    }
  }
  1L
}

gwas_heatmap_metadata_theme <- function(show_y = FALSE, show_x = FALSE) {
  ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = if (show_x) ggplot2::element_text(angle = 45, hjust = 1) else ggplot2::element_blank(),
      axis.text.y = if (show_y) ggplot2::element_text(size = 7) else ggplot2::element_blank(),
      axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(t = 10, r = 8, b = 2, l = 2),
      strip.clip = "off",
      strip.placement = "outside",
      strip.text.x = ggplot2::element_text(angle = 45, hjust = 0, size = 8, margin = ggplot2::margin(b = 2)),
      strip.text.x.bottom = ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2))
    )
}

#' Order GWAS score plot ids
#'
#' Order GWAS rows and feature columns for clustered score heatmaps.
#'
#' @param score_data GWAS-by-feature score tibble with GWAS ID, cluster/feature,
#'   and score columns.
#' @param row_metadata GWAS metadata tibble used to annotate and order GWAS rows.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param score_col Numeric score column used for clustering/order calculation.
#' @param compartments_patterns Optional named regex patterns used to group
#'   feature columns into compartments before plotting.
#' @return List with ordered GWAS IDs, ordered feature IDs, and feature metadata.
#' @keywords internal

order_GWAS_score_plot_ids <- function(score_data, row_metadata, cluster_col = "cluster", score_col = "median_score", compartments_patterns = NULL) {
  score_mat <- score_data |>
    dplyr::select(GWAS_ID, dplyr::all_of(c(cluster_col, score_col))) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(cluster_col), values_from = dplyr::all_of(score_col), values_fill = 0) |>
    tibble::column_to_rownames("GWAS_ID") |>
    as.matrix()

  cluster_metadata <- get_compartment_metadata(
    colnames(score_mat),
    compartments_patterns,
    type_col = cluster_col,
    default_compartment = "Cluster"
  )

  cluster_order <- cluster_metadata |>
    dplyr::arrange(compartment, .data[[cluster_col]]) |>
    dplyr::group_split(compartment, .keep = FALSE) |>
    purrr::map(\(group_df) {
      clusters <- group_df[[cluster_col]]
      if (length(clusters) > 2) {
        clusters[stats::hclust(stats::dist(t(score_mat[, clusters, drop = FALSE])))$order]
      } else {
        clusters
      }
    }) |>
    unlist(use.names = FALSE)

  row_order <- row_metadata |>
    dplyr::arrange(Category, GWAS_ID) |>
    dplyr::pull(GWAS_ID)

  list(row_order = row_order, cluster_order = cluster_order, cluster_metadata = cluster_metadata)
}

get_ordered_GWAS_score_plot_data <- function(data_per_GWAS_and_cluster_df, compartments_patterns, score_col = "median_score") {
  row_metadata <- data_per_GWAS_and_cluster_df |>
    dplyr::distinct(Category, GWAS_ID) |>
    dplyr::distinct(GWAS_ID, .keep_all = TRUE)
  axes <- order_GWAS_score_plot_ids(
    data_per_GWAS_and_cluster_df,
    row_metadata,
    score_col = score_col,
    compartments_patterns = compartments_patterns
  )
  row_levels <- rev(axes$row_order)
  cluster_levels <- axes$cluster_order
  cluster_support <- data_per_GWAS_and_cluster_df |>
    dplyr::select(cluster, dplyr::any_of(c("n_cells", "n_counts", "n_features", "counts_per_feature"))) |>
    dplyr::distinct(cluster, .keep_all = TRUE)

  list(
    metadata = row_metadata |> dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = row_levels)) |> dplyr::arrange(GWAS_ID),
    scores = data_per_GWAS_and_cluster_df |> dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = row_levels), cluster = factor(cluster, levels = cluster_levels)),
    clusters = axes$cluster_metadata |>
      dplyr::left_join(cluster_support, by = "cluster") |>
      dplyr::mutate(cluster = factor(cluster, levels = cluster_levels)),
    row_levels = row_levels,
    cluster_levels = cluster_levels
  )
}

get_plot_group_breaks <- function(ordered_values) {
  group_lengths <- rle(as.character(ordered_values))$lengths
  if (length(group_lengths) <= 1) {
    return(numeric())
  }
  cumsum(group_lengths)[seq_len(length(group_lengths) - 1)] + 0.5
}

#' Plot GWAS feature heatmap
#'
#' Draw a GWAS-by-feature heatmap with score mapped directly to cell fill.
#'
#' @param score_plot_data Ordered score tibble containing GWAS IDs, feature IDs,
#'   and fill values.
#' @param feature_metadata Metadata for plotted features, including compartment
#'   ordering used for vertical separators.
#' @param feature_col Feature/cluster column plotted on the x axis.
#' @param fill_col Numeric column mapped to tile fill.
#' @param fill_label Legend label for the fill scale.
#' @param fill_midpoint Midpoint for the diverging fill scale.
#' @param title Optional plot title.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_feature_heatmap <- function(
  score_plot_data,
  feature_metadata,
  feature_col,
  fill_col,
  fill_label,
  fill_midpoint = 0,
  title = NULL,
  support_label_col = NULL
) {
  row_categories <- score_plot_data |>
    dplyr::distinct(GWAS_ID, Category) |>
    dplyr::mutate(GWAS_ID = as.character(GWAS_ID))
  row_breaks <- get_plot_group_breaks(row_categories$Category[match(levels(score_plot_data$GWAS_ID), row_categories$GWAS_ID)])
  feature_breaks <- feature_metadata |>
    dplyr::arrange(.data[[feature_col]]) |>
    dplyr::pull(compartment) |>
    get_plot_group_breaks()

  heatmap <- score_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[feature_col]], y = GWAS_ID, fill = .data[[fill_col]])) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.15)

  if (!is.null(support_label_col) && support_label_col %in% colnames(score_plot_data)) {
    heatmap <- heatmap +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[[support_label_col]]),
        size = 2.8,
        color = "grey10",
        na.rm = TRUE
      )
  }

  heatmap +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = feature_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = fill_midpoint,
      name = fill_label,
      guide = ggplot2::guide_colorbar(title.position = "top")
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.justification = "right",
      legend.box.just = "right",
      strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 2))
    )
}

plot_GWAS_feature_support_tracks <- function(feature_metadata, feature_col = "cluster") {
  if (!"n_cells" %in% colnames(feature_metadata) || all(is.na(feature_metadata$n_cells))) {
    return(patchwork::plot_spacer())
  }

  feature_metadata |>
    dplyr::arrange(.data[[feature_col]]) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[feature_col]], y = n_cells)) +
    ggplot2::geom_col(fill = "grey45", width = 0.8, na.rm = TRUE) +
    ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_y_log10(
      position = "right",
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::annotation_logticks(sides = "r") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = "Nuclei") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(angle = 90),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 10, r = 8, b = 0, l = 2)
    )
}

#' Plot GWAS metadata tracks
#'
#' Plot GWAS category, finemapping method, sample-size, loci, and ancestry tracks.
#'
#' @param ordered_metadata GWAS metadata tibble already ordered/factored by
#'   `GWAS_ID`, with category, method, loci, sample-size, and ancestry columns.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_GWAS_metadata_tracks <- function(ordered_metadata) {
  method_colors <- c("SuSie" = "#238B45", "SuSiE-inf" = "#41B6C4", "PICS" = "#F16913")
  endpoint_breaks <- function(limits) {
    limits <- limits[is.finite(limits)]
    if (length(limits) == 0) {
      return(numeric())
    }
    unique(range(limits))
  }
  row_levels <- if (is.factor(ordered_metadata$GWAS_ID)) levels(ordered_metadata$GWAS_ID) else unique(as.character(ordered_metadata$GWAS_ID))
  row_categories <- ordered_metadata |>
    dplyr::distinct(GWAS_ID, Category) |>
    dplyr::mutate(GWAS_ID = as.character(GWAS_ID))
  row_breaks <- get_plot_group_breaks(row_categories$Category[match(row_levels, row_categories$GWAS_ID)])
  bar_plot_data <- ordered_metadata |>
    dplyr::transmute(GWAS_ID, Loci = n_credible_set_loci, Samples = sample_size) |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "track", values_to = "value")
  ancestry_plot_data <- ordered_metadata |>
    dplyr::select(GWAS_ID, dplyr::matches("^ancestry_(EUR|EAS|AFR|AMR|SAS|OTH)$")) |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "ancestry_group", values_to = "fraction") |>
    dplyr::mutate(ancestry_group = stringr::str_remove(ancestry_group, "^ancestry_"))

  category_plot <- ordered_metadata |>
    ggplot2::ggplot(ggplot2::aes(x = 1, y = GWAS_ID, fill = Category)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = make_named_heatmap_palette(ordered_metadata$Category, "Set3"),
      name = "Category",
      guide = ggplot2::guide_legend(ncol = 1, title.position = "top")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme(show_y = TRUE)

  method_plot <- ordered_metadata |>
    ggplot2::ggplot(ggplot2::aes(x = 1, y = GWAS_ID, fill = finemappingMethod)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = method_colors,
      name = "Method",
      guide = ggplot2::guide_legend(ncol = 1, title.position = "top")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme()

  bar_plot <- bar_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = value, y = GWAS_ID)) +
    ggplot2::geom_col(fill = "grey45", width = 0.8, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::facet_grid(. ~ track, scales = "free_x", switch = "x") +
    ggplot2::scale_x_continuous(breaks = endpoint_breaks, labels = scales::label_number(scale_cut = scales::cut_short_scale()), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme(show_x = TRUE)

  ancestry_plot <- ancestry_plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = fraction, y = GWAS_ID, fill = ancestry_group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = row_breaks, color = "grey25", linewidth = 0.35) +
    ggplot2::scale_x_continuous(breaks = c(0, 1), labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = c("EUR" = "#4DAF4A", "EAS" = "#377EB8", "AFR" = "#984EA3", "AMR" = "#FF7F00", "SAS" = "#E41A1C", "OTH" = "#999999"),
      name = "Ancestry",
      guide = ggplot2::guide_legend(ncol = min(2L, get_heatmap_legend_ncol(ancestry_plot_data$ancestry_group, max_row_chars = 35)), byrow = TRUE, title.position = "top")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    gwas_heatmap_metadata_theme(show_x = TRUE)

  patchwork::wrap_plots(category_plot, method_plot, bar_plot, ancestry_plot, nrow = 1, widths = c(0.55, 0.55, 1.4, 1.45), guides = "collect") &
    ggplot2::theme(legend.justification = "left", legend.box.just = "left")
}

#' Plot GWAS by cluster heatmap
#'
#' Combine GWAS metadata tracks with the ordered cluster score heatmap.
#'
#' @param data_per_GWAS_and_cluster_df Score summary tibble with one row per
#'   GWAS/cluster combination.
#' @param GWAS_metadata_tracks_plot Patchwork/ggplot metadata track aligned to
#'   the same GWAS ordering.
#' @param compartments_patterns Optional named regex patterns used to group
#'   clusters into compartments.
#' @param scaled Logical; when `TRUE`, use 0.5 as the diverging color midpoint
#'   for min-max scaled scores.
#' @param fill_col Numeric column mapped to heatmap fill.
#' @param fill_label Legend label for the heatmap fill.
#' @param support_label_col Optional text column drawn on top of heatmap tiles.
#' @param show_feature_support Logical; when `TRUE`, draw a nuclei-count support
#'   annotation if `n_cells` is available.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_GWAS_by_cluster_heatmap <- function(
  data_per_GWAS_and_cluster_df,
  GWAS_metadata_tracks_plot,
  compartments_patterns = NULL,
  scaled = FALSE,
  fill_col = "median_score",
  fill_label = "Score",
  support_label_col = NULL,
  show_feature_support = TRUE
) {
  if (nrow(data_per_GWAS_and_cluster_df) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  ordered_data <- get_ordered_GWAS_score_plot_data(
    data_per_GWAS_and_cluster_df,
    compartments_patterns,
    score_col = fill_col
  )

  heatmap <- plot_GWAS_feature_heatmap(
    score_plot_data = ordered_data$scores,
    feature_metadata = ordered_data$clusters,
    feature_col = "cluster",
    fill_col = fill_col,
    fill_label = fill_label,
    fill_midpoint = if (scaled) 0.5 else 0,
    support_label_col = support_label_col
  )

  if (isTRUE(show_feature_support) && "n_cells" %in% colnames(ordered_data$clusters) && any(!is.na(ordered_data$clusters$n_cells))) {
    left_panel <- patchwork::plot_spacer() / GWAS_metadata_tracks_plot + patchwork::plot_layout(heights = c(0.18, 1))
    right_panel <- plot_GWAS_feature_support_tracks(ordered_data$clusters) / heatmap + patchwork::plot_layout(heights = c(0.18, 1))
    return(patchwork::wrap_plots(left_panel, right_panel, nrow = 1, widths = c(5.8, 9)))
  }

  patchwork::wrap_plots(GWAS_metadata_tracks_plot, heatmap, nrow = 1, widths = c(5.8, 9))
}

#' Plot grouped GWAS by cluster heatmaps
#'
#' Split GWAS cluster heatmaps by a grouping column and return a named plot list.
#'
#' @param data_per_GWAS_and_cluster_df Score summary tibble containing `split_col`
#'   in addition to GWAS and cluster score fields.
#' @param split_col Column used to split the data into one heatmap per value.
#' @param GWAS_metadata_tracks_plot Patchwork/ggplot metadata track aligned to
#'   the same GWAS ordering.
#' @param name_suffix Optional suffix appended to names of returned plot-list
#'   elements.
#' @param compartments_patterns Optional named regex patterns used to group
#'   clusters into compartments.
#' @param scaled Logical; when `TRUE`, use scaled-score color midpoint behavior.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_grouped_GWAS_by_cluster_heatmaps <- function(
  data_per_GWAS_and_cluster_df,
  split_col,
  GWAS_metadata_tracks_plot,
  name_suffix = NULL,
  compartments_patterns = NULL,
  scaled = FALSE
) {
  if (nrow(data_per_GWAS_and_cluster_df) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  grouped_heatmap_data <- dplyr::group_by(data_per_GWAS_and_cluster_df, .data[[split_col]])
  plot_names <- dplyr::group_keys(grouped_heatmap_data)[[split_col]]
  if (!is.null(name_suffix)) {
    plot_names <- stringr::str_c(plot_names, "_", name_suffix)
  }

  grouped_heatmap_data |>
    dplyr::group_split() |>
    purrr::set_names(plot_names) |>
    purrr::map(\(group_data) {
      plot_GWAS_by_cluster_heatmap(
        group_data |> dplyr::select(-dplyr::all_of(split_col)),
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = compartments_patterns,
        scaled = scaled
      )
    })
}

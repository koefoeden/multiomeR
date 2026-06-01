Michaelis_Menten_model <- function(x, k, K_m) {
  unique_reads <- (k * x) / (x + K_m)
  return(unique_reads)
}


poisson_capture_model <- function(n, C, r = 1) {
  C * (1 - exp(-n * r / C))
}


binom_model <- function(x, k) {
  frac <- (k - 1) / k
  E_X <- k * (1 - frac^x)
  return(E_X)
}


#' Fit a saturation model and extrapolate beyond current depth.
#'
#' Fits either a single-pool Poisson capture model or a Michaelis-Menten model
#' to observed (reads_per_cell, unique_reads_per_cell) data, then evaluates the
#' fitted curve on an evenly-spaced grid extending to `max_predict_reads`.
#' Returns a model-agnostic `predict_fn` closure so downstream code does not
#' need to know which model was used.
#'
#' @param observed_reads_tibble Tibble with columns `reads_per_cell` and `unique_reads_per_cell`.
#' @param model One of `"poisson"` (Lander-Waterman with efficiency, 2 parameters: C and r) or `"menten"` (Michaelis-Menten, 2 parameters: k and K_m).
#' @param max_predict_reads Upper limit of the extrapolation grid in reads/cell.
#' @param n_extrap_points Number of evenly-spaced points in the extrapolation grid.
#' @param predict_increment_reads Increment in reads/cell for the extrapolation grid.
#' @param current_unique_reads Current median unique fragments per cell (used to compute `relative_unique_gain`).
#' @return A list with elements `model_name`, `coefs`, `predict_fn`, `nls_fit`, and
#'   `saturation_tibble` (a tibble with `type` column: `"observed"`, `"fitted"` within observed range, or `"extrapolated"` beyond).
fit_saturation_model <- function(
  observed_reads_tibble,
  current_reads,
  current_unique_reads,
  model = c("poisson", "menten"),
  max_predict_reads = 2e9,
  predict_increment_reads = 1e6
) {
  model <- match.arg(model)

  if (model == "poisson") {
    nls_fit <- nls.multstart::nls_multstart(
      formula = unique_reads ~ poisson_capture_model(reads, C, r),
      data = observed_reads_tibble,
      iter = 500,
      start_lower = c(C = current_unique_reads * 1.01, r = 0.01),
      start_upper = c(C = current_unique_reads * 100, r = 1)
    )
    coefs <- list(C = stats::coef(nls_fit)[["C"]], r = stats::coef(nls_fit)[["r"]])
    predict_fn <- function(n) poisson_capture_model(n, coefs$C, coefs$r)
  } else {
    nls_fit <- nls.multstart::nls_multstart(
      formula = unique_reads ~ Michaelis_Menten_model(reads, k, K_m),
      data = observed_reads_tibble,
      iter = 250,
      start_lower = c(k = 1e3, K_m = 1e3),
      start_upper = c(k = 1e5, K_m = 1e5)
    )
    coefs <- list(k = stats::coef(nls_fit)[["k"]], K_m = stats::coef(nls_fit)[["K_m"]])
    predict_fn <- function(n) Michaelis_Menten_model(n, coefs$k, coefs$K_m)
  }

  fitted_reads_tibble <- tibble::tibble(
    reads = seq(0, current_reads, by = predict_increment_reads),
    unique_reads = predict_fn(reads),
    delta_unique_reads = unique_reads - dplyr::lag(unique_reads, default = predict_fn(current_reads) - predict_fn(current_reads - predict_increment_reads))
  ) %>%
    dplyr::mutate(type = "fitted")

  extrapolated_reads_tibble <- tibble::tibble(
    reads = seq(current_reads, max_predict_reads, by = predict_increment_reads),
    unique_reads = predict_fn(reads),
    delta_unique_reads = unique_reads - dplyr::lag(unique_reads, default = predict_fn(current_reads - predict_increment_reads))
  ) %>%
    dplyr::mutate(type = "extrapolated")

  saturation_tibble <- dplyr::bind_rows(
    fitted_reads_tibble,
    extrapolated_reads_tibble,
    observed_reads_tibble %>% dplyr::mutate(type = "observed")
  ) %>%
    dplyr::mutate(
      marginal_discovery_rate = delta_unique_reads / predict_increment_reads,
      relative_unique_gain = delta_unique_reads / current_unique_reads
    )

  list(
    model_name = model,
    coefs = coefs,
    predict_fn = predict_fn,
    nls_fit = nls_fit,
    saturation_tibble = saturation_tibble
  )
}


#' Get sequencing allocation tibble
#'
#' Allocate extra sequencing reads across reactions from saturation-model predictions.
#'
#' @param combined_saturation_per_reaction_tibble One row per reaction with
#'   current read/unique counts, TSS enrichment, donor count, fitted prediction
#'   function, and nested extrapolated saturation steps.
#' @param predict_increment_reads Read increment represented by each candidate
#'   extrapolation step and by each budget allocation step.
#' @param sequencing_budget Total extra read budget available across reactions.
#' @param TSS_enrichment_threshold TSS enrichment value at which the quality
#'   weight is capped at 1; lower values reduce priority scores.
#' @param min_predicted_extra_unique_threshold Minimum predicted extra unique
#'   reads required for a reaction to keep allocated budget.
#' @param max_allocation_passes Maximum reallocation passes after dropping
#'   reactions below the predicted-gain threshold.
#' @return A tibble with per-reaction allocated reads, predicted final reads and
#'   unique reads, allocation fraction, and `over_threshold` status.
#' @keywords internal

get_sequencing_allocation_tibble <- function(
  combined_saturation_per_reaction_tibble,
  predict_increment_reads,
  sequencing_budget,
  TSS_enrichment_threshold,
  min_predicted_extra_unique_threshold,
  max_allocation_passes = 10
) {
  reaction_info_tibble <- combined_saturation_per_reaction_tibble %>%
    dplyr::distinct(map_TENX_reaction_ID, .keep_all = TRUE)

  pooled_steps_tibble <- reaction_info_tibble %>%
    tidyr::unnest(saturation_tibble) %>%
    dplyr::filter(type == "extrapolated")

  ranked_steps_tibble <- pooled_steps_tibble %>%
    dplyr::mutate(
      quality_weight = pmin(tss_enrichment / TSS_enrichment_threshold, 1),
      priority_score = marginal_discovery_rate *
        quality_weight *
        relative_unique_gain *
        n_donors
    ) %>%
    dplyr::arrange(dplyr::desc(priority_score))

  allocate_once <- function(candidate_steps_tibble) {
    allocated_steps_tibble <- candidate_steps_tibble %>%
      dplyr::mutate(
        step_reads = predict_increment_reads,
        cumulative_reads = cumsum(step_reads),
        cumulative_unique_reads = cumsum(delta_unique_reads)
      ) %>%
      dplyr::filter(cumulative_reads <= sequencing_budget)

    allocated_steps_tibble %>%
      dplyr::summarise(
        allocated_extra_reads = sum(step_reads),
        .by = map_TENX_reaction_ID
      ) %>%
      dplyr::left_join(reaction_info_tibble, by = "map_TENX_reaction_ID") %>%
      dplyr::mutate(
        predicted_final_reads = current_reads + allocated_extra_reads,
        predicted_final_unique_reads = purrr::map2_dbl(
          predict_fn,
          predicted_final_reads,
          ~ .x(.y)
        ),
        predicted_extra_unique_reads = predicted_final_unique_reads - current_unique_reads
      )
  }

  kept_reactions_vec <- ranked_steps_tibble %>%
    dplyr::distinct(map_TENX_reaction_ID) %>%
    dplyr::pull(map_TENX_reaction_ID)

  allocation_summary_w_sat_tibble <- allocate_once(ranked_steps_tibble %>% dplyr::slice(0))

  # Reallocate budget after removing reactions that do not clear the final thresholds.
  for (pass_idx in seq_len(max_allocation_passes)) {
    allocation_summary_w_sat_tibble <- ranked_steps_tibble %>%
      dplyr::filter(map_TENX_reaction_ID %in% kept_reactions_vec) %>%
      allocate_once()

    next_kept_reactions_vec <- allocation_summary_w_sat_tibble %>%
      dplyr::filter(
        predicted_extra_unique_reads >= min_predicted_extra_unique_threshold
      ) %>%
      dplyr::pull(map_TENX_reaction_ID)

    if (base::setequal(next_kept_reactions_vec, kept_reactions_vec)) {
      kept_reactions_vec <- next_kept_reactions_vec
      break
    }

    kept_reactions_vec <- next_kept_reactions_vec
  }

  allocation_summary_w_sat_tibble <- allocation_summary_w_sat_tibble %>%
    dplyr::filter(map_TENX_reaction_ID %in% kept_reactions_vec) %>%
    dplyr::mutate(
      allocation_fraction = allocated_extra_reads / sum(allocated_extra_reads)
    )

  unallocated_reactions_vec <- setdiff(reaction_info_tibble$map_TENX_reaction_ID, allocation_summary_w_sat_tibble$map_TENX_reaction_ID)
  unallocated_summary_w_sat_tibble <- reaction_info_tibble %>%
    dplyr::filter(map_TENX_reaction_ID %in% unallocated_reactions_vec) %>%
    dplyr::mutate(
      over_threshold = FALSE,
      allocated_extra_reads = 0,
      allocation_fraction = 0,
      predicted_extra_unique_reads = 0,
      predicted_final_reads = current_reads,
      predicted_final_unique_reads = current_unique_reads
    )

  allocation_summary_w_sat_tibble %>%
    dplyr::mutate(over_threshold = TRUE) %>%
    dplyr::bind_rows(
      unallocated_summary_w_sat_tibble
    )
}

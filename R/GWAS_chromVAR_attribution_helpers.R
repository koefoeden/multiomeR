#' Map credible-set variant weights to ATAC peaks
#'
#' Allocate each capped trait-level peak weight back to its overlapping
#' credible-set variants in proportion to their configured posterior weights.
#'
#' @param GWAS_input_record One trait-level GWAS input record.
#' @param peak_ranges Ordered ATAC peak ranges.
#' @param posterior_probability_cutoff Minimum posterior probability retained.
#' @param posterior_probability_weighting_function Optional variant weighting function.
#' @return One row per peak-variant overlap with weights summing to the capped
#'   trait-level weight for each peak.
#' @keywords internal

get_GWAS_chromVAR_peak_variant_weight_tibble <- function(
  GWAS_input_record,
  peak_ranges,
  posterior_probability_cutoff = NULL,
  posterior_probability_weighting_function = NULL
) {
  variant_GRanges <- GWAS_input_record$credible_set_GRanges |>
    filter_credible_set_variants(
      posterior_probability_cutoff = posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    )
  overlap_hits <- GenomicRanges::findOverlaps(peak_ranges, variant_GRanges)
  if (length(overlap_hits) == 0L) {
    stop("No credible-set variants overlap ATAC peaks for GWAS_ID: ", GWAS_input_record$GWAS_ID)
  }

  peak_idx <- S4Vectors::queryHits(overlap_hits)
  variant_idx <- S4Vectors::subjectHits(overlap_hits)
  peak_tibble <- tibble::tibble(
    peak_name = get_peak_names_from_GRanges(peak_ranges)[peak_idx],
    peak_chromosome = as.character(GenomicRanges::seqnames(peak_ranges))[peak_idx],
    peak_start = GenomicRanges::start(peak_ranges)[peak_idx],
    peak_end = GenomicRanges::end(peak_ranges)[peak_idx]
  )
  variant_tibble <- GenomicRanges::mcols(variant_GRanges)[variant_idx, , drop = FALSE] |>
    as.data.frame() |>
    tibble::as_tibble()

  dplyr::bind_cols(peak_tibble, variant_tibble) |>
    dplyr::mutate(GWAS_ID = GWAS_input_record$GWAS_ID, .before = 1) |>
    dplyr::mutate(
      uncapped_peak_weight = sum(.data$posteriorProbability),
      peak_weight = pmin(.data$uncapped_peak_weight, 1),
      peak_weight_scale = .data$peak_weight / .data$uncapped_peak_weight,
      peak_variant_weight = .data$posteriorProbability * .data$peak_weight_scale,
      .by = c(GWAS_ID, peak_name)
    ) |>
    dplyr::select(
      GWAS_ID,
      peak_name,
      peak_chromosome,
      peak_start,
      peak_end,
      peak_weight,
      uncapped_peak_weight,
      peak_variant_weight,
      peak_weight_scale,
      dplyr::everything()
    )
}

#' Decompose cell-type chromVAR deviations into peak contributions
#'
#' @param chromVAR_deviation_record Cell-type pseudobulk deviation record.
#' @param psbulk_ATAC_data_matrix Peak-by-cell-type count matrix.
#' @param chromVAR_obj Template chromVAR object with peak ranges.
#' @param annotation_matrix Peak-by-GWAS trait weight matrix.
#' @param GWAS_inputs_tibble GWAS metadata.
#' @return Long tibble whose peak contributions sum to the raw deviation,
#'   relative deviation, and z score for each GWAS and cell type.
#' @keywords internal

get_GWAS_chromVAR_peak_attribution_tibble <- function(
  chromVAR_deviation_record,
  psbulk_ATAC_data_matrix,
  chromVAR_obj,
  annotation_matrix,
  GWAS_inputs_tibble
) {
  peak_names <- chromVAR_deviation_record$peak_names
  counts_matrix <- psbulk_ATAC_data_matrix[peak_names, , drop = FALSE]
  if (inherits(counts_matrix, "IterableMatrix")) {
    counts_matrix <- methods::as(counts_matrix, "dgCMatrix")
  }
  annotation_matrix <- annotation_matrix[peak_names, , drop = FALSE]
  deviation_SE <- chromVAR_deviation_record$deviation_SE
  background <- chromVAR_deviation_record$background
  peak_ranges <- SummarizedExperiment::rowRanges(chromVAR_obj)
  peak_ranges <- peak_ranges[match(peak_names, get_peak_names_from_GRanges(peak_ranges))]
  peak_metadata_tibble <- GenomicRanges::as.data.frame(peak_ranges) |>
    tibble::as_tibble() |>
    dplyr::transmute(
      peak_name = peak_names,
      peak_chromosome = as.character(.data$seqnames),
      peak_start = .data$start,
      peak_end = .data$end
    )

  GWAS_IDs <- rownames(SummarizedExperiment::assay(deviation_SE, "deviations"))
  cluster_names <- colnames(counts_matrix)
  background_expectation_matrix <- methods::slot(background, "E")
  background_variance_matrix <- methods::slot(background, "V")
  peak_background_bin_idx <- methods::slot(background, "peak2bin")
  peak_expectation_vec <- methods::slot(background, "expectation")
  cluster_depth_vec <- methods::slot(background, "depth")

  purrr::map_dfr(GWAS_IDs, \(GWAS_ID) {
    peak_weight_vec <- as.numeric(annotation_matrix[, GWAS_ID])
    weighted_peak_idx <- which(peak_weight_vec != 0)
    peak_weight_vec <- peak_weight_vec[weighted_peak_idx]
    weighted_peak_bin_idx <- peak_background_bin_idx[weighted_peak_idx]
    observed_count_matrix <- as.matrix(counts_matrix[weighted_peak_idx, , drop = FALSE])
    expected_count_matrix <- background_expectation_matrix[weighted_peak_bin_idx, , drop = FALSE]
    deviation_numerator_matrix <- sweep(
      observed_count_matrix - expected_count_matrix,
      MARGIN = 1,
      STATS = peak_weight_vec,
      FUN = "*"
    )

    global_weighted_expectation <- sum(peak_weight_vec * peak_expectation_vec[weighted_peak_idx])
    deviation_denominator_vec <- global_weighted_expectation * cluster_depth_vec / sum(peak_expectation_vec)
    deviation_contribution_matrix <- sweep(
      deviation_numerator_matrix,
      MARGIN = 2,
      STATS = deviation_denominator_vec,
      FUN = "/"
    )

    bin_weight_vec <- Matrix::sparseMatrix(
      i = weighted_peak_bin_idx,
      j = rep.int(1L, length(weighted_peak_bin_idx)),
      x = peak_weight_vec,
      dims = c(nrow(background_variance_matrix), 1L)
    ) |>
      Matrix::rowSums()
    z_denominator_vec <- sqrt(as.numeric(Matrix::crossprod(bin_weight_vec, background_variance_matrix)))
    z_contribution_matrix <- sweep(
      deviation_numerator_matrix,
      MARGIN = 2,
      STATS = z_denominator_vec,
      FUN = "/"
    )

    deviation_vec <- SummarizedExperiment::assay(deviation_SE, "deviations")[GWAS_ID, ]
    z_vec <- SummarizedExperiment::assay(deviation_SE, "z")[GWAS_ID, ]
    deviation_sd <- stats::sd(deviation_vec)
    relative_deviation_vec <- scale_within_vector(deviation_vec)
    relative_contribution_matrix <- if (is.na(deviation_sd) || deviation_sd == 0) {
      deviation_contribution_matrix * 0
    } else {
      sweep(deviation_contribution_matrix, 1, rowMeans(deviation_contribution_matrix), FUN = "-") / deviation_sd
    }

    n_weighted_peaks <- length(weighted_peak_idx)
    cluster_idx <- rep(seq_along(cluster_names), each = n_weighted_peaks)
    peak_idx <- rep(seq_len(n_weighted_peaks), times = length(cluster_names))
    GWAS_metadata <- GWAS_inputs_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID) |>
      dplyr::slice_head(n = 1)

    peak_metadata_tibble[weighted_peak_idx[peak_idx], ] |>
      dplyr::mutate(
        GWAS_ID = GWAS_ID,
        Category = GWAS_metadata$Category[[1]],
        variant_weighting_mode = GWAS_metadata$variant_weighting_mode[[1]],
        cluster = cluster_names[cluster_idx],
        peak_weight = peak_weight_vec[peak_idx],
        background_bin = weighted_peak_bin_idx[peak_idx],
        observed_count = as.vector(observed_count_matrix),
        background_expected_count = as.vector(expected_count_matrix),
        deviation_numerator_contribution = as.vector(deviation_numerator_matrix),
        deviation_contribution = as.vector(deviation_contribution_matrix),
        relative_deviation_contribution = as.vector(relative_contribution_matrix),
        z_contribution = as.vector(z_contribution_matrix),
        deviation = deviation_vec[cluster_idx],
        relative_deviation = relative_deviation_vec[cluster_idx],
        z = z_vec[cluster_idx],
        .before = 1
      )
  })
}

#' Allocate peak contributions to credible-set variants
#'
#' @param peak_attribution_tibble Exact peak-level chromVAR contributions.
#' @param peak_variant_weight_tibble Peak-variant weights.
#' @return Variant-level attribution rows, summed across overlapping peaks.
#' @keywords internal

get_GWAS_chromVAR_variant_attribution_tibble <- function(peak_attribution_tibble, peak_variant_weight_tibble) {
  contribution_cols <- c(
    "deviation_numerator_contribution",
    "deviation_contribution",
    "relative_deviation_contribution",
    "z_contribution"
  )

  peak_variant_attribution_tibble <- peak_attribution_tibble |>
    dplyr::rename(attribution_peak_weight = peak_weight) |>
    dplyr::inner_join(
      peak_variant_weight_tibble |>
        dplyr::select(-dplyr::any_of(c("Category", "variant_weighting_mode"))) |>
        dplyr::rename(mapped_peak_weight = peak_weight),
      by = c("GWAS_ID", "peak_name", "peak_chromosome", "peak_start", "peak_end"),
      relationship = "many-to-many"
    )
  peak_weight_error <- max(
    abs(peak_variant_attribution_tibble$attribution_peak_weight - peak_variant_attribution_tibble$mapped_peak_weight),
    na.rm = TRUE
  )
  if (peak_weight_error > 1e-10) {
    stop("Peak-variant weights do not reproduce the trait-level annotation matrix; maximum error: ", peak_weight_error)
  }

  peak_variant_attribution_tibble |>
    dplyr::mutate(
      variant_share_of_peak_weight = .data$peak_variant_weight / .data$mapped_peak_weight,
      dplyr::across(dplyr::all_of(contribution_cols), \(value) value * .data$variant_share_of_peak_weight)
    ) |>
    dplyr::summarise(
      Category = dplyr::first(.data$Category),
      variant_weighting_mode = dplyr::first(.data$variant_weighting_mode),
      credibleSetIndex = dplyr::first(.data$credibleSetIndex),
      finemappingMethod = dplyr::first(.data$finemappingMethod),
      chromosome = dplyr::first(.data$chromosome),
      position = dplyr::first(.data$position),
      locusStart = dplyr::first(.data$locusStart),
      locusEnd = dplyr::first(.data$locusEnd),
      posteriorProbability = dplyr::first(.data$posteriorProbability),
      posteriorProbability_raw = dplyr::first(.data$posteriorProbability_raw),
      peak_variant_weight = sum(.data$peak_variant_weight),
      n_peaks = dplyr::n_distinct(.data$peak_name),
      peak_names = list(sort(unique(.data$peak_name))),
      peak_start = min(.data$peak_start),
      peak_end = max(.data$peak_end),
      dplyr::across(dplyr::all_of(contribution_cols), sum),
      deviation = dplyr::first(.data$deviation),
      relative_deviation = dplyr::first(.data$relative_deviation),
      z = dplyr::first(.data$z),
      .by = c(GWAS_ID, cluster, studyLocusId, variantId)
    ) |>
    dplyr::relocate(
      GWAS_ID,
      Category,
      variant_weighting_mode,
      cluster,
      studyLocusId,
      variantId
    )
}

#' Sum credible-set variant contributions by locus
#'
#' @param variant_attribution_tibble Variant-level chromVAR attribution.
#' @return Exact locus-level contributions with lead-variant labels.
#' @keywords internal

get_GWAS_chromVAR_locus_attribution_tibble <- function(variant_attribution_tibble) {
  lead_variant_tibble <- variant_attribution_tibble |>
    dplyr::arrange(dplyr::desc(.data$posteriorProbability), .data$variantId) |>
    dplyr::slice_head(n = 1, by = c(GWAS_ID, studyLocusId)) |>
    dplyr::select(
      GWAS_ID,
      studyLocusId,
      lead_variantId = variantId,
      lead_variant_position = position,
      lead_variant_PIP = posteriorProbability
    )

  variant_attribution_tibble |>
    dplyr::summarise(
      Category = dplyr::first(.data$Category),
      variant_weighting_mode = dplyr::first(.data$variant_weighting_mode),
      chromosome = dplyr::first(.data$chromosome),
      locus_start = min(c(.data$locusStart, .data$peak_start, .data$position), na.rm = TRUE),
      locus_end = max(c(.data$locusEnd, .data$peak_end, .data$position), na.rm = TRUE),
      n_variants = dplyr::n_distinct(.data$variantId),
      n_peaks = length(unique(unlist(.data$peak_names))),
      posterior_probability_sum = sum(.data$posteriorProbability),
      deviation_numerator_contribution = sum(.data$deviation_numerator_contribution),
      deviation_contribution = sum(.data$deviation_contribution),
      relative_deviation_contribution = sum(.data$relative_deviation_contribution),
      z_contribution = sum(.data$z_contribution),
      deviation = dplyr::first(.data$deviation),
      relative_deviation = dplyr::first(.data$relative_deviation),
      z = dplyr::first(.data$z),
      .by = c(GWAS_ID, cluster, studyLocusId)
    ) |>
    dplyr::left_join(lead_variant_tibble, by = c("GWAS_ID", "studyLocusId")) |>
    dplyr::mutate(
      locus_label = .data$lead_variantId,
      attribution_rank = dplyr::min_rank(dplyr::desc(abs(.data$relative_deviation_contribution))),
      .by = c(GWAS_ID, cluster)
    ) |>
    dplyr::relocate(
      GWAS_ID,
      Category,
      variant_weighting_mode,
      cluster,
      studyLocusId,
      locus_label,
      lead_variantId
    )
}

#' Validate additive GWAS chromVAR attribution
#'
#' @param peak_attribution_tibble Peak-level attribution.
#' @param variant_attribution_tibble Variant-level attribution.
#' @param locus_attribution_tibble Locus-level attribution.
#' @return Reconciliation errors for each GWAS and cell type.
#' @keywords internal

get_GWAS_chromVAR_attribution_QC_tibble <- function(
  peak_attribution_tibble,
  variant_attribution_tibble,
  locus_attribution_tibble
) {
  summarize_level <- function(attribution_tibble, level) {
    attribution_tibble |>
      dplyr::summarise(
        deviation_contribution_sum = sum(.data$deviation_contribution),
        relative_deviation_contribution_sum = sum(.data$relative_deviation_contribution),
        z_contribution_sum = sum(.data$z_contribution),
        deviation = dplyr::first(.data$deviation),
        relative_deviation = dplyr::first(.data$relative_deviation),
        z = dplyr::first(.data$z),
        .by = c(GWAS_ID, cluster)
      ) |>
      dplyr::mutate(level = level)
  }

  dplyr::bind_rows(
    summarize_level(peak_attribution_tibble, "peak"),
    summarize_level(variant_attribution_tibble, "variant"),
    summarize_level(locus_attribution_tibble, "locus")
  ) |>
    dplyr::mutate(
      deviation_error = .data$deviation_contribution_sum - .data$deviation,
      relative_deviation_error = .data$relative_deviation_contribution_sum - .data$relative_deviation,
      z_error = .data$z_contribution_sum - .data$z,
      attribution_reconciles = abs(.data$deviation_error) < 1e-8 &
        abs(.data$relative_deviation_error) < 1e-8 &
        abs(.data$z_error) < 1e-8
    ) |>
    dplyr::relocate(GWAS_ID, cluster, level)
}

collapse_GWAS_locus_attribution_for_plot <- function(locus_attribution_tibble, n_top_loci = 12L) {
  top_locus_ids <- locus_attribution_tibble |>
    dplyr::summarise(
      max_abs_contribution = max(abs(.data$relative_deviation_contribution)),
      locus_label = dplyr::first(.data$locus_label),
      .by = studyLocusId
    ) |>
    dplyr::slice_max(.data$max_abs_contribution, n = n_top_loci, with_ties = FALSE) |>
    dplyr::arrange(dplyr::desc(.data$max_abs_contribution))

  top_tibble <- locus_attribution_tibble |>
    dplyr::semi_join(top_locus_ids, by = "studyLocusId") |>
    dplyr::mutate(plot_locus = .data$locus_label)
  other_tibble <- locus_attribution_tibble |>
    dplyr::anti_join(top_locus_ids, by = "studyLocusId") |>
    dplyr::mutate(plot_locus = dplyr::if_else(.data$relative_deviation_contribution >= 0, "Other positive", "Other negative")) |>
    dplyr::summarise(
      relative_deviation_contribution = sum(.data$relative_deviation_contribution),
      .by = c(cluster, plot_locus)
    )
  locus_order <- c(top_locus_ids$locus_label, "Other positive", "Other negative")

  dplyr::bind_rows(
    top_tibble |>
      dplyr::select(cluster, plot_locus, relative_deviation_contribution),
    other_tibble
  ) |>
    tidyr::complete(
      cluster,
      plot_locus = locus_order,
      fill = list(relative_deviation_contribution = 0)
    ) |>
    dplyr::mutate(plot_locus = factor(.data$plot_locus, levels = locus_order))
}

#' Plot cell-type-by-locus chromVAR contributions
#'
#' @param locus_attribution_tibble Exact locus-level attribution.
#' @param chromVAR_deviation_tibble Original cell-type-by-GWAS heatmap data.
#' @param n_top_loci Number of individually labelled loci per GWAS.
#' @return Named list with one contribution heatmap per GWAS.
#' @keywords internal

plot_GWAS_locus_attribution_heatmaps <- function(
  locus_attribution_tibble,
  chromVAR_deviation_tibble,
  n_top_loci = 12L
) {
  GWAS_IDs <- unique(locus_attribution_tibble$GWAS_ID)
  plots <- purrr::map(GWAS_IDs, \(GWAS_ID) {
    locus_tibble <- locus_attribution_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID)
    heatmap_tibble <- collapse_GWAS_locus_attribution_for_plot(locus_tibble, n_top_loci = n_top_loci)
    total_tibble <- chromVAR_deviation_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID)
    cluster_levels <- rev(unique(total_tibble$cluster))
    heatmap_tibble <- heatmap_tibble |>
      dplyr::mutate(cluster = factor(.data$cluster, levels = cluster_levels))
    total_tibble <- total_tibble |>
      dplyr::mutate(cluster = factor(.data$cluster, levels = cluster_levels))

    heatmap_plot <- heatmap_tibble |>
      ggplot2::ggplot(ggplot2::aes(x = .data$plot_locus, y = .data$cluster, fill = .data$relative_deviation_contribution)) +
      ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        name = "Relative deviation\ncontribution"
      ) +
      ggplot2::labs(x = "Credible-set locus (lead variant)", y = NULL) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    total_plot <- total_tibble |>
      ggplot2::ggplot(ggplot2::aes(x = .data$relative_deviation, y = .data$cluster, fill = .data$relative_deviation >= 0)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::geom_vline(xintercept = 0, color = "grey45", linewidth = 0.3) +
      ggplot2::geom_text(ggplot2::aes(label = .data$support_label), size = 3, hjust = -0.3) +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#B40426", `FALSE` = "#3B4CC0"), guide = "none") +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::labs(x = "Total relative deviation", y = NULL) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
      )

    combined_plot <- patchwork::wrap_plots(heatmap_plot, total_plot, nrow = 1, widths = c(5, 1.2)) +
      patchwork::plot_annotation(title = GWAS_ID)
    combined_plot$labels$title <- ggplot2::waiver()
    combined_plot
  })
  names(plots) <- GWAS_IDs
  plots
}

prepare_GWAS_locus_attribution_waterfall_tibble <- function(locus_tibble, n_top_loci = 12L) {
  top_tibble <- locus_tibble |>
    dplyr::slice_max(abs(.data$relative_deviation_contribution), n = n_top_loci, with_ties = FALSE)
  other_tibble <- locus_tibble |>
    dplyr::anti_join(top_tibble, by = "studyLocusId") |>
    dplyr::mutate(locus_label = dplyr::if_else(.data$relative_deviation_contribution >= 0, "Other positive", "Other negative")) |>
    dplyr::summarise(
      relative_deviation_contribution = sum(.data$relative_deviation_contribution),
      .by = locus_label
    )
  step_tibble <- dplyr::bind_rows(
    top_tibble |>
      dplyr::select(locus_label, relative_deviation_contribution),
    other_tibble
  ) |>
    dplyr::filter(.data$relative_deviation_contribution != 0) |>
    dplyr::mutate(is_other = stringr::str_starts(.data$locus_label, "Other ")) |>
    dplyr::arrange(
      dplyr::desc(.data$relative_deviation_contribution >= 0),
      .data$is_other,
      dplyr::desc(abs(.data$relative_deviation_contribution))
    ) |>
    dplyr::mutate(
      plot_index = dplyr::row_number(),
      start = dplyr::lag(cumsum(.data$relative_deviation_contribution), default = 0),
      end = .data$start + .data$relative_deviation_contribution,
      direction = dplyr::if_else(.data$relative_deviation_contribution >= 0, "Positive", "Negative")
    )
  total <- sum(step_tibble$relative_deviation_contribution)

  dplyr::bind_rows(
    step_tibble,
    tibble::tibble(
      locus_label = "Total",
      relative_deviation_contribution = total,
      plot_index = nrow(step_tibble) + 1L,
      start = 0,
      end = total,
      direction = "Total"
    )
  )
}

#' Plot locus-attribution waterfalls
#'
#' @param locus_attribution_tibble Exact locus-level attribution.
#' @param n_top_loci Number of individually labelled loci per plot.
#' @return Named list with one waterfall per GWAS and cell type.
#' @keywords internal

plot_GWAS_locus_attribution_waterfalls <- function(locus_attribution_tibble, n_top_loci = 12L) {
  split_tibbles <- split(
    locus_attribution_tibble,
    interaction(locus_attribution_tibble$GWAS_ID, locus_attribution_tibble$cluster, drop = TRUE, lex.order = TRUE)
  )

  purrr::imap(split_tibbles, \(locus_tibble, plot_name) {
    waterfall_tibble <- prepare_GWAS_locus_attribution_waterfall_tibble(locus_tibble, n_top_loci = n_top_loci)
    connector_tibble <- waterfall_tibble |>
      dplyr::filter(.data$locus_label != "Total")
    connector_tibble <- connector_tibble |>
      dplyr::slice_head(n = max(0L, nrow(connector_tibble) - 1L))

    ggplot2::ggplot(waterfall_tibble) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$plot_index - 0.42,
          xmax = .data$plot_index + 0.42,
          ymin = pmin(.data$start, .data$end),
          ymax = pmax(.data$start, .data$end),
          fill = .data$direction
        )
      ) +
      ggplot2::geom_segment(
        data = connector_tibble,
        ggplot2::aes(
          x = .data$plot_index + 0.42,
          xend = .data$plot_index + 1 - 0.42,
          y = .data$end,
          yend = .data$end
        ),
        color = "grey55",
        linewidth = 0.3
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "grey35", linewidth = 0.3) +
      ggplot2::scale_x_continuous(
        breaks = waterfall_tibble$plot_index,
        labels = waterfall_tibble$locus_label,
        expand = ggplot2::expansion(add = 0.6)
      ) +
      ggplot2::scale_fill_manual(values = c(Positive = "#B40426", Negative = "#3B4CC0", Total = "grey25"), name = NULL) +
      ggplot2::labs(
        x = "Credible-set locus (lead variant)",
        y = "Cumulative relative deviation",
        title = stringr::str_glue("{locus_tibble$GWAS_ID[[1]]} - {locus_tibble$cluster[[1]]}")
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.major.x = ggplot2::element_blank(),
        legend.position = "bottom"
      )
  }) |>
    rlang::set_names(stringr::str_replace_all(names(split_tibbles), "[^A-Za-z0-9_.-]+", "_"))
}

make_GWAS_variant_attribution_track <- function(variant_tibble, region) {
  region <- BPCells:::normalize_ranges(region)
  max_abs_contribution <- max(abs(variant_tibble$relative_deviation_contribution))
  variant_tibble <- variant_tibble |>
    dplyr::mutate(
      direction = dplyr::if_else(.data$relative_deviation_contribution >= 0, "Positive", "Negative"),
      variant_label = dplyr::if_else(
        dplyr::min_rank(dplyr::desc(abs(.data$relative_deviation_contribution))) <= 3,
        .data$variantId,
        ""
      )
    )

  BPCells:::wrap_trackplot(
    ggplot2::ggplot(variant_tibble, ggplot2::aes(x = .data$position, y = .data$relative_deviation_contribution)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data$position, y = 0, yend = .data$relative_deviation_contribution, color = .data$direction),
        linewidth = 0.4
      ) +
      ggplot2::geom_point(ggplot2::aes(size = .data$posteriorProbability, color = .data$direction)) +
      ggrepel::geom_text_repel(
        data = dplyr::filter(variant_tibble, .data$variant_label != ""),
        ggplot2::aes(label = .data$variant_label),
        size = 2,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      ggplot2::scale_color_manual(values = c(Positive = "#B40426", Negative = "#3B4CC0"), guide = "none") +
      ggplot2::scale_size_continuous(range = c(1.5, 4), name = "PIP") +
      ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels = scales::label_number()) +
      ggplot2::scale_y_continuous(limits = c(-max_abs_contribution, max_abs_contribution) * 1.15, expand = c(0, 0)) +
      ggplot2::labs(x = "Genomic position (bp)", y = "Relative deviation contribution") +
      BPCells:::trackplot_theme(),
    ggplot2::unit(2, "null"),
    region = region
  )
}

#' Plot variant-level attribution details for top loci
#'
#' @param variant_attribution_tibble Variant-level attribution.
#' @param locus_attribution_tibble Locus-level attribution used to select loci.
#' @param consensus_peak_GRanges Consensus ATAC peaks.
#' @param fragments BPCells fragment object.
#' @param metadata_tibble Cell metadata with cell-type labels.
#' @param n_top_cell_types Number of highest-deviation cell types detailed per GWAS.
#' @param n_top_loci Number of loci plotted per GWAS and cell type.
#' @param flank Bases added around the selected locus.
#' @return Named list of variant, accessibility, and peak track plots.
#' @keywords internal

plot_GWAS_variant_attribution_details <- function(
  variant_attribution_tibble,
  locus_attribution_tibble,
  consensus_peak_GRanges,
  fragments,
  metadata_tibble,
  n_top_cell_types = 1L,
  n_top_loci = 3L,
  flank = 25000L,
  group_cells_by_col = "PCA_harmony_SNN_cluster_cell_type"
) {
  selected_cell_type_tibble <- locus_attribution_tibble |>
    dplyr::distinct(.data$GWAS_ID, .data$cluster, .data$relative_deviation) |>
    dplyr::slice_max(
      .data$relative_deviation,
      n = n_top_cell_types,
      with_ties = FALSE,
      by = GWAS_ID
    )
  selected_locus_tibble <- locus_attribution_tibble |>
    dplyr::semi_join(selected_cell_type_tibble, by = c("GWAS_ID", "cluster")) |>
    dplyr::slice_max(
      abs(.data$relative_deviation_contribution),
      n = n_top_loci,
      with_ties = FALSE,
      by = c(GWAS_ID, cluster)
    ) |>
    dplyr::mutate(
      detail_rank = dplyr::row_number(),
      .by = c(GWAS_ID, cluster)
    )

  fragments <- BPCells::select_cells(fragments, metadata_tibble$barcode_w_prefix)
  fragment_cell_names <- BPCells::cellNames(fragments)
  metadata <- metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(
      .data$barcode_w_prefix %in% fragment_cell_names,
      !is.na(.data[[group_cells_by_col]])
    ) |>
    dplyr::arrange(match(.data$barcode_w_prefix, fragment_cell_names))
  fragments <- BPCells::select_cells(fragments, metadata$barcode_w_prefix)
  cell_read_counts <- if ("atac_fragments" %in% colnames(metadata)) metadata$atac_fragments else metadata$nCount_ATAC
  group_values <- stringr::str_replace_all(as.character(metadata[[group_cells_by_col]]), "_", "-")

  locus_tibble_list <- selected_locus_tibble |>
    dplyr::group_by(.data$GWAS_ID, .data$cluster, .data$studyLocusId, .data$detail_rank) |>
    dplyr::group_split()
  plots <- locus_tibble_list |>
    purrr::map(\(locus_tibble) {
      GWAS_ID <- locus_tibble$GWAS_ID[[1]]
      cluster <- locus_tibble$cluster[[1]]
      studyLocusId <- locus_tibble$studyLocusId[[1]]
      variant_tibble <- variant_attribution_tibble |>
        dplyr::filter(
          .data$GWAS_ID == .env$GWAS_ID,
          .data$cluster == .env$cluster,
          .data$studyLocusId == .env$studyLocusId
        )
      region <- GenomicRanges::GRanges(
        seqnames = stringr::str_c("chr", locus_tibble$chromosome[[1]]),
        ranges = IRanges::IRanges(
          start = max(1L, locus_tibble$locus_start[[1]] - flank),
          end = locus_tibble$locus_end[[1]] + flank
        )
      )
      group_levels <- c(cluster, setdiff(gtools::mixedsort(unique(group_values)), cluster))
      coverage_tibble <- BPCells::trackplot_coverage(
        fragments = fragments,
        region = region,
        groups = factor(group_values, levels = group_levels),
        cell_read_counts = cell_read_counts,
        group_order = group_levels,
        bins = 500,
        return_data = TRUE
      )
      coverage_colors <- c("#B40426", rep("grey75", length(group_levels) - 1L))

      plot_title <- stringr::str_glue(
        "{GWAS_ID} - {cluster} - {locus_tibble$locus_label[[1]]}; ",
        "locus contribution={round(locus_tibble$relative_deviation_contribution[[1]], 3)}"
      )
      detail_plot <- BPCells::trackplot_combine(
        tracks = list(
          make_GWAS_variant_attribution_track(variant_tibble, region),
          make_BPCells_ATAC_coverage_track_from_tibble(coverage_tibble, region, colors = coverage_colors),
          make_consensus_peak_locus_track(consensus_peak_GRanges, region)
        ),
        title = plot_title
      )
      detail_plot$labels$title <- ggplot2::waiver()
      detail_plot
    })
  plot_names <- purrr::map_chr(locus_tibble_list, \(locus_tibble) {
    stringr::str_c(
      locus_tibble$GWAS_ID[[1]],
      locus_tibble$cluster[[1]],
      sprintf("%02d", locus_tibble$detail_rank[[1]]),
      locus_tibble$locus_label[[1]],
      sep = "__"
    ) |>
      stringr::str_replace_all("[^A-Za-z0-9_.-]+", "_")
  })
  rlang::set_names(plots, plot_names)
}

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

scale_within_vector <- function(x) {
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / x_sd)
}

#' Decompose cell-type chromVAR scores into peak contributions
#'
#' @param chromVAR_background_record Cell-type pseudobulk background record.
#' @param psbulk_ATAC_data_matrix Peak-by-cell-type count matrix.
#' @param chromVAR_obj Template chromVAR object with peak ranges.
#' @param annotation_matrix Peak-by-GWAS trait weight matrix.
#' @param GWAS_inputs_tibble GWAS metadata.
#' @return Long tibble whose peak contributions sum to the raw deviation,
#'   relative deviation, and z score for each GWAS and cell type.
#' @keywords internal

get_GWAS_chromVAR_peak_contribution_tibble <- function(
  chromVAR_background_record,
  psbulk_ATAC_data_matrix,
  chromVAR_obj,
  annotation_matrix,
  GWAS_inputs_tibble
) {
  peak_names <- chromVAR_background_record$peak_names
  counts_matrix <- psbulk_ATAC_data_matrix[peak_names, , drop = FALSE]
  if (inherits(counts_matrix, "IterableMatrix")) {
    counts_matrix <- methods::as(counts_matrix, "dgCMatrix")
  }
  annotation_matrix <- annotation_matrix[peak_names, , drop = FALSE]
  background <- chromVAR_background_record$background
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

  GWAS_IDs <- colnames(annotation_matrix)
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

    deviation_vec <- colSums(deviation_contribution_matrix)
    z_vec <- colSums(z_contribution_matrix)
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

#' Summarize peak contribution into the cell-type GWAS heatmap table
#'
#' @param peak_contribution_tibble Exact peak-level chromVAR contributions.
#' @param cell_type_support_tibble Optional nuclei and ATAC-depth support by cell type.
#' @return One row per GWAS and cell type with summed deviations, z-score support,
#'   and optional cell-count annotations.
#' @keywords internal

summarize_GWAS_chromVAR_peak_contributions <- function(
  peak_contribution_tibble,
  cell_type_support_tibble = NULL
) {
  out <- peak_contribution_tibble |>
    dplyr::distinct(
      .data$GWAS_ID,
      .data$Category,
      .data$variant_weighting_mode,
      .data$cluster,
      .data$deviation,
      .data$relative_deviation,
      .data$z
    ) |>
    dplyr::mutate(
      z_p = stats::pnorm(.data$z, lower.tail = FALSE),
      z_q = stats::p.adjust(.data$z_p, method = "BH"),
      support_label = dplyr::case_when(
        .data$z >= 3 ~ "**",
        .data$z >= 2 ~ "*",
        .default = ""
      )
    )

  if (!is.null(cell_type_support_tibble)) {
    out <- out |>
      dplyr::left_join(cell_type_support_tibble, by = "cluster")
  }

  out |>
    dplyr::relocate(
      GWAS_ID,
      Category,
      variant_weighting_mode,
      cluster,
      deviation,
      relative_deviation,
      z,
      z_p,
      z_q,
      support_label,
      dplyr::any_of(c("n_cells", "n_counts", "n_features", "counts_per_feature"))
    )
}

#' Allocate peak contributions to credible-set variants
#'
#' @param peak_contribution_tibble Exact peak-level chromVAR contributions.
#' @param peak_variant_weight_tibble Peak-variant weights.
#' @return Variant-level contribution rows, summed across overlapping peaks.
#' @keywords internal

get_GWAS_chromVAR_variant_contribution_tibble <- function(peak_contribution_tibble, peak_variant_weight_tibble) {
  contribution_cols <- c(
    "deviation_numerator_contribution",
    "deviation_contribution",
    "relative_deviation_contribution",
    "z_contribution"
  )

  peak_variant_contribution_tibble <- peak_contribution_tibble |>
    dplyr::rename(contribution_peak_weight = peak_weight) |>
    dplyr::inner_join(
      peak_variant_weight_tibble |>
        dplyr::select(-dplyr::any_of(c("Category", "variant_weighting_mode"))) |>
        dplyr::rename(mapped_peak_weight = peak_weight),
      by = c("GWAS_ID", "peak_name", "peak_chromosome", "peak_start", "peak_end"),
      relationship = "many-to-many"
    )
  peak_weight_error <- max(
    abs(peak_variant_contribution_tibble$contribution_peak_weight - peak_variant_contribution_tibble$mapped_peak_weight),
    na.rm = TRUE
  )
  if (peak_weight_error > 1e-10) {
    stop("Peak-variant weights do not reproduce the trait-level annotation matrix; maximum error: ", peak_weight_error)
  }

  peak_variant_contribution_tibble |>
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

#' Get Open Targets locus-to-gene predictions
#'
#' Read and rank Open Targets L2G predictions for the credible sets represented
#' in the chromVAR attribution input.
#'
#' @param GWAS_locus_tibble Distinct GWAS and credible-set locus identifiers.
#' @param open_targets_gwas_credible_sets_evidence_dataset_path Local Open
#'   Targets `evidence_gwas_credible_sets` Parquet dataset.
#' @param open_targets_target_dataset_path Local Open Targets `target` Parquet
#'   dataset.
#' @return One row per credible-set gene prediction, ranked by decreasing L2G
#'   score within each GWAS locus.
#' @keywords internal

get_open_targets_GWAS_locus_to_gene_tibble <- function(
  GWAS_locus_tibble,
  open_targets_gwas_credible_sets_evidence_dataset_path,
  open_targets_target_dataset_path
) {
  required_cols <- c("GWAS_ID", "studyId", "studyLocusId", "open_targets_release")
  missing_cols <- setdiff(required_cols, colnames(GWAS_locus_tibble))
  if (length(missing_cols) > 0L) {
    stop("GWAS_locus_tibble is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  selected_loci_tibble <- GWAS_locus_tibble |>
    dplyr::distinct(dplyr::across(dplyr::all_of(required_cols)))
  empty_tibble <- tibble::tibble(
    GWAS_ID = character(),
    studyId = character(),
    studyLocusId = character(),
    targetId = character(),
    gene_symbol = character(),
    gene_name = character(),
    L2G_score = numeric(),
    L2G_rank = integer(),
    open_targets_release = character()
  )
  if (nrow(selected_loci_tibble) == 0L) {
    return(empty_tibble)
  }

  evidence_release <- basename(dirname(open_targets_gwas_credible_sets_evidence_dataset_path))
  target_release <- basename(dirname(open_targets_target_dataset_path))
  selected_releases <- unique(selected_loci_tibble$open_targets_release)
  if (
    length(selected_releases) != 1L ||
      !identical(selected_releases, evidence_release) ||
      !identical(selected_releases, target_release)
  ) {
    stop(
      "Open Targets releases must match across GWAS loci, L2G evidence, and target metadata: ",
      "loci=", paste(selected_releases, collapse = ", "),
      ", evidence=", evidence_release,
      ", target=", target_release
    )
  }

  evidence_tibble <- arrow::open_dataset(open_targets_gwas_credible_sets_evidence_dataset_path) |>
    dplyr::filter(studyLocusId %in% selected_loci_tibble$studyLocusId) |>
    dplyr::select(studyLocusId, targetId, resourceScore, score) |>
    dplyr::collect()
  score_delta <- abs(evidence_tibble$score - evidence_tibble$resourceScore)
  if (any(score_delta > sqrt(.Machine$double.eps), na.rm = TRUE)) {
    stop("Open Targets evidence_gwas_credible_sets score and resourceScore are not equivalent.")
  }

  evidence_tibble <- evidence_tibble |>
    dplyr::filter(!is.na(.data$score), .data$score >= 0.05) |>
    dplyr::transmute(
      studyLocusId = .data$studyLocusId,
      targetId = stringr::str_remove(.data$targetId, "\\.[0-9]+$"),
      L2G_score = .data$score
    ) |>
    dplyr::distinct()
  if (nrow(evidence_tibble) == 0L) {
    return(empty_tibble)
  }

  target_metadata_tibble <- arrow::open_dataset(open_targets_target_dataset_path) |>
    dplyr::filter(id %in% evidence_tibble$targetId) |>
    dplyr::select(id, approvedSymbol, approvedName) |>
    dplyr::collect() |>
    dplyr::transmute(
      targetId = .data$id,
      gene_symbol = dplyr::na_if(stringr::str_squish(.data$approvedSymbol), ""),
      gene_name = dplyr::na_if(stringr::str_squish(.data$approvedName), "")
    )
  missing_target_ids <- setdiff(evidence_tibble$targetId, target_metadata_tibble$targetId)
  if (length(missing_target_ids) > 0L) {
    stop("Open Targets target metadata is missing L2G target IDs: ", paste(missing_target_ids, collapse = ", "))
  }

  evidence_tibble |>
    dplyr::inner_join(selected_loci_tibble, by = "studyLocusId", relationship = "many-to-many") |>
    dplyr::left_join(target_metadata_tibble, by = "targetId", relationship = "many-to-one") |>
    dplyr::mutate(gene_symbol = dplyr::coalesce(.data$gene_symbol, .data$targetId)) |>
    dplyr::arrange(.data$GWAS_ID, .data$studyLocusId, dplyr::desc(.data$L2G_score), .data$gene_symbol, .data$targetId) |>
    dplyr::mutate(L2G_rank = dplyr::row_number(), .by = c(GWAS_ID, studyLocusId)) |>
    dplyr::select(
      GWAS_ID,
      studyId,
      studyLocusId,
      targetId,
      gene_symbol,
      gene_name,
      L2G_score,
      L2G_rank,
      open_targets_release
    )
}

#' Sum credible-set variant contributions by locus
#'
#' @param variant_contribution_tibble Variant-level chromVAR contribution.
#' @param locus_to_gene_tibble Ranked Open Targets L2G predictions.
#' @return Exact locus-level contributions with lead-variant labels.
#' @keywords internal

get_GWAS_chromVAR_locus_contribution_tibble <- function(variant_contribution_tibble, locus_to_gene_tibble) {
  lead_variant_tibble <- variant_contribution_tibble |>
    dplyr::arrange(dplyr::desc(.data$posteriorProbability), .data$variantId) |>
    dplyr::slice_head(n = 1, by = c(GWAS_ID, studyLocusId)) |>
    dplyr::select(
      GWAS_ID,
      studyLocusId,
      lead_variantId = variantId,
      lead_variant_position = position,
      lead_variant_PIP = posteriorProbability
    )

  locus_gene_summary_tibble <- locus_to_gene_tibble |>
    dplyr::summarise(
      top_L2G_gene = dplyr::first(.data$gene_symbol),
      top_L2G_score = dplyr::first(.data$L2G_score),
      likely_genes_label = paste0(
        .data$gene_symbol[.data$L2G_rank <= 3L],
        " (",
        sprintf("%.2f", .data$L2G_score[.data$L2G_rank <= 3L]),
        ")",
        collapse = ", "
      ),
      n_L2G_genes = dplyr::n(),
      .by = c(GWAS_ID, studyLocusId)
    )

  variant_contribution_tibble |>
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
    dplyr::left_join(locus_gene_summary_tibble, by = c("GWAS_ID", "studyLocusId"), relationship = "many-to-one") |>
    dplyr::mutate(
      locus_label = .data$lead_variantId,
      likely_genes_label = dplyr::coalesce(.data$likely_genes_label, "No L2G prediction >= 0.05"),
      n_L2G_genes = dplyr::coalesce(.data$n_L2G_genes, 0L),
      contribution_rank = dplyr::min_rank(dplyr::desc(abs(.data$relative_deviation_contribution))),
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

#' Validate additive GWAS chromVAR contribution
#'
#' @param peak_contribution_tibble Peak-level contribution.
#' @param variant_contribution_tibble Variant-level contribution.
#' @param locus_contribution_tibble Locus-level contribution.
#' @return Reconciliation errors for each GWAS and cell type.
#' @keywords internal

get_GWAS_chromVAR_contribution_reconciliation_tibble <- function(
  peak_contribution_tibble,
  variant_contribution_tibble,
  locus_contribution_tibble
) {
  summarize_level <- function(contribution_tibble, level) {
    contribution_tibble |>
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
    summarize_level(peak_contribution_tibble, "peak"),
    summarize_level(variant_contribution_tibble, "variant"),
    summarize_level(locus_contribution_tibble, "locus")
  ) |>
    dplyr::mutate(
      deviation_error = .data$deviation_contribution_sum - .data$deviation,
      relative_deviation_error = .data$relative_deviation_contribution_sum - .data$relative_deviation,
      z_error = .data$z_contribution_sum - .data$z,
      contribution_reconciles = abs(.data$deviation_error) < 1e-8 &
        abs(.data$relative_deviation_error) < 1e-8 &
        abs(.data$z_error) < 1e-8
    ) |>
    dplyr::relocate(GWAS_ID, cluster, level)
}

collapse_GWAS_locus_contribution_for_plot <- function(locus_contribution_tibble, n_top_loci = 15L) {
  top_locus_ids <- locus_contribution_tibble |>
    dplyr::summarise(
      max_abs_contribution = max(abs(.data$relative_deviation_contribution)),
      locus_label = dplyr::first(.data$locus_label),
      top_L2G_gene = dplyr::first(.data$top_L2G_gene),
      .by = studyLocusId
    ) |>
    dplyr::slice_max(.data$max_abs_contribution, n = n_top_loci, with_ties = FALSE) |>
    dplyr::arrange(dplyr::desc(.data$max_abs_contribution)) |>
    dplyr::mutate(
      plot_locus = dplyr::if_else(
        is.na(.data$top_L2G_gene),
        .data$locus_label,
        stringr::str_c(.data$locus_label, .data$top_L2G_gene, sep = "\n")
      ),
      plot_locus = dplyr::if_else(
        duplicated(.data$plot_locus) | duplicated(.data$plot_locus, fromLast = TRUE),
        stringr::str_c(.data$plot_locus, .data$studyLocusId, sep = "\n"),
        .data$plot_locus
      )
    )

  top_tibble <- locus_contribution_tibble |>
    dplyr::semi_join(top_locus_ids, by = "studyLocusId") |>
    dplyr::left_join(
      top_locus_ids |>
        dplyr::select(studyLocusId, plot_locus),
      by = "studyLocusId",
      relationship = "many-to-one"
    )
  other_tibble <- locus_contribution_tibble |>
    dplyr::anti_join(top_locus_ids, by = "studyLocusId") |>
    dplyr::mutate(plot_locus = dplyr::if_else(.data$relative_deviation_contribution >= 0, "Other positive", "Other negative")) |>
    dplyr::summarise(
      relative_deviation_contribution = sum(.data$relative_deviation_contribution),
      .by = c(cluster, plot_locus)
    )
  locus_order <- c(top_locus_ids$plot_locus, "Other positive", "Other negative")

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
#' @param locus_contribution_tibble Exact locus-level contribution.
#' @param chromVAR_deviation_tibble Original cell-type-by-GWAS heatmap data.
#' @param n_top_loci Number of individually labelled loci per GWAS.
#' @return Named list with one contribution heatmap per GWAS.
#' @keywords internal

plot_GWAS_locus_contribution_heatmaps <- function(
  locus_contribution_tibble,
  chromVAR_deviation_tibble,
  n_top_loci = 15L
) {
  GWAS_IDs <- unique(locus_contribution_tibble$GWAS_ID)
  plots <- purrr::map(GWAS_IDs, \(GWAS_ID) {
    locus_tibble <- locus_contribution_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID)
    heatmap_tibble <- collapse_GWAS_locus_contribution_for_plot(locus_tibble, n_top_loci = n_top_loci)
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
      patchwork::plot_annotation(title = paste("Locus contributions to GWAS-linked accessibility:", GWAS_ID),
        subtitle = stringr::str_wrap("Look for enrichment dominated by a few loci versus distributed support; the adjacent bar shows the total.", width = 100),
        caption = stringr::str_wrap(paste("Loci are ranked by their largest absolute relative-deviation contribution across cell types; up to", n_top_loci,
          "are shown individually. Remaining contributions are summed as Other positive/negative. Contributions sum to the total deviation standardized across cell types. Total-bar stars mark background Z >= 2 (*) or >= 3 (**), not FDR. Open Targets gene labels are prioritizations, not causal assignments."), width = 110))
    combined_plot$labels$title <- ggplot2::waiver()
    combined_plot
  })
  names(plots) <- GWAS_IDs
  plots
}

prepare_GWAS_locus_contribution_waterfall_tibble <- function(locus_tibble, n_top_loci = 15L) {
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
      dplyr::select(
        studyLocusId,
        locus_label,
        top_L2G_gene,
        top_L2G_score,
        likely_genes_label,
        relative_deviation_contribution
      ),
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
      direction = dplyr::if_else(.data$relative_deviation_contribution >= 0, "Positive", "Negative"),
      gene_label = dplyr::if_else(
        is.na(.data$top_L2G_gene),
        "",
        stringr::str_glue("{.data$top_L2G_gene}\nL2G {sprintf('%.2f', .data$top_L2G_score)}")
      )
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
      direction = "Total",
      gene_label = ""
    )
  )
}

#' Plot one locus-contribution waterfall
#'
#' @param waterfall_tibble Minimal plotting data returned by
#'   `prepare_GWAS_locus_contribution_waterfall_tibble()`.
#' @param title Plot title identifying the GWAS and cell type.
#' @return A ggplot ready for saving or composition.
#' @keywords internal

plot_GWAS_locus_contribution_waterfall <- function(waterfall_tibble, title) {
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
    ggrepel::geom_label_repel(
      data = dplyr::filter(waterfall_tibble, .data$gene_label != ""),
      ggplot2::aes(x = .data$plot_index, y = .data$end, label = .data$gene_label),
      direction = "y",
      seed = 1,
      size = 2.2,
      min.segment.length = 0,
      box.padding = 0.2,
      point.padding = 0.1,
      max.overlaps = Inf,
      fill = scales::alpha("white", 0.9),
      color = "grey15",
      linewidth = 0.15
    ) +
    ggplot2::scale_x_continuous(
      breaks = waterfall_tibble$plot_index,
      labels = waterfall_tibble$locus_label,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::scale_fill_manual(values = c(Positive = "#B40426", Negative = "#3B4CC0", Total = "grey25"), name = NULL) +
    ggplot2::labs(
      x = "Credible-set locus (lead variant)",
      y = "Cumulative relative deviation",
      title = paste("Locus contributions:", title),
      subtitle = stringr::str_wrap("Follow positive and negative steps to see which loci drive or offset the total; the ordering is not genomic.", width = 100),
      caption = stringr::str_wrap("Step height is a locus's relative-deviation contribution; the final bar is their sum. Selected loci are shown individually; others are pooled by sign. Labels show the top Open Targets L2G gene and score, not a proven causal gene. These are descriptive contributions to a pooled cell-type score.", width = 110)
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

#' Plot locus-contribution waterfalls
#'
#' @param locus_contribution_tibble Exact locus-level contribution.
#' @param n_top_loci Number of individually labelled loci per plot.
#' @return Named list with one waterfall per GWAS and cell type.
#' @keywords internal

plot_GWAS_locus_contribution_waterfalls <- function(locus_contribution_tibble, n_top_loci = 15L) {
  split_tibbles <- split(
    locus_contribution_tibble,
    interaction(locus_contribution_tibble$GWAS_ID, locus_contribution_tibble$cluster, drop = TRUE, lex.order = TRUE)
  )

  plot_names <- if (dplyr::n_distinct(locus_contribution_tibble$GWAS_ID) == 1L) {
    vapply(split_tibbles, \(x) as.character(x$cluster[[1]]), character(1))
  } else names(split_tibbles)

  purrr::map(split_tibbles, \(locus_tibble) {
    plot_GWAS_locus_contribution_waterfall(
      waterfall_tibble = prepare_GWAS_locus_contribution_waterfall_tibble(
        locus_tibble,
        n_top_loci = n_top_loci
      ),
      title = stringr::str_glue("{locus_tibble$GWAS_ID[[1]]} - {locus_tibble$cluster[[1]]}")
    )
  }) |>
    rlang::set_names(stringr::str_replace_all(plot_names, "[^A-Za-z0-9_.-]+", "_"))
}

#' Prepare variant-level contribution detail records for top loci
#'
#' @param variant_contribution_tibble Variant-level contribution.
#' @param locus_contribution_tibble Locus-level contribution used to select loci.
#' @param absolute_effect_locus_tibble Loci ranked under absolute-effect weighting.
#' @param consensus_peak_GRanges Consensus ATAC peaks.
#' @param fragments BPCells fragment object.
#' @param metadata_tibble Cell metadata with cell-type labels.
#' @param n_top_cell_types Number of highest-deviation cell types detailed per GWAS.
#' @param n_top_loci Number of loci plotted per GWAS and cell type.
#' @param flank Bases added around the selected locus.
#' @return Named list of plot-ready variant, accessibility, and peak-track records.
#' @keywords internal

prepare_GWAS_variant_contribution_detail_records <- function(
  variant_contribution_tibble,
  locus_contribution_tibble,
  absolute_effect_locus_tibble,
  consensus_peak_GRanges,
  fragments,
  metadata_tibble,
  n_top_cell_types = 1L,
  n_top_loci = 3L,
  flank = 25000L,
  group_cells_by_col = "PCA_harmony_SNN_cluster_cell_type"
) {
  selected_locus_tibble <- select_GWAS_detail_loci(locus_contribution_tibble,
    absolute_effect_locus_tibble, n_top_cell_types, n_top_loci)

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
  plot_records <- locus_tibble_list |>
    purrr::map(\(locus_tibble) {
      GWAS_ID <- locus_tibble$GWAS_ID[[1]]
      cluster <- locus_tibble$cluster[[1]]
      studyLocusId <- locus_tibble$studyLocusId[[1]]
      variant_tibble <- variant_contribution_tibble |>
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
      coverage_colors <- stats::setNames(
        c("#B40426", rep("grey75", length(group_levels) - 1L)),
        group_levels
      )
      coverage_tibble <- BPCells::trackplot_coverage(
        fragments = fragments,
        region = region,
        groups = factor(group_values, levels = group_levels),
        cell_read_counts = cell_read_counts,
        group_order = group_levels,
        colors = coverage_colors,
        bins = 500,
        return_data = TRUE
      )

      plot_title <- stringr::str_glue(
        "{GWAS_ID} - {cluster} - {locus_tibble$locus_label[[1]]}; ",
        "locus contribution={round(locus_tibble$relative_deviation_contribution[[1]], 3)}"
      )
      list(
        GWAS_ID = GWAS_ID,
        selection = locus_tibble$selection[[1]],
        variant_tibble = variant_tibble,
        coverage_tibble = coverage_tibble,
        coverage_colors = coverage_colors,
        consensus_peak_GRanges = IRanges::subsetByOverlaps(consensus_peak_GRanges, region),
        region = region,
        plot_title = plot_title
      )
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
  rlang::set_names(plot_records, plot_names)
}

#' Plot aligned variant-detail tracks for one GWAS and cell type
#'
#' @param plot_records Cached locus coverage and contribution records.
#' @param variants Effect metadata from the original credible sets.
#' @param gene_GRanges Reference gene bodies.
#' @return A patchwork with loci as columns and genomic tracks as rows.
plot_GWAS_variant_contribution_detail_panels <- function(plot_records, variants, gene_GRanges) {
  plot_records <- plot_records[order(vapply(plot_records, \(record)
    abs(sum(record$variant_tibble$relative_deviation_contribution)), numeric(1)), decreasing = TRUE)]
  loci <- purrr::imap_chr(unname(plot_records), \(record, index) {
    variant <- record$variant_tibble
    lead_id <- variants$lead_id[match(variant$studyLocusId[[1]], variants$studyLocusId)]
    gene <- variants$top_L2G_gene[match(variant$studyLocusId[[1]], variants$studyLocusId)]
    paste0("Locus ", index, ": ", lead_id,
      if (!is.na(gene) && nzchar(gene)) paste0(" — L2G: ", gene) else "",
      "\nContribution = ", signif(sum(variant$relative_deviation_contribution), 3),
      if (!is.null(record$selection)) paste0("; selected: ", record$selection) else "")
  })
  bounds <- purrr::map2_dfr(plot_records, loci, \(record, locus) tibble::tibble(
    locus = locus, position = c(GenomicRanges::start(record$region), GenomicRanges::end(record$region))))
  contributions <- purrr::map2_dfr(plot_records, loci, \(record, locus) {
    record$variant_tibble |> dplyr::mutate(locus = locus)
  }) |>
    dplyr::left_join(variants, by = c("studyLocusId", "variantId"), relationship = "many-to-one")
  coverage <- purrr::map2_dfr(plot_records, loci, \(record, locus) {
    ymax <- as.numeric(stats::quantile(record$coverage_tibble$normalized_insertions, 0.999))
    record$coverage_tibble |> dplyr::mutate(locus = locus,
      coverage = if (ymax > 0) pmin(normalized_insertions, ymax) / ymax else 0,
      range_label = paste0("0–", signif(ymax, 3)))
  })
  genes <- purrr::map2_dfr(plot_records, loci, \(record, locus) {
    region <- record$region
    selected <- IRanges::subsetByOverlaps(gene_GRanges, region)
    tibble::tibble(locus = locus, start = pmax(GenomicRanges::start(selected), GenomicRanges::start(region)),
      end = pmin(GenomicRanges::end(selected), GenomicRanges::end(region)),
      gene = as.character(selected$gene_name), strand = as.character(GenomicRanges::strand(selected))) |>
      dplyr::arrange(start, end) |>
      dplyr::mutate(lane = IRanges::disjointBins(IRanges::IRanges(start, end)))
  })
  peaks <- purrr::map2_dfr(plot_records, loci, \(record, locus) {
    ranges <- record$consensus_peak_GRanges
    tibble::tibble(locus = locus,
      start = pmax(GenomicRanges::start(ranges), GenomicRanges::start(record$region)),
      end = pmin(GenomicRanges::end(ranges), GenomicRanges::end(record$region)))
  })
  bounds$locus <- factor(bounds$locus, levels = loci)
  contributions$locus <- factor(contributions$locus, levels = loci)
  coverage$locus <- factor(coverage$locus, levels = loci)
  genes$locus <- factor(genes$locus, levels = loci)
  peaks$locus <- factor(peaks$locus, levels = loci)
  facet <- ggplot2::facet_grid(cols = ggplot2::vars(locus), scales = "free_x")
  track_theme <- ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, colour = "grey65"),
    panel.spacing.x = grid::unit(1, "lines"),
    strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey65"),
    strip.text = ggplot2::element_text(face = "bold"), legend.position = "top")
  # Training every track on identical bounds keeps genomic coordinates aligned.
  base_track <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = bounds, ggplot2::aes(position, 0)) + facet + track_theme +
    ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::label_number())
  effect_max <- max(c(0, contributions$effect_magnitude), na.rm = TRUE)
  contribution_track <- base_track +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_segment(data = contributions, ggplot2::aes(x = position, xend = position,
      y = 0, yend = relative_deviation_contribution, colour = effect_magnitude), linewidth = 0.4) +
    ggplot2::geom_point(data = contributions, ggplot2::aes(position, relative_deviation_contribution,
      size = posteriorProbability_raw, colour = effect_magnitude, shape = effect_source), stroke = 0.8, show.legend = TRUE) +
    ggrepel::geom_text_repel(data = contributions |>
      dplyr::slice_max(abs(relative_deviation_contribution), n = 3, with_ties = FALSE, by = locus),
      ggplot2::aes(position, relative_deviation_contribution, label = variantId),
      size = 2.3, seed = 1, min.segment.length = 0, max.overlaps = Inf) +
    ggplot2::scale_size_area(max_size = 4, limits = c(0, 1), breaks = c(0.1, 0.5, 1), name = "PIP") +
    ggplot2::scale_colour_viridis_c(option = "C", limits = c(0, if (effect_max > 0) effect_max else 1),
      na.value = "grey60", breaks = scales::breaks_pretty(n = 3), name = "Effect magnitude |β|") +
    ggplot2::scale_shape_manual(values = c(Variant = 16, `Lead variant proxy` = 1, Unavailable = 16),
      limits = c("Variant", "Lead variant proxy", "Unavailable"), drop = FALSE, name = "Effect source") +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = c("black", "black", "grey60")))) +
    ggplot2::labs(x = NULL, y = "Relative\ndeviation\ncontribution") +
    ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5))
  gene_track <- base_track +
    ggplot2::geom_segment(data = genes, ggplot2::aes(
      x = ifelse(strand == "-", end, start), xend = ifelse(strand == "-", start, end),
      y = -lane, yend = -lane), arrow = grid::arrow(length = grid::unit(1.3, "mm")), linewidth = 0.5) +
    ggplot2::geom_text(data = genes, ggplot2::aes((start + end) / 2, -lane, label = gene),
      vjust = -0.7, size = 2.5, check_overlap = TRUE) +
    ggplot2::geom_text(data = dplyr::filter(bounds, !locus %in% genes$locus) |>
      dplyr::summarise(position = mean(position), .by = locus),
      ggplot2::aes(position, 0, label = "No annotated genes overlap"), size = 3) +
    ggplot2::labs(x = NULL, y = "Gene bodies") + ggplot2::guides(y = "none") +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  coverage_track <- base_track +
    ggplot2::geom_area(data = coverage, ggplot2::aes(pos, coverage, fill = group)) +
    ggplot2::geom_text(data = dplyr::distinct(coverage, locus, group, range_label),
      ggplot2::aes(x = -Inf, y = 1, label = range_label), hjust = -0.1, vjust = 1.1, size = 2.2) +
    ggplot2::facet_grid(rows = ggplot2::vars(group), cols = ggplot2::vars(locus), scales = "free_x") +
    ggplot2::scale_fill_manual(values = plot_records[[1]]$coverage_colors, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "ATAC coverage (RPKM; range shown)") + ggplot2::guides(y = "none") +
    ggplot2::theme(strip.text.y.right = ggplot2::element_text(angle = 0, size = 8),
      panel.grid = ggplot2::element_blank(), panel.spacing.y = grid::unit(0.1, "lines"))
  peak_track <- base_track +
    ggplot2::geom_segment(data = peaks, ggplot2::aes(start, 0, xend = end, yend = 0), linewidth = 2) +
    ggplot2::labs(x = "Genomic position (bp)", y = "Consensus peaks") + ggplot2::guides(y = "none") +
    ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, size = 9),
      panel.grid = ggplot2::element_blank())
  tracks <- list(contribution_track, gene_track, coverage_track, peak_track)
  tracks[1:3] <- purrr::map(tracks[1:3], \(plot) plot + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()))
  tracks[2:4] <- purrr::map(tracks[2:4], \(plot) plot + ggplot2::theme(strip.text.x = ggplot2::element_blank()))
  plot <- patchwork::wrap_plots(tracks, ncol = 1,
    heights = grid::unit.c(grid::unit(c(1, max(1, 0.3 * max(c(0, genes$lane))),
      0.45 * nlevels(coverage$group)), "null"), grid::unit(4, "mm")), guides = "collect") +
    patchwork::plot_annotation(
      title = paste("Variant contributions:", plot_records[[1]]$GWAS_ID, "—", plot_records[[1]]$variant_tibble$cluster[[1]]),
      subtitle = "Compare loci on a shared contribution scale; point area shows PIP and colour shows effect magnitude.\nGene bodies and accessibility provide context, not proof of a causal variant or target gene.",
      caption = stringr::str_wrap(paste(
        "Loci are the union of the leading ordinary and absolute-effect-weighted contributions for the focal cell type; plotted stems remain ordinary contributions. PIP is the original fine-mapping probability, before effect weighting. Filled points use the variant's own |β|; hollow points inherit the locus lead variant's |β|; grey means unavailable.",
        "The lead is the recorded lead variant, or the highest-PIP variant if no lead is recorded. No signed effects are inherited. Effect units are study-specific.",
        "Facet genes are the top available Open Targets L2G predictions, not established causal genes. Gene arrows show strand; gene bodies are clipped to each window and overlapping labels may be omitted. Coverage uses cached depth-normalized 500-bin profiles, clipped at each locus's 99.9th percentile; displayed coverage ranges differ between loci."), 180),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))) &
    ggplot2::theme(legend.position = "top")
  # The composite already has a title; prevent list saving from titling its last track.
  plot$labels$title <- ggplot2::waiver()
  plot
}

#' Plot one aligned locus figure per GWAS and focal cell type
#' @param plot_records Cached locus records for one GWAS.
#' @param GWAS_input_records Original credible sets, including effect metadata.
#' @param gene_GRanges Reference gene-body annotations.
#' @param locus_contribution_tibble Locus contributions with top L2G predictions.
#' @return Named list of patchworks, one per focal cell type.
plot_GWAS_variant_contribution_details <- function(plot_records, GWAS_input_records, gene_GRanges, locus_contribution_tibble) {
  GWAS_ID <- plot_records[[1]]$GWAS_ID
  input <- purrr::keep(GWAS_input_records, \(record) identical(record$GWAS_ID, GWAS_ID))[[1]]
  variants <- tibble::as_tibble(as.data.frame(S4Vectors::mcols(input$credible_set_GRanges))) |>
    dplyr::mutate(beta = dplyr::if_else(is.finite(beta), beta, NA_real_)) |>
    dplyr::group_by(studyLocusId) |>
    dplyr::arrange(dplyr::desc(posteriorProbability_raw), variantId, .by_group = TRUE) |>
    dplyr::mutate(
      lead_id = dplyr::coalesce(dplyr::first(lead_variantId), dplyr::first(variantId)),
      lead_beta = beta[match(lead_id, variantId)],
      effect_magnitude = abs(dplyr::coalesce(beta, lead_beta)),
      effect_source = dplyr::case_when(!is.na(beta) ~ "Variant", !is.na(lead_beta) ~ "Lead variant proxy", TRUE ~ "Unavailable")) |>
    dplyr::ungroup() |>
    dplyr::select(studyLocusId, variantId, lead_id, effect_magnitude, effect_source) |>
    dplyr::left_join(locus_contribution_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID) |>
      dplyr::distinct(studyLocusId, top_L2G_gene), by = "studyLocusId", relationship = "many-to-one")
  groups <- split(plot_records, purrr::map_chr(plot_records, \(record) record$variant_tibble$cluster[[1]]))
  purrr::map(groups, plot_GWAS_variant_contribution_detail_panels, variants = variants, gene_GRanges = gene_GRanges) |>
    rlang::set_names(stringr::str_replace_all(names(groups), "[^A-Za-z0-9_.-]+", "_"))
}

#' Plot compact locus contributions across cell types for one GWAS
#'
#' @param locus_contribution_tibble Exact contributions for one GWAS.
#' @param n_top_loci Number of individually labelled loci per cell type.
#' @return A faceted ggplot with a shared contribution scale.
plot_GWAS_locus_contribution_bars <- function(locus_contribution_tibble, n_top_loci = 5L) {
  plot_data <- locus_contribution_tibble |>
    dplyr::group_by(cluster) |>
    dplyr::group_modify(\(data, key) {
      steps <- prepare_GWAS_locus_contribution_waterfall_tibble(data, n_top_loci)
      total <- steps$relative_deviation_contribution[steps$locus_label == "Total"]
      bars <- steps |> dplyr::filter(locus_label != "Total")
      if (nrow(bars) == 0L) bars <- tibble::tibble(
        locus_label = "No nonzero contributions", relative_deviation_contribution = 0,
        top_L2G_gene = NA_character_, is_other = TRUE, direction = "Positive")
      bars |>
        dplyr::arrange(is_other, dplyr::desc(abs(relative_deviation_contribution))) |>
        dplyr::mutate(
          panel = paste0(key$cluster, " | Total = ", signif(total, 3)),
          label = dplyr::if_else(!is.na(top_L2G_gene) & nzchar(top_L2G_gene),
            paste0(top_L2G_gene, " · ", locus_label), locus_label))
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(row = factor(dplyr::row_number(), levels = rev(seq_len(dplyr::n()))),
      panel = factor(panel, levels = unique(panel)))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = relative_deviation_contribution, y = row, fill = direction)) +
    ggplot2::geom_blank() +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_col(data = \(data) dplyr::filter(data, !is_other), width = 0.7) +
    ggplot2::geom_col(data = \(data) dplyr::filter(data, is_other),
      ggplot2::aes(colour = direction), fill = NA, width = 0.7, linewidth = 0.65, show.legend = FALSE) +
    ggplot2::scale_y_discrete(labels = stats::setNames(plot_data$label, plot_data$row)) +
    ggplot2::scale_fill_manual(values = c(Positive = "#B40426", Negative = "#3B4CC0"), name = NULL) +
    ggplot2::scale_colour_manual(values = c(Positive = "#B40426", Negative = "#3B4CC0"), guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(panel), ncol = 3, scales = "free_y") +
    ggplot2::labs(title = paste("Loci contributing to GWAS-linked accessibility:", unique(locus_contribution_tibble$GWAS_ID)),
      subtitle = "Compare bar lengths on the shared scale: positive loci drive the total; negative loci offset it.\nEach cell type has its own leading loci; separate remainder bars preserve cancellation.",
      x = "Relative-deviation contribution", y = NULL,
      caption = stringr::str_wrap(paste("Up to", n_top_loci,
        "loci per cell type are selected by absolute contribution; remaining loci are summed separately by sign and drawn as unfilled outlines. Bars sum to the total in each panel heading. Zero contributions are omitted except in all-zero panels. Gene labels are top Open Targets L2G assignments, not proven causal genes. Contributions describe pooled cell-type chromVAR relative deviations, not significance or cumulative sums."), 150)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom",
      panel.border = ggplot2::element_rect(colour = "grey55", fill = NA, linewidth = 0.5),
      panel.spacing = grid::unit(1.2, "lines"),
      strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey55", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = 11,
        margin = ggplot2::margin(7, 5, 7, 5)), plot.caption = ggplot2::element_text(hjust = 0))
}

#' Select the union of leading loci under ordinary and absolute-effect weighting
#' @param locus_contribution_tibble Ordinary locus contributions for one or more GWAS.
#' @param absolute_effect_locus_tibble Absolute-effect contributions, possibly empty.
#' @param n_top_cell_types Number of focal cell types selected by ordinary deviation.
#' @param n_top_loci Number of loci retained per ranking and focal cell type.
#' @return Unique loci with selection provenance and coordinates from the ordinary analysis.
select_GWAS_detail_loci <- function(locus_contribution_tibble, absolute_effect_locus_tibble,
  n_top_cell_types = 1L, n_top_loci = 3L) {
  focal <- locus_contribution_tibble |>
    dplyr::distinct(GWAS_ID, cluster, relative_deviation) |>
    dplyr::slice_max(relative_deviation, n = n_top_cell_types, with_ties = FALSE, by = GWAS_ID)
  leading <- function(data, ranking) {
    data |> dplyr::semi_join(focal, by = c("GWAS_ID", "cluster")) |>
      dplyr::slice_max(abs(relative_deviation_contribution), n = n_top_loci,
        with_ties = FALSE, by = c(GWAS_ID, cluster)) |>
      dplyr::transmute(GWAS_ID, cluster, studyLocusId, selection = ranking)
  }
  selected <- dplyr::bind_rows(leading(locus_contribution_tibble, "Ordinary"),
    leading(absolute_effect_locus_tibble, "Absolute effect")) |>
    dplyr::summarise(selection = paste(selection, collapse = " + "), .by = c(GWAS_ID, cluster, studyLocusId))
  locus_contribution_tibble |>
    dplyr::inner_join(selected, by = c("GWAS_ID", "cluster", "studyLocusId"), relationship = "one-to-one") |>
    dplyr::mutate(detail_rank = dplyr::row_number(), .by = c(GWAS_ID, cluster))
}

#' Plot absolute-effect locus contributions using the ordinary faceted layout
#' @param locus_contribution_tibble Exact absolute-effect-weighted contributions.
#' @return One named plot per eligible GWAS, or an empty plot list.
plot_GWAS_absolute_effect_locus_bars <- function(locus_contribution_tibble) {
  if (nrow(locus_contribution_tibble) == 0L) return(structure(list(), class = c("empty_plot_list", "list")))
  split(locus_contribution_tibble, locus_contribution_tibble$GWAS_ID) |>
    purrr::map(\(data) plot_GWAS_locus_contribution_bars(data, n_top_loci = 5L) +
      ggplot2::labs(title = paste("Effect-weighted loci contributing to GWAS-linked accessibility:", data$GWAS_ID[[1]]),
        subtitle = "Compare leading loci after PIP × |β| weighting; bar lengths share a scale across cell types.\nPositive loci drive the effect-weighted total and negative loci offset it; genetic effect direction is discarded.",
        caption = stringr::str_wrap(paste(
          "Weights match the absolute-effect heatmap: variant |β| when all variants have effects, otherwise the highest-PIP available |β| per locus. Only eligible GWAS are included.",
          "Weights preserve retained PIP mass before peak weights are capped at one; contributions are recomputed and sum to the effect-weighted relative deviation.",
          "Five loci per cell type are selected by absolute contribution. Unfilled remainder bars pool other positive and negative loci separately. Gene labels are top Open Targets L2G predictions, not proven causal genes."), 150)))
}

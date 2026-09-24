#' Get BPCells pseudobulk matrix
#'
#' Sum a feature-by-cell matrix into cluster-by-donor pseudobulk samples.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param threads Number of threads passed to BPCells, HNSW, or matrix-stat routines.
#' @return A named matrix-like object with rows and columns aligned to the input feature/cell identifiers.
#' @keywords internal

get_BPCells_pseudobulk_matrix <- function(feature_matrix, metadata_tibble, cluster_col, threads = 1) {
  cell_metadata <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::transmute(
      barcode_w_prefix,
      cluster = stringr::str_replace_all(as.character(.data[[cluster_col]]), "_", "-"),
      donor_id = stringr::str_replace_all(as.character(donor_id), "_", "-")
    ) |>
    dplyr::filter(.data$barcode_w_prefix %in% colnames(feature_matrix)) |>
    dplyr::arrange(match(.data$barcode_w_prefix, colnames(feature_matrix)))

  aligned_feature_matrix <- feature_matrix[, cell_metadata$barcode_w_prefix, drop = FALSE]
  cell_groups <- stringr::str_c(cell_metadata$cluster, "_", cell_metadata$donor_id)

  BPCells::pseudobulk_matrix(
    mat = aligned_feature_matrix,
    cell_groups = cell_groups,
    method = "sum",
    threads = threads
  )
}

get_BPCells_group_pseudobulk_matrix <- function(feature_matrix, metadata_tibble, group_col, threads = 1) {
  cell_metadata <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::transmute(
      barcode_w_prefix,
      group = stringr::str_replace_all(as.character(.data[[group_col]]), "_", "-")
    ) |>
    dplyr::filter(.data$barcode_w_prefix %in% colnames(feature_matrix)) |>
    dplyr::arrange(match(.data$barcode_w_prefix, colnames(feature_matrix)))

  aligned_feature_matrix <- feature_matrix[, cell_metadata$barcode_w_prefix, drop = FALSE]

  BPCells::pseudobulk_matrix(
    mat = aligned_feature_matrix,
    cell_groups = cell_metadata$group,
    method = "sum",
    threads = threads
  )
}

pseudobulk_rowSums <- function(x) {
  if (inherits(x, "IterableMatrix")) {
    BPCells::rowSums(x)
  } else {
    Matrix::rowSums(x)
  }
}

pseudobulk_colSums <- function(x) {
  if (inherits(x, "IterableMatrix")) {
    BPCells::colSums(x)
  } else {
    Matrix::colSums(x)
  }
}

get_pseudobulk_depth_tibble <- function(pseudobulk_data_matrix, count_col = "n_counts") {
  n_counts <- as.numeric(pseudobulk_colSums(pseudobulk_data_matrix))
  n_features <- as.numeric(pseudobulk_colSums(pseudobulk_data_matrix > 0))

  get_pseudobulk_sample_tibble(pseudobulk_data_matrix) |>
    dplyr::mutate(
      !!count_col := n_counts[match(sample_name, colnames(pseudobulk_data_matrix))],
      n_features = n_features[match(sample_name, colnames(pseudobulk_data_matrix))],
      counts_per_feature = .data[[count_col]] / n_features,
      log10_n_counts = log10(.data[[count_col]]),
      log10_n_features = log10(n_features),
      log10_counts_per_feature = log10(counts_per_feature)
    )
}

get_group_pseudobulk_depth_tibble <- function(pseudobulk_data_matrix, group_col = "cluster", count_col = "n_counts") {
  n_counts <- as.numeric(pseudobulk_colSums(pseudobulk_data_matrix))
  n_features <- as.numeric(pseudobulk_colSums(pseudobulk_data_matrix > 0))

  tibble::tibble(
    !!group_col := colnames(pseudobulk_data_matrix),
    !!count_col := n_counts,
    n_features = n_features,
    counts_per_feature = n_counts / n_features,
    log10_n_counts = log10(n_counts),
    log10_n_features = log10(n_features),
    log10_counts_per_feature = log10(n_counts / n_features)
  )
}

get_group_cell_count_tibble <- function(metadata_tibble, group_col, output_col = "cluster") {
  metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::transmute(!!output_col := stringr::str_replace_all(as.character(.data[[group_col]]), "_", "-")) |>
    dplyr::count(.data[[output_col]], name = "n_cells")
}

#' Plot pseudobulk depth distribution
#'
#' Plot pseudobulk sample depth, detected features, and counts per feature.
#'
#' @param pseudobulk_depth_tibble Depth summary tibble with sample, cluster, and
#'   `n_counts`, `n_features`, and `counts_per_feature` metrics.
#' @param min_ATAC_sample_counts Optional vertical threshold shown on count
#'   panels for ATAC pseudobulk sample filtering.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_pseudobulk_depth_distribution <- function(pseudobulk_depth_tibble, min_ATAC_sample_counts = NULL) {
  plot_tibble <- pseudobulk_depth_tibble |>
    tidyr::pivot_longer(
      cols = c(n_counts, n_features, counts_per_feature),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c("n_counts", "n_features", "counts_per_feature"),
        labels = c("Total counts", "Detected features", "Counts per detected feature")
      )
    )

  plot <- plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(cluster, value, median), y = value, color = cluster)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 1) +
    ggplot2::scale_y_log10(labels = scales::label_number()) +
    ggplot2::facet_grid(metric ~ modality, scales = "free_y") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      x = "Cluster",
      y = NULL,
      title = "Pseudobulk depth per cluster-donor sample",
      subtitle = stringr::str_wrap("Look for sparse donor-group samples before fitting models; pooled depth can conceal uneven donor support.", width = 100),
      caption = stringr::str_wrap("Points represent cluster-donor samples; boxes show medians and interquartile ranges, with whiskers to 1.5 x IQR. Rows use separate logarithmic depth scales. Detected features have positive counts; counts per feature use only detected features.", width = 110)
    )

  if (!is.null(min_ATAC_sample_counts)) {
    plot <- plot +
      ggplot2::geom_hline(
        data = tibble::tibble(
          modality = "ATAC",
          metric = factor("Total counts", levels = levels(plot_tibble$metric)),
          value = min_ATAC_sample_counts
        ),
        ggplot2::aes(yintercept = value),
        linetype = "dashed",
        inherit.aes = FALSE
      ) +
      ggplot2::labs(caption = stringr::str_wrap(paste(plot$labels$caption,
        stringr::str_glue("Dashed line marks motif_family_accessibility minimum ATAC count threshold: {min_ATAC_sample_counts}.")), width = 110))
  }

  plot
}

normalize_pseudobulk_feature_models <- function(cfg_pseudobulk_feature_models) {
  normalize_differential_models(cfg_pseudobulk_feature_models) |>
    purrr::map(function(model) {
      if (!is.null(model$GEM_well_IDs)) stop("Feature pseudobulks already pool wells; GEM_well_IDs filtering must occur before aggregation.")
      model$cell_type_subset <- normalize_pseudobulk_cell_type_subset(model$cell_type_subset)
      model
    })
}


normalize_pseudobulk_cell_type_subset <- function(cell_type_subset) {
  if (is.null(cell_type_subset)) {
    return(NULL)
  }

  stringr::str_replace_all(as.character(cell_type_subset), "_", "-")
}


get_pseudobulk_cell_type_subset_label <- function(cell_type_subset) {
  if (is.null(cell_type_subset)) {
    return("all")
  }

  stringr::str_c(cell_type_subset, collapse = "+")
}


get_pseudobulk_sample_tibble <- function(pseudobulk_data_matrix) {
  sample_tibble <- tibble::tibble(sample_name = colnames(pseudobulk_data_matrix)) |>
    tidyr::extract(
      col = sample_name,
      into = c("cluster", "donor_id"),
      regex = "^(.+)_([^_]+)$",
      remove = FALSE
    )

  assert_with_info(
    !anyNA(sample_tibble$cluster) && !anyNA(sample_tibble$donor_id),
    glue_info = "Pseudobulk sample names must follow '<cluster>_<donor_id>'."
  )

  sample_tibble
}


#' Get donor metadata columns required by configured differential analyses
#'
#' Collect donor-level variables from pseudobulk model formulas, pairing and
#' correlation fields, and differential cell-type composition settings.
#' @keywords internal

get_differential_analysis_metadata_columns <- function(models, abundance_models) {
  all_models <- c(models, abundance_models)
  unique(unlist(lapply(all_models, function(model) {
    c(get_differential_model_metadata_columns(model), model$plot_phenotype_vars, model$color_by)
  }), use.names = FALSE))
}

#' Filter pseudobulk data matrix
#'
#' Align and filter a pseudobulk feature matrix before model fitting.
#'
#' @param pseudobulk_data_matrix Feature-by-pseudobulk-sample count/activity matrix.
#' @param pseudobulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @param extended_donor_id_metadata_tibble Donor metadata after GEM well/sample-level covariates have been added for model design.
#' @param sample_depth_tibble Optional sample depth tibble used to remove low
#'   count pseudobulk columns before modeling.
#' @param min_sample_counts Minimum sample counts required when
#'   `sample_depth_tibble` is supplied.
#' @return Filtered feature matrix with columns aligned to model samples and
#'   features passing expression/support filters.
#' @keywords internal

filter_pseudobulk_data_matrix <- function(
  pseudobulk_data_matrix,
  pseudobulk_feature_dynamic_tibble,
  extended_donor_id_metadata_tibble,
  sample_depth_tibble = NULL,
  min_sample_counts = NULL
) {
  model_list <- pseudobulk_feature_dynamic_tibble$model[[1]]
  cell_type_subset <- model_list$cell_type_subset

  sample_tibble <- get_pseudobulk_sample_tibble(pseudobulk_data_matrix)
  cohort <- get_differential_model_cohort(extended_donor_id_metadata_tibble, model_list, sample_tibble$donor_id)
  donor_ids_keep <- stringr::str_replace_all(cohort$donor_id[cohort$included], "_", "-")
  sample_tibble <- dplyr::filter(sample_tibble, donor_id %in% donor_ids_keep)

  if (!is.null(cell_type_subset)) {
    sample_tibble <- sample_tibble |>
      dplyr::filter(.data$cluster %in% cell_type_subset)
  }

  if (!is.null(sample_depth_tibble) && !is.null(min_sample_counts)) {
    sample_tibble <- sample_tibble |>
      dplyr::inner_join(sample_depth_tibble |> dplyr::select(sample_name, n_counts), by = "sample_name") |>
      dplyr::filter(.data$n_counts >= min_sample_counts)
  }

  surviving_col_names <- dplyr::pull(sample_tibble, sample_name)

  assert_with_info(
    length(surviving_col_names) > 0,
    glue_info = "No pseudobulk samples remain after donor/model and cell-type-subset filtering."
  )

  filtered_pseudobulk_data_matrix <- pseudobulk_data_matrix[, surviving_col_names, drop = FALSE]
  if (inherits(filtered_pseudobulk_data_matrix, "IterableMatrix")) {
    filtered_pseudobulk_data_matrix <- methods::as(filtered_pseudobulk_data_matrix, "dgCMatrix")
  }
  assertthat::assert_that(is.matrix(filtered_pseudobulk_data_matrix) || inherits(filtered_pseudobulk_data_matrix, "dgCMatrix"))
  filtered_pseudobulk_data_matrix
}

#' Get pseudobulk chromVAR background record
#'
#' Fit the betterChromVAR background model for a pseudobulk count matrix.
#'
#' The expected accessibility defines the reference that deviations are
#' measured against and the accessibility on which background peaks are
#' matched. By default it pools the columns by depth, so a column holding most
#' reads is largely compared with itself. With equal column weights, each
#' column is depth-normalized before averaging, so every column is compared
#' with the average column profile.
#'
#' @param pseudobulk_ATAC_data_matrix Peak-by-pseudobulk-sample ATAC count matrix.
#' @param chromVAR_obj chromVAR SummarizedExperiment containing deviations, annotations, and background metadata.
#' @param equal_column_weights Whether each column contributes equally to the
#'   expected accessibility instead of in proportion to its depth.
#' @return A list containing the fitted background and retained peak names.
#' @keywords internal

get_pseudobulk_chromVAR_background_record <- function(
  pseudobulk_ATAC_data_matrix,
  chromVAR_obj,
  equal_column_weights = FALSE
) {
  peaks_above_cut_off_names <- pseudobulk_ATAC_data_matrix |>
    pseudobulk_rowSums() |>
    magrittr::is_greater_than(0) |>
    which() |>
    names()

  rowranges <- SummarizedExperiment::rowRanges(chromVAR_obj)
  peak_range_idx <- match(peaks_above_cut_off_names, get_peak_names_from_GRanges(rowranges))
  if (anyNA(peak_range_idx)) {
    stop("Some pseudobulk ATAC peaks were not found in single-cell chromVAR row ranges.")
  }
  filtered_rowranges <- rowranges[peak_range_idx]

  peak_filtered_pseudobulk_ATAC_data_matrix <- pseudobulk_ATAC_data_matrix[peaks_above_cut_off_names, , drop = FALSE]
  if (inherits(peak_filtered_pseudobulk_ATAC_data_matrix, "IterableMatrix")) {
    peak_filtered_pseudobulk_ATAC_data_matrix <- methods::as(peak_filtered_pseudobulk_ATAC_data_matrix, "dgCMatrix")
  }
  pseudobulk_chromVAR_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_filtered_pseudobulk_ATAC_data_matrix),
    rowRanges = filtered_rowranges
  )
  SummarizedExperiment::rowData(pseudobulk_chromVAR_obj) <- SummarizedExperiment::rowData(chromVAR_obj)[peak_range_idx, , drop = FALSE]

  expectation_vec <- betterChromVAR::getExpectation(pseudobulk_chromVAR_obj,
    grouping = if (equal_column_weights) colnames(pseudobulk_chromVAR_obj), normalize = TRUE)
  background_bins <- betterChromVAR::getBackgroundBins(
    x = expectation_vec,
    bias = SummarizedExperiment::rowData(pseudobulk_chromVAR_obj)$bias,
    flbias = SummarizedExperiment::rowData(pseudobulk_chromVAR_obj)$flbias,
    verbose = FALSE
  )
  background <- betterChromVAR::computeBackgrounds(
    object = peak_filtered_pseudobulk_ATAC_data_matrix,
    bins = background_bins,
    expectation = expectation_vec,
    verbose = FALSE
  )
  list(
    background = background,
    peak_names = rownames(peak_filtered_pseudobulk_ATAC_data_matrix)
  )
}

#' Pseudobulk motif-family accessibility
#'
#' Compute analytic chromVAR z-scores of motif families on pseudobulk ATAC
#' counts, centre each sample, and quantile-normalize across samples.
#' @keywords internal
get_pseudobulk_motif_family_accessibility_matrix <- function(
  pseudobulk_ATAC_data_matrix,
  chromVAR_obj,
  chromVAR_motif_family_matrix
) {
  background_record <- get_pseudobulk_chromVAR_background_record(
    pseudobulk_ATAC_data_matrix = pseudobulk_ATAC_data_matrix,
    chromVAR_obj = chromVAR_obj
  )
  peak_names <- background_record$peak_names
  counts_matrix <- pseudobulk_ATAC_data_matrix[peak_names, , drop = FALSE]
  if (inherits(counts_matrix, "IterableMatrix")) {
    counts_matrix <- methods::as(counts_matrix, "dgCMatrix")
  }
  rowranges <- SummarizedExperiment::rowRanges(chromVAR_obj)
  peak_range_idx <- match(peak_names, get_peak_names_from_GRanges(rowranges))
  pseudobulk_chromVAR_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowRanges = rowranges[peak_range_idx]
  )
  SummarizedExperiment::rowData(pseudobulk_chromVAR_obj) <- SummarizedExperiment::rowData(chromVAR_obj)[peak_range_idx, , drop = FALSE]
  chromVAR_z_scores <- betterChromVAR::computeDeviationsAnalytic(
    object = pseudobulk_chromVAR_obj,
    background = background_record$background,
    annotations = chromVAR_motif_family_matrix[peak_names, , drop = FALSE],
    verbose = FALSE,
    retSE = TRUE,
    compute = c("deviations", "z")
  ) |>
    SummarizedExperiment::assay("z")

  chromVAR_z_scores |>
    sweep(2, colMeans(chromVAR_z_scores), FUN = "-") |>
    limma::normalizeBetweenArrays(method = "quantile")
}

get_pseudobulk_cell_type_design <- function(sample_tibble, formula_chr) {
  design_matrix <- stats::model.matrix(stats::as.formula(formula_chr), data = sample_tibble)
  colnames(design_matrix) <- colnames(design_matrix) |>
    stringr::str_replace_all(":", ".") |>
    stringr::str_replace_all("-", "_")
  rownames(design_matrix) <- sample_tibble$ID
  design_matrix
}

#' Estimate the consensus intra-block correlation on evenly spaced features
#'
#' `limma::duplicateCorrelation()` fits a mixed model to every feature only to
#' average one consensus correlation, which takes hours for hundreds of thousands
#' of peaks. Evenly spaced features estimate it with negligible error; the linear
#' models are still fitted to every feature.
duplicate_correlation_on_feature_subset <- function(object, design, block, weights = NULL, max_features = 50000L) {
  n_features <- nrow(limma::getEAWP(object)$exprs)
  if (n_features <= max_features) {
    return(limma::duplicateCorrelation(object, design, block = block, weights = weights))
  }
  rows <- unique(round(seq(1, n_features, length.out = max_features)))
  limma::duplicateCorrelation(object[rows, , drop = FALSE], design, block = block,
    weights = if (is.matrix(weights)) weights[rows, , drop = FALSE] else weights)
}

#' `edgeR::voomLmFit()` with the intra-block correlation estimated on a feature subset
#'
#' voomLmFit calls duplicateCorrelation() by name, so a copy evaluated in a child
#' of its namespace uses the subset estimate without modifying edgeR code.
voom_lm_fit_subset_correlation <- function(counts, ..., max_correlation_features = 50000L) {
  fit <- edgeR::voomLmFit
  environment(fit) <- list2env(list(duplicateCorrelation = function(object, design, block, weights = NULL) {
    duplicate_correlation_on_feature_subset(object, design, block, weights, max_features = max_correlation_features)
  }), parent = environment(edgeR::voomLmFit))
  fit(counts, ...)
}

fit_pseudobulk_cell_type_matrix <- function(
  pseudobulk_feature_matrix,
  sample_tibble,
  formula_chr,
  correlation_block = NULL,
  voom_fit = voom_lm_fit_subset_correlation
) {
  design_matrix <- get_pseudobulk_cell_type_design(sample_tibble, formula_chr)
  is_count_data <- is_count_matrix(pseudobulk_feature_matrix)

  has_correlation_block <- !is.null(correlation_block) && nzchar(correlation_block)
  if (has_correlation_block) {
    assert_with_info(
      correlation_block %in% names(sample_tibble) && !anyNA(sample_tibble[[correlation_block]]),
      glue_info = "Configured cell-type correlation block is missing or incomplete: {correlation_block}."
    )
  }

  if (is_count_data) {
    gene_expression_list <- edgeR::DGEList(
      counts = pseudobulk_feature_matrix,
      samples = tibble::column_to_rownames(sample_tibble, var = "ID")
    )
    sample_cols_expressed <- pseudobulk_colSums(gene_expression_list$counts) > 0
    gene_expression_list <- gene_expression_list[, sample_cols_expressed, keep.lib.sizes = FALSE]
    design_matrix <- design_matrix[sample_cols_expressed, , drop = FALSE]
    feature_rows_expressed <- edgeR::filterByExpr(gene_expression_list, design = design_matrix)
    gene_expression_list <- edgeR::normLibSizes(gene_expression_list[feature_rows_expressed, , keep.lib.sizes = FALSE])

    assert_with_info(
      nrow(gene_expression_list) > 0 && nrow(design_matrix) > ncol(design_matrix),
      glue_info = "No valid cell-type-specific count model remains after filtering."
    )

    fit <- voom_fit(
      gene_expression_list,
      design = design_matrix,
      block = if (has_correlation_block) factor(gene_expression_list$samples[[correlation_block]]) else NULL,
      normalize.method = "none",
      keep.EList = TRUE
    )
    retained_samples <- gene_expression_list$samples
  } else {
    retained_samples <- tibble::column_to_rownames(sample_tibble, var = "ID")
    expression_object <- structure(
      list(E = pseudobulk_feature_matrix, samples = retained_samples),
      class = "EList"
    )
    if (has_correlation_block) {
      block <- factor(retained_samples[[correlation_block]])
      correlation_fit <- duplicate_correlation_on_feature_subset(
        expression_object,
        design = design_matrix,
        block = block
      )
      fit <- limma::lmFit(
        expression_object,
        design = design_matrix,
        block = block,
        correlation = correlation_fit$consensus.correlation
      )
    } else {
      fit <- limma::lmFit(expression_object, design = design_matrix)
    }
    fit$EList <- expression_object
  }

  fit$samples <- retained_samples
  list(
    fit = fit,
    design = design_matrix,
    samples = retained_samples,
    correlation_block = correlation_block,
    block_correlation = fit$correlation %||% NA_real_
  )
}

get_pseudobulk_standardized_residuals <- function(cell_type_record, feature_ids, donor_ids, random_effect) {
  sample_idx <- match(donor_ids, cell_type_record$samples[[random_effect]])
  feature_idx <- match(feature_ids, rownames(cell_type_record$fit$coefficients))
  expression_matrix <- cell_type_record$fit$EList$E[feature_idx, sample_idx, drop = FALSE]
  fitted_matrix <- cell_type_record$fit$coefficients[feature_idx, , drop = FALSE] %*%
    t(cell_type_record$design[sample_idx, , drop = FALSE])
  residual_matrix <- expression_matrix - fitted_matrix
  weights_matrix <- cell_type_record$fit$EList$weights
  if (!is.null(weights_matrix)) {
    residual_matrix <- residual_matrix * sqrt(weights_matrix[feature_idx, sample_idx, drop = FALSE])
  }
  residual_matrix - rowMeans(residual_matrix)
}

estimate_pseudobulk_cell_type_residual_correlations <- function(cell_type_records, random_effect, max_features = 2000L) {
  cell_types <- names(cell_type_records)
  correlation_matrix <- diag(length(cell_types))
  dimnames(correlation_matrix) <- list(cell_types, cell_types)
  common_feature_ids <- Reduce(
    intersect,
    purrr::map(cell_type_records, ~ rownames(.x$fit$coefficients))
  ) |>
    sort()

  assert_with_info(
    length(common_feature_ids) > 0,
    glue_info = "Cell-type-specific fits do not share any testable features."
  )

  if (length(common_feature_ids) > max_features) {
    common_feature_ids <- common_feature_ids[
      unique(round(seq(1, length(common_feature_ids), length.out = max_features)))
    ]
  }

  for (i in seq_len(length(cell_types) - 1L)) {
    for (j in seq.int(i + 1L, length(cell_types))) {
      cell_type_i <- cell_types[[i]]
      cell_type_j <- cell_types[[j]]
      common_donor_ids <- intersect(
        cell_type_records[[cell_type_i]]$samples[[random_effect]],
        cell_type_records[[cell_type_j]]$samples[[random_effect]]
      )
      assert_with_info(
        length(common_donor_ids) > ncol(cell_type_records[[cell_type_i]]$design) + 2L,
        glue_info = "Too few paired donors remain to estimate cross-cell-type residual correlation."
      )

      residuals_i <- get_pseudobulk_standardized_residuals(
        cell_type_records[[cell_type_i]], common_feature_ids, common_donor_ids, random_effect
      )
      residuals_j <- get_pseudobulk_standardized_residuals(
        cell_type_records[[cell_type_j]], common_feature_ids, common_donor_ids, random_effect
      )
      denominator <- sqrt(rowSums(residuals_i^2) * rowSums(residuals_j^2))
      feature_correlations <- rowSums(residuals_i * residuals_j) / denominator
      feature_correlations <- feature_correlations[is.finite(feature_correlations)]
      fisher_correlations <- atanh(pmax(-0.999999, pmin(0.999999, feature_correlations)))
      consensus_correlation <- tanh(mean(fisher_correlations, trim = 0.15))
      correlation_matrix[i, j] <- consensus_correlation
      correlation_matrix[j, i] <- consensus_correlation
    }
  }

  correlation_matrix
}

get_pseudobulk_stacked_cell_type_design <- function(cell_type_records, joint_design_matrix) {
  retained_sample_ids <- purrr::map(cell_type_records, ~ rownames(.x$design)) |>
    unlist(use.names = FALSE)
  separate_coefficient_names <- purrr::imap(
    cell_type_records,
    ~ paste(.y, colnames(.x$design), sep = "::")
  ) |>
    unlist(use.names = FALSE)
  stacked_design_matrix <- matrix(
    0,
    nrow = length(retained_sample_ids),
    ncol = length(separate_coefficient_names),
    dimnames = list(retained_sample_ids, separate_coefficient_names)
  )

  column_offset <- 0L
  for (cell_type in names(cell_type_records)) {
    cell_type_design <- cell_type_records[[cell_type]]$design
    column_idx <- seq.int(column_offset + 1L, column_offset + ncol(cell_type_design))
    stacked_design_matrix[rownames(cell_type_design), column_idx] <- cell_type_design
    column_offset <- max(column_idx)
  }

  retained_joint_design <- joint_design_matrix[retained_sample_ids, , drop = FALSE]
  joint_from_separate_matrix <- qr.solve(retained_joint_design, stacked_design_matrix)
  rownames(joint_from_separate_matrix) <- colnames(retained_joint_design)
  colnames(joint_from_separate_matrix) <- colnames(stacked_design_matrix)

  assert_with_info(
    max(abs(retained_joint_design %*% joint_from_separate_matrix - stacked_design_matrix)) < 1e-8,
    glue_info = "The joint and cell-type-specific pseudobulk formulas do not span the same model."
  )

  list(
    joint_design = retained_joint_design,
    stacked_design = stacked_design_matrix,
    joint_from_separate = joint_from_separate_matrix
  )
}

fit_pseudobulk_feature_matrix_by_cell_type <- function(
  pseudobulk_feature_matrix,
  final_sample_tibble,
  model_list,
  design_matrix,
  basis_matrix
) {
  pairing_variable <- model_list$pairing_variable %||% model_list$random_effect
  correlation_block <- model_list$correlation_block %||% NULL
  assert_with_info(
    !is.null(pairing_variable) && pairing_variable != "",
    glue_info = "cell_type_formula requires a configured donor/sample pairing_variable."
  )

  cell_types <- sort(unique(final_sample_tibble$cluster))
  assert_with_info(
    length(cell_types) > 1L,
    glue_info = "cell_type_formula requires at least two retained cell types."
  )

  fit_cell_type <- function(cell_type) {
    sample_idx <- which(final_sample_tibble$cluster == cell_type)
    cell_type_samples <- final_sample_tibble[sample_idx, , drop = FALSE]
    assert_with_info(
      !anyDuplicated(cell_type_samples[[pairing_variable]]),
      glue_info = "Cell-type-specific pseudobulk fits require at most one sample per pairing unit and cell type."
    )
    fit_pseudobulk_cell_type_matrix(
      pseudobulk_feature_matrix = pseudobulk_feature_matrix[, sample_idx, drop = FALSE],
      sample_tibble = cell_type_samples,
      formula_chr = model_list$cell_type_formula,
      correlation_block = correlation_block
    )
  }
  n_parallel_cell_type_fits <- min(6L, length(cell_types))

  cell_type_records <- if (n_parallel_cell_type_fits == 1L) {
    purrr::map(cell_types, fit_cell_type)
  } else {
    parallel::mclapply(
      cell_types,
      function(cell_type) {
        # Each fork makes many small BLAS calls; a multithreaded BLAS per fork
        # oversubscribes the cores and made these fits about 50 times slower.
        RhpcBLASctl::blas_set_num_threads(1L)
        fit_cell_type(cell_type)
      },
      mc.cores = n_parallel_cell_type_fits,
      mc.preschedule = TRUE
    )
  }
  child_errors <- purrr::keep(cell_type_records, ~ inherits(.x, "try-error"))
  assert_with_info(
    length(child_errors) == 0L,
    glue_info = "One or more parallel cell-type fits failed: {paste(child_errors, collapse = '; ')}"
  )
  cell_type_records <- cell_type_records |>
    purrr::set_names(cell_types)

  residual_correlations <- estimate_pseudobulk_cell_type_residual_correlations(
    cell_type_records = cell_type_records,
    random_effect = pairing_variable
  )
  stacked_design <- get_pseudobulk_stacked_cell_type_design(cell_type_records, design_matrix)
  cell_type_fits <- purrr::map(cell_type_records, function(record) {
    record$fit$EList <- NULL
    record$fit
  })
  retained_sample_ids <- rownames(stacked_design$joint_design)
  samples <- final_sample_tibble[match(retained_sample_ids, final_sample_tibble$ID), , drop = FALSE] |>
    tibble::column_to_rownames(var = "ID")

  structure(
    list(
      cell_type_fits = cell_type_fits,
      residual_correlations = residual_correlations,
      correlation_block = correlation_block,
      block_correlations = purrr::map_dbl(cell_type_records, "block_correlation"),
      n_parallel_cell_type_fits = n_parallel_cell_type_fits,
      joint_from_separate = stacked_design$joint_from_separate,
      design = stacked_design$joint_design,
      basis_matrix = basis_matrix,
      samples = samples,
      analysis_type = if (is.null(correlation_block) || !nzchar(correlation_block)) {
        "limma_cell_type_paired"
      } else {
        "limma_cell_type_paired_blocked"
      }
    ),
    class = c("pseudobulk_cell_type_fit", "list")
  )
}

#' Fit pseudobulk feature matrix model
#'
#' Fit the configured pseudobulk model for one feature matrix branch.
#'
#' @param pseudobulk_feature_matrix Feature-by-sample matrix to fit with voom/edgeR.
#' @param extended_donor_id_metadata_tibble Donor/sample metadata containing the
#'   covariates referenced by the configured model.
#' @param pseudobulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @return An edgeR `DGEGLM`, limma `MArrayLM`, or paired-cell-type
#'   `pseudobulk_cell_type_fit`, with sample/design information and an
#'   `analysis_type` field.
#' @keywords internal

fit_pseudobulk_feature_matrix_model <- function(pseudobulk_feature_matrix, extended_donor_id_metadata_tibble, pseudobulk_feature_dynamic_tibble) {
  formatted_extended_donor_id_metadata_tibble <- extended_donor_id_metadata_tibble |>
    dplyr::mutate(donor_id = stringr::str_replace_all(donor_id, "_", "-"))

  # set up sample tibble
  final_sample_tibble <- colnames(pseudobulk_feature_matrix) |>
    tibble::enframe(value = "ID", name = NULL) |>
    tidyr::separate(ID, into = c("cluster", "donor_id"), sep = "_", remove = FALSE) |>
    dplyr::left_join(formatted_extended_donor_id_metadata_tibble, by = "donor_id")

  # Set up design
  model_list <- pseudobulk_feature_dynamic_tibble$model[[1]]
  design_and_basis_matrices <- get_design_and_basis_matrices_from_model_list(sample_tibble = final_sample_tibble, model_list = model_list)
  design_matrix <- design_and_basis_matrices$design_matrix
  basis_matrix <- design_and_basis_matrices$basis_matrix

  rownames(design_matrix) <- final_sample_tibble$ID

  if (!is.null(model_list$cell_type_formula) && model_list$cell_type_formula != "") {
    return(fit_pseudobulk_feature_matrix_by_cell_type(
      pseudobulk_feature_matrix = pseudobulk_feature_matrix,
      final_sample_tibble = final_sample_tibble,
      model_list = model_list,
      design_matrix = design_matrix,
      basis_matrix = basis_matrix
    ))
  }

  # Determine analysis type based on data characteristics and model configuration
  is_count_data <- is_count_matrix(pseudobulk_feature_matrix)
  random_effect <- model_list$random_effect
  analysis_type <- if (!is.null(random_effect) && random_effect != "") {
    # Random effects require limma with duplicateCorrelation (works with any data type)
    "limma_duplicateCorrelation"
  } else if (is_count_data && ncol(design_matrix) > 1) {
    # EdgeR requires integer count matrix and at least 2 parameters
    "edgeR"
  } else {
    # Basic limma as fallback for continuous data or insufficient parameters
    "limma_basic"
  }

  if (is_count_data && analysis_type != "limma_basic") {
    gene_expression_list <- edgeR::DGEList(counts = pseudobulk_feature_matrix, samples = final_sample_tibble)
    sample_cols_expressed <- colSums(gene_expression_list$counts) > 0
    gene_expression_list <- gene_expression_list[, sample_cols_expressed, keep.lib.sizes = FALSE]
    design_matrix <- design_matrix[sample_cols_expressed, , drop = FALSE]
    feature_rows_expressed <- edgeR::filterByExpr(gene_expression_list, design = design_matrix)
    gene_expression_list <- edgeR::normLibSizes(gene_expression_list[feature_rows_expressed, , keep.lib.sizes = FALSE]) # consider using method="TMMwsp" that is designed for 0-inflated data. However, needs quantification of 0-inflatedness.
    if (!nrow(gene_expression_list$counts)) {
      stop("No features remain after expression filtering. Please check count depth and model complexity.")
    }
  } else {
    gene_expression_list <- structure(
      list(E = pseudobulk_feature_matrix, samples = tibble::column_to_rownames(final_sample_tibble, var = "ID")),
      class = "EList"
    )
  }
  validate_differential_design(design_matrix)

  gene_expression_list_fit <- if (analysis_type == "limma_duplicateCorrelation") {
    # Apply voom transformation only for count data; use EList directly for non-count data
    expression_object <- if (is_count_data) {
      limma::voom(gene_expression_list, design = design_matrix, plot = FALSE)
    } else {
      gene_expression_list
    }
    block <- gene_expression_list$samples[[random_effect]]
    corfit <- duplicate_correlation_on_feature_subset(expression_object, design = design_matrix, block = block)
    limma::lmFit(expression_object, design = design_matrix, correlation = corfit$consensus, block = block)
  } else if (analysis_type == "edgeR") {
    gene_expression_list |>
      edgeR::estimateDisp.DGEList(design = design_matrix) |>
      edgeR::glmQLFit(design = design_matrix, robust = TRUE)
  } else {
    limma::lmFit(gene_expression_list, design = design_matrix)
  }
  # add basis matrix possibly required for complex contrast functions
  gene_expression_list_fit |>
    magrittr::inset2("basis_matrix", basis_matrix) |>
    magrittr::inset2("analysis_type", analysis_type) |>
    magrittr::inset2("samples", gene_expression_list$samples)
}

#' Get design and basis matrices from model list
#'
#' Construct design and optional basis matrices from a model specification.
#'
#' @param sample_tibble Sample metadata used by the model formula or custom
#'   design-matrix function.
#' @param model_list Model specification list containing a formula or custom
#'   design-matrix function name.
#' @return A list with `design_matrix` and optional `basis_matrix`; formula-based
#'   designs sanitize interaction and dash characters in column names.
#' @keywords internal

get_design_and_basis_matrices_from_model_list <- function(sample_tibble, model_list) {
  # design matrix
  design_matrix_func_str <- model_list$design_matrix_func_name

  if (is.null(design_matrix_func_str) || design_matrix_func_str == "") {
    # Design matrix will be generated from formula.
    design_formula <- stats::as.formula(model_list$formula)
    design_matrix <- stats::model.matrix(object = design_formula, data = sample_tibble)
    colnames(design_matrix) <- colnames(design_matrix) |>
      stringr::str_replace_all(":", ".") |>
      stringr::str_replace_all("-", "_")
    design_and_basis_matrices <- list(design_matrix = design_matrix, basis_matrix = NULL)
  } else {
    # Design matrix will be generated from function.
    design_matrix_func <- get(design_matrix_func_str)
    design_and_basis_matrices <- design_matrix_func(sample_tibble)
  }

  return(design_and_basis_matrices)
}

normalize_contrast_result <- function(contrast_result, contrast_name) {
  if (is.numeric(contrast_result) && is.atomic(contrast_result)) {
    return(purrr::set_names(list(contrast_result), contrast_name))
  }

  if (!is.list(contrast_result)) {
    stop("Contrast function must return a numeric contrast vector or a list of numeric contrast vectors.")
  }

  if (is.null(names(contrast_result))) {
    if (length(contrast_result) == 1) {
      names(contrast_result) <- contrast_name
    } else {
      stop("Contrast functions returning multiple contrasts must return a named list.")
    }
  }

  contrast_result
}

#' Get contrast vec list
#'
#' Build model contrast vectors from literal specs or custom contrast functions.
#'
#' @param pseudobulk_feature_matrix_fit Fitted pseudobulk model object with design
#'   matrix and sample metadata.
#' @param model_list Model specification list containing literal contrast specs
#'   and/or custom contrast function names.
#' @return A named list whose elements preserve downstream branch or plot labels where applicable.
#' @keywords internal

get_contrast_vec_list <- function(pseudobulk_feature_matrix_fit, model_list) {
  contrast_specs <- model_list$contrast_specs_vec
  suppress_warnings_matching(
    expr = {
      contrast_vec_list <- list()

      if (!is.null(contrast_specs)) {
        sanitized_contrast_specs <- contrast_specs |>
          stringr::str_replace_all(":", ".") |>
          stringr::str_replace_all("([:alnum:]+)-([:alnum:]+)", "\\1_\\2") # TODO: Very rough... - replace dashes in cluster_names, but don't replace actual minus signs (which should have space in between)

        contrast_vec_list <- c(
          contrast_vec_list,
          limma::makeContrasts(contrasts = sanitized_contrast_specs, levels = pseudobulk_feature_matrix_fit$design) |>
            as.data.frame() |>
            as.list() |>
            purrr::set_names(names(contrast_specs))
        )
      }

      if (!is.null(model_list$contrast_functions)) {
        contrast_vec_list <- c(
          contrast_vec_list,
          model_list$contrast_functions |>
          as.list() |>
          purrr::imap(\(function_name, contrast_name) {
            get(function_name)(pseudobulk_feature_matrix_fit) |>
              normalize_contrast_result(contrast_name = contrast_name)
          }) |>
          purrr::flatten()
        )
      }

      assert_with_info(
        length(contrast_vec_list) > 0,
        glue_info = "Model must define contrast_specs_vec, contrast_functions, or both."
      )

      duplicated_contrast_names <- names(contrast_vec_list)[duplicated(names(contrast_vec_list))]
      if (length(duplicated_contrast_names) > 0) {
        stop("Duplicated contrast name(s): ", stringr::str_c(unique(duplicated_contrast_names), collapse = ", "))
      }

      contrast_vec_list
    },
    pattern = "Renaming (Intercept) to Intercept"
  )
}

#' Get pseudobulk feature model contrast support
#'
#' Count samples and donors supporting each side of a model contrast.
#'
#' @param pseudobulk_feature_matrix_fit Fitted pseudobulk model object with design
#'   matrix and sample metadata.
#' @param contrast_vec Named numeric contrast vector aligned to the model coefficient columns.
#' @return One-row tibble describing contrast support: sample counts, donor
#'   counts, positive/negative side counts, paired donors, and analysis type.
#' @keywords internal

get_pseudobulk_feature_model_contrast_support <- function(pseudobulk_feature_matrix_fit, contrast_vec) {
  design_matrix <- pseudobulk_feature_matrix_fit$design
  samples_tibble <- pseudobulk_feature_matrix_fit$samples |>
    tibble::as_tibble(rownames = ".sample_id")
  if (!"ID" %in% names(samples_tibble)) {
    samples_tibble <- samples_tibble |>
      dplyr::mutate(ID = .data$.sample_id)
  }
  contrast_score <- as.numeric(design_matrix %*% contrast_vec)
  eps <- sqrt(.Machine$double.eps)

  support_tibble <- samples_tibble |>
    dplyr::mutate(
      positive_side = contrast_score > eps,
      negative_side = contrast_score < -eps
    )

  donor_support_tibble <- support_tibble |>
    dplyr::summarise(
      has_positive_side = any(positive_side),
      has_negative_side = any(negative_side),
      .by = donor_id
    )
  n_positive_samples <- sum(support_tibble$positive_side)
  n_negative_samples <- sum(support_tibble$negative_side)

  tibble::tibble(
    n_samples = nrow(support_tibble),
    n_donors = dplyr::n_distinct(support_tibble$donor_id),
    n_positive_samples = n_positive_samples,
    n_negative_samples = n_negative_samples,
    min_group_n_samples = min(n_positive_samples, n_negative_samples),
    n_paired_donors = sum(donor_support_tibble$has_positive_side & donor_support_tibble$has_negative_side),
    analysis_type = pseudobulk_feature_matrix_fit$analysis_type %||% NA_character_
  )
}

get_pseudobulk_cell_type_contrasts <- function(pseudobulk_feature_matrix_fit, joint_contrast_vec) {
  separate_contrast_vec <- drop(t(pseudobulk_feature_matrix_fit$joint_from_separate) %*% joint_contrast_vec)
  names(separate_contrast_vec) <- colnames(pseudobulk_feature_matrix_fit$joint_from_separate)

  purrr::map(names(pseudobulk_feature_matrix_fit$cell_type_fits), function(cell_type) {
    coefficient_prefix <- paste0(cell_type, "::")
    coefficient_idx <- startsWith(names(separate_contrast_vec), coefficient_prefix)
    cell_type_contrast <- separate_contrast_vec[coefficient_idx]
    names(cell_type_contrast) <- substring(names(cell_type_contrast), nchar(coefficient_prefix) + 1L)
    cell_type_contrast
  }) |>
    purrr::set_names(names(pseudobulk_feature_matrix_fit$cell_type_fits)) |>
    purrr::keep(~ any(abs(.x) > sqrt(.Machine$double.eps)))
}

get_pseudobulk_cell_type_contrast_statistics <- function(cell_type_fit, contrast_vec) {
  contrast_fit <- cell_type_fit |>
    limma::contrasts.fit(contrast = contrast_vec) |>
    limma::eBayes()
  standard_error <- drop(contrast_fit$stdev.unscaled) * sqrt(contrast_fit$s2.post)

  tibble::tibble(
    feature_id = rownames(contrast_fit$coefficients),
    logFC = drop(contrast_fit$coefficients),
    AveExpr = rep_len(contrast_fit$Amean %||% NA_real_, nrow(contrast_fit)),
    standard_error = standard_error,
    df_total = rep_len(contrast_fit$df.total, nrow(contrast_fit)),
    t = drop(contrast_fit$t),
    PValue = drop(contrast_fit$p.value),
    B = drop(contrast_fit$lods)
  )
}

get_pseudobulk_paired_cell_type_contrast_statistics <- function(pseudobulk_feature_matrix_fit, cell_type_contrasts) {
  active_cell_types <- names(cell_type_contrasts)
  reference_contrast <- cell_type_contrasts[[1]]
  reference_unit <- reference_contrast / sqrt(sum(reference_contrast^2))
  contrast_scales <- purrr::map_dbl(cell_type_contrasts, ~ sum(.x * reference_unit))
  proportional_residuals <- purrr::map2_dbl(
    cell_type_contrasts,
    contrast_scales,
    ~ max(abs(.x - .y * reference_unit))
  )
  assert_with_info(
    max(proportional_residuals) < 1e-8,
    glue_info = "Cross-cell-type contrasts must compare the same within-cell-type coefficient combination."
  )

  cell_type_statistics <- purrr::map2(
    pseudobulk_feature_matrix_fit$cell_type_fits[active_cell_types],
    cell_type_contrasts,
    get_pseudobulk_cell_type_contrast_statistics
  )
  common_feature_ids <- Reduce(intersect, purrr::map(cell_type_statistics, ~ .x$feature_id))
  assert_with_info(
    length(common_feature_ids) > 0,
    glue_info = "No features are independently testable in every cell type required by the contrast."
  )
  cell_type_statistics <- purrr::map(
    cell_type_statistics,
    ~ .x[match(common_feature_ids, .x$feature_id), , drop = FALSE]
  )

  logFC <- Reduce(`+`, purrr::map(cell_type_statistics, ~ .x$logFC))
  variance <- Reduce(`+`, purrr::map(cell_type_statistics, ~ .x$standard_error^2))
  for (i in seq_len(length(active_cell_types) - 1L)) {
    for (j in seq.int(i + 1L, length(active_cell_types))) {
      correlation <- pseudobulk_feature_matrix_fit$residual_correlations[
        active_cell_types[[i]], active_cell_types[[j]]
      ]
      covariance_sign <- sign(contrast_scales[[i]] * contrast_scales[[j]])
      variance <- variance +
        2 * correlation * covariance_sign *
          cell_type_statistics[[i]]$standard_error * cell_type_statistics[[j]]$standard_error
    }
  }
  assert_with_info(
    all(is.finite(variance) & variance > 0),
    glue_info = "Paired cell-type contrast produced non-positive or non-finite variance."
  )

  standard_error <- sqrt(variance)
  t_statistic <- logFC / standard_error
  df_total <- Reduce(pmin, purrr::map(cell_type_statistics, ~ .x$df_total))
  PValue <- 2 * stats::pt(-abs(t_statistic), df = df_total)
  AveExpr <- Reduce(`+`, purrr::map(cell_type_statistics, ~ .x$AveExpr)) / length(cell_type_statistics)

  tibble::tibble(
    feature_id = common_feature_ids,
    logFC = logFC,
    AveExpr = AveExpr,
    standard_error = standard_error,
    df_total = df_total,
    t = t_statistic,
    PValue = PValue,
    B = NA_real_
  )
}

get_pseudobulk_cell_type_contrast_support <- function(
  pseudobulk_feature_matrix_fit,
  cell_type_contrasts,
  pairing_variable
) {
  active_cell_types <- names(cell_type_contrasts)
  sample_tibbles <- purrr::map(
    pseudobulk_feature_matrix_fit$cell_type_fits[active_cell_types],
    "samples"
  )
  assert_with_info(
    all(purrr::map_lgl(sample_tibbles, ~ pairing_variable %in% names(.x))),
    glue_info = "Pairing variable is missing from one or more cell-type fits: {pairing_variable}."
  )

  reference_contrast <- cell_type_contrasts[[1]]
  reference_unit <- reference_contrast / sqrt(sum(reference_contrast^2))
  contrast_scales <- purrr::map_dbl(
    cell_type_contrasts,
    ~ sum(.x * reference_unit)
  )
  sample_counts <- purrr::map_int(sample_tibbles, nrow)
  donor_ids <- purrr::map(sample_tibbles, ~ as.character(.x[[pairing_variable]]))

  tibble::tibble(
    n_samples = sum(sample_counts),
    n_donors = dplyr::n_distinct(unlist(donor_ids, use.names = FALSE)),
    n_positive_samples = sum(sample_counts[contrast_scales > 0]),
    n_negative_samples = sum(sample_counts[contrast_scales < 0]),
    min_group_n_samples = min(
      sum(sample_counts[contrast_scales > 0]),
      sum(sample_counts[contrast_scales < 0])
    ),
    n_paired_donors = if (length(donor_ids) > 1L) {
      length(Reduce(intersect, donor_ids))
    } else {
      0L
    },
    analysis_type = pseudobulk_feature_matrix_fit$analysis_type %||% NA_character_
  )
}

get_pseudobulk_cell_type_model_results <- function(
  pseudobulk_feature_matrix_fit,
  contrast_vec_list,
  pairing_variable = "donor_id"
) {
  purrr::imap(
    contrast_vec_list,
    function(contrast_vec, contrast_name) {
      cell_type_contrasts <- get_pseudobulk_cell_type_contrasts(pseudobulk_feature_matrix_fit, contrast_vec)
      statistics <- if (length(cell_type_contrasts) == 1L) {
        cell_type <- names(cell_type_contrasts)[[1]]
        get_pseudobulk_cell_type_contrast_statistics(
          pseudobulk_feature_matrix_fit$cell_type_fits[[cell_type]],
          cell_type_contrasts[[cell_type]]
        )
      } else {
        get_pseudobulk_paired_cell_type_contrast_statistics(
          pseudobulk_feature_matrix_fit,
          cell_type_contrasts
        )
      }

      statistics |>
        dplyr::mutate(
          FDR = stats::p.adjust(PValue, method = "BH"),
          contrast = contrast_name,
          n_tested_features = dplyr::n()
        ) |>
        dplyr::arrange(PValue) |>
        dplyr::bind_cols(
          get_pseudobulk_cell_type_contrast_support(
            pseudobulk_feature_matrix_fit = pseudobulk_feature_matrix_fit,
            cell_type_contrasts = cell_type_contrasts,
            pairing_variable = pairing_variable
          )
        )
    }
  )
}

#' Get pseudobulk feature model results
#'
#' Run configured contrasts and collect pseudobulk model test results.
#'
#' @param pseudobulk_feature_matrix_fit Fitted pseudobulk model object, either an
#'   edgeR GLM fit or limma/voom fit.
#' @param pseudobulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @return Long differential-result tibble with one row per feature/contrast,
#'   model metadata, and contrast-support diagnostics.
#' @keywords internal

get_pseudobulk_feature_model_results <- function(pseudobulk_feature_matrix_fit, pseudobulk_feature_dynamic_tibble) {
  model_list <- pseudobulk_feature_dynamic_tibble$model[[1]]
  contrast_vec_list <- get_contrast_vec_list(pseudobulk_feature_matrix_fit, model_list)

  test_table_list <- if (inherits(pseudobulk_feature_matrix_fit, "pseudobulk_cell_type_fit")) {
    get_pseudobulk_cell_type_model_results(
      pseudobulk_feature_matrix_fit,
      contrast_vec_list,
      pairing_variable = model_list$pairing_variable %||% model_list$random_effect %||% "donor_id"
    )
  } else if (class(pseudobulk_feature_matrix_fit) == "DGEGLM") {
    purrr::imap(
      contrast_vec_list,
      ~ edgeR::glmQLFTest(pseudobulk_feature_matrix_fit, contrast = .x)$table |>
        dplyr::rename(AveExpr = logCPM) |>
        tibble::rownames_to_column("feature_id") |>
        dplyr::mutate(FDR = stats::p.adjust(PValue, method = "BH"), contrast = .y) |>
        dplyr::bind_cols(get_pseudobulk_feature_model_contrast_support(pseudobulk_feature_matrix_fit, .x))
    )
  } else if (class(pseudobulk_feature_matrix_fit) == "MArrayLM") {
    purrr::imap(
      contrast_vec_list,
      ~ pseudobulk_feature_matrix_fit |>
        limma::contrasts.fit(contrast = .x) |>
        limma::eBayes() |>
        limma::topTable(coef = 1, number = Inf) |>
        dplyr::rename(PValue = P.Value, FDR = adj.P.Val) |>
        tibble::rownames_to_column("feature_id") |>
        dplyr::mutate(contrast = .y) |>
        dplyr::bind_cols(get_pseudobulk_feature_model_contrast_support(pseudobulk_feature_matrix_fit, .x))
    ) # format to EdgeR style
  }

  formatted_out <- test_table_list |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      cell_type_subset = pseudobulk_feature_dynamic_tibble$cell_type_subset,
      model = pseudobulk_feature_dynamic_tibble$model_name
    ) |>
    tibble::as_tibble()

  return(formatted_out)
}

#' Get pseudobulk differential significant elements tibble
#'
#' Extract significant pseudobulk differential features at an FDR threshold.
#'
#' @param combined_pseudobulk_differential_results_tibble Combined differential-results tibble across models and contrasts.
#' @param FDR_threshold Adjusted-P-value cutoff used to classify features as significant.
#' @param FDR_col Name of the adjusted-P-value column to threshold, commonly `FDR`.
#' @return A tibble with stable identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_pseudobulk_differential_significant_elements_tibble <- function(combined_pseudobulk_differential_results_tibble, FDR_threshold = 0.05, FDR_col = "FDR") {
  pseudobulk_differential_results_tibble <- combined_pseudobulk_differential_results_tibble |>
    dplyr::bind_rows()

  if (nrow(pseudobulk_differential_results_tibble) == 0) {
    return(tibble::tibble(model = character(), contrast = character(), direction = character(), n_significant = integer(), n_signed = integer()))
  }

  assert_with_info(
    FDR_col %in% names(pseudobulk_differential_results_tibble),
    glue_info = "pseudobulk_differential_results_tibble must contain the requested FDR column: {FDR_col}."
  )

  significant_elements_tibble <- pseudobulk_differential_results_tibble |>
    dplyr::filter(.data[[FDR_col]] < FDR_threshold, logFC != 0) |>
    dplyr::mutate(direction = dplyr::if_else(logFC > 0, "up", "down")) |>
    dplyr::count(model, contrast, direction, name = "n_significant")

  pseudobulk_differential_results_tibble |>
    dplyr::distinct(model, contrast) |>
    tidyr::expand_grid(direction = c("down", "up")) |>
    dplyr::left_join(significant_elements_tibble, by = c("model", "contrast", "direction")) |>
    dplyr::mutate(
      n_significant = tidyr::replace_na(n_significant, 0L),
      n_signed = dplyr::if_else(direction == "down", -n_significant, n_significant)
    )
}

plot_pseudobulk_differential_significant_elements <- function(significant_elements_tibble, modality = "features") {
  contrast_levels <- significant_elements_tibble |>
    dplyr::summarise(total_significant = sum(n_significant), .by = contrast) |>
    dplyr::arrange(total_significant, contrast) |>
    dplyr::pull(contrast)

  significant_elements_tibble |>
    dplyr::mutate(
      model = get_mixsorted_factor(model),
      contrast = factor(contrast, levels = contrast_levels)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = contrast, y = n_signed, fill = model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = \(x) abs(x)) +
    ggplot2::labs(
      title = paste("Differential", modality, "by contrast"),
      subtitle = "Bars to the right count increases; bars to the left count decreases in the positive contrast direction.\nCompare direction and yield across models; more discoveries alone do not imply a stronger biological effect.",
      caption = stringr::str_wrap("Counts include features with BH FDR < 0.05 and nonzero fitted effect, corrected within each model and contrast.\nTick labels show absolute counts; contrasts are ordered by total discoveries across models. Feature sets and sample support may differ between models.", width = 110),
      x = NULL, y = "Significant features (left: decreases; right: increases)", fill = "Model") +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot pseudobulk differential significant elements modality distribution
#'
#' Plot the modality/category composition of significant differential features.
#'
#' @param significant_elements_modality_distribution_tibble Summary tibble counting significant features by contrast/model/modality.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_pseudobulk_differential_significant_elements_modality_distribution <- function(significant_elements_modality_distribution_tibble) {
  if (nrow(significant_elements_modality_distribution_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  significant_elements_modality_distribution_tibble |>
    group_split_by("model") |>
    purrr::imap(\(model_tibble, model_name) {
      contrast_levels <- model_tibble |>
        dplyr::summarise(total_significant = sum(n_significant), .by = contrast) |>
        dplyr::arrange(total_significant, contrast) |>
        dplyr::pull(contrast)

      model_tibble |>
        dplyr::mutate(
          contrast = factor(contrast, levels = contrast_levels),
          modality = get_mixsorted_factor(stringr::str_replace_all(modality, "_", " "))
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = contrast, y = prop_significant, fill = modality)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
        ggplot2::geom_text(
          ggplot2::aes(label = dplyr::if_else(n_significant > 0, as.character(n_significant), "")),
          position = ggplot2::position_dodge(width = 0.8),
          hjust = -0.15,
          size = 3
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(labels = \(x) paste0(round(100 * x), "%"), expand = ggplot2::expansion(mult = c(0, 0.15))) +
        ggplot2::labs(
          x = NULL,
          y = "Share of significant elements within modality",
          fill = "Modality",
          title = paste("Distribution of discoveries across contrasts:", model_name),
          subtitle = stringr::str_wrap("Compare which contrasts account for most discoveries within each modality; this is not the fraction of tested features found significant.", width = 100),
          caption = stringr::str_wrap("Numerator: FDR-significant features in this contrast. Denominator: discoveries summed across contrasts within the same model and modality. Labels give counts; features can contribute to several contrasts. A modality with no discoveries is displayed as zero.", width = 110)
        ) +
        ggplot2::theme(legend.position = "bottom")
    })
}


get_cameraPR_statistic <- function(contrast_statistics) {
  signed_normal_score <- sign(contrast_statistics$logFC) * stats::qnorm(
    pmax(pmin(contrast_statistics$PValue, 1), .Machine$double.xmin) / 2,
    lower.tail = FALSE
  )
  if ("t" %in% names(contrast_statistics)) {
    dplyr::if_else(is.finite(contrast_statistics$t), contrast_statistics$t, signed_normal_score)
  } else {
    signed_normal_score
  }
}

#' Get gene-set enrichment results
#'
#' Run competitive limma cameraPR gene-set tests for pseudobulk contrasts.
#'
#' @param results_tibble Contrast results from `get_pseudobulk_feature_model_results()`.
#' @param gene_sets Named list of detected gene identifiers per gene set.
#' @param min_genes_per_set Minimum number of contrast-tested genes required per gene set.
#' @return Competitive cameraPR result tibble for each contrast and gene set.
#' @keywords internal

get_gene_set_enrichment_results <- function(results_tibble, gene_sets, min_genes_per_set = 10L) {
  contrast_statistics <- results_tibble |>
    dplyr::mutate(cameraPR_statistic = get_cameraPR_statistic(dplyr::pick(dplyr::everything()))) |>
    dplyr::filter(is.finite(cameraPR_statistic))

  contrast_statistics |>
    dplyr::group_split(contrast) |>
    purrr::map_dfr(function(statistics) {
      gene_indices <- limma::ids2indices(
        gene.sets = gene_sets,
        identifiers = statistics$feature_id,
        remove.empty = FALSE
      )
      gene_indices <- gene_indices[lengths(gene_indices) >= min_genes_per_set]
      if (!length(gene_indices)) {
        return(tibble::tibble())
      }

      limma::cameraPR(
        statistic = statistics$cameraPR_statistic,
        index = gene_indices,
        inter.gene.cor = 0.01,
        sort = FALSE
      ) |>
        tibble::as_tibble(rownames = "ID") |>
        dplyr::mutate(
          FDR = stats::p.adjust(PValue, method = "BH"),
          contrast = statistics$contrast[[1]],
          method = "cameraPR",
          `-log10(PValue)` = -log10(PValue),
          `-log10(FDR)` = -log10(FDR),
          cell_type_subset = statistics$cell_type_subset[[1]],
          model = statistics$model[[1]],
          color_category = dplyr::case_when(
            FDR < 0.05 & Direction == "Up" ~ "SigUP",
            FDR < 0.05 & Direction == "Down" ~ "SigDown",
            TRUE ~ "NS"
          )
        )
    })
}

#' Plot pseudobulk gene expression volcano
#'
#' Plot a pseudobulk differential gene-expression volcano for one contrast.
#' Optional feature labels affect displayed text only; joins retain feature IDs.
#'
#' @param pseudobulk_differential_results_tibble Differential result tibble with model,
#'   contrast, feature, statistics, and significance columns.
#' @param x_val Column name mapped to the volcano plot x-axis.
#' @param y_val Column name mapped to the volcano plot y-axis, usually a P-value/FDR-derived statistic.
#' @param pseudobulk_differential_top_features_tibble Optional top-feature annotations used for
#'   volcano labels.
#' @param pseudobulk_differential_top_feature_open_targets_evidence_tibble Optional Open Targets evidence
#'   annotations for labelled genes.
#' @param feature_labels Named character vector of display labels keyed by feature ID.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_pseudobulk_feature_volcano <- function(
  pseudobulk_differential_results_tibble,
  x_val = "logFC",
  y_val = c("log10pvalue", "log10FDR")[1],
  pseudobulk_differential_top_features_tibble = tibble::tibble(feature_id = character()),
  pseudobulk_differential_top_feature_open_targets_evidence_tibble = tibble::tibble(feature_id = character(), OT_GWAS_evidence = logical()),
  feature_labels = character()
) {
  test_results_formatted_tibble <- pseudobulk_differential_results_tibble |>
    dplyr::mutate(
      log10FDR = -log10(FDR),
      log10pvalue = -log10(PValue),
      log10FDR_effect = dplyr::case_when(logFC < 0 ~ -log10FDR, .default = log10FDR),
      log10pvalue_effect = dplyr::case_when(logFC < 0 ~ -log10pvalue, .default = log10pvalue)
    ) |>
    dplyr::ungroup()

  top_features <- pseudobulk_differential_top_features_tibble |>
    dplyr::filter(.data$contrast %in% unique(pseudobulk_differential_results_tibble$contrast)) |>
    dplyr::pull(feature_id) |>
    unique()

  OT_evidence_tibble <- pseudobulk_differential_top_feature_open_targets_evidence_tibble |>
    dplyr::select(feature_id, OT_GWAS_evidence) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  has_OT_evidence <- any(!is.na(OT_evidence_tibble$OT_GWAS_evidence))

  test_results_w_GWAS_evidence_tibble <- test_results_formatted_tibble |>
    dplyr::left_join(OT_evidence_tibble, by = "feature_id") |>
    dplyr::mutate(
      OT_GWAS_evidence = tidyr::replace_na(as.character(OT_GWAS_evidence), "Not tested"),
      feature_label = dplyr::coalesce(unname(feature_labels[as.character(feature_id)]), as.character(feature_id))
    )

  plot_tibble <- test_results_w_GWAS_evidence_tibble |>
    dplyr::select(dplyr::all_of(c(x_val, y_val)), FDR, OT_GWAS_evidence, feature_label, feature_id)
  label_tibble <- plot_tibble |>
    dplyr::filter(feature_id %in% top_features, FDR < 0.05)

  build_pseudobulk_feature_volcano_plot(
    plot_tibble, label_tibble, x_val, y_val, has_OT_evidence,
    title = stringr::str_c(unique(test_results_formatted_tibble$model), ": ", unique(test_results_formatted_tibble$contrast))
  )
}

build_pseudobulk_feature_volcano_plot <- function(plot_tibble, label_tibble, x_val, y_val, has_OT_evidence, title) {
  force(plot_tibble)
  force(label_tibble)
  force(x_val)
  force(y_val)
  force(has_OT_evidence)
  force(title)
  ggplot2::ggplot(plot_tibble, ggplot2::aes(x = .data[[x_val]], y = .data[[y_val]], shape = FDR < 0.05)) +
    {
      if (has_OT_evidence) {
        ggplot2::geom_point(ggplot2::aes(color = OT_GWAS_evidence), alpha = 0.5)
      } else {
        ggplot2::geom_point(alpha = 0.5)
      }
    } +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
    ggplot2::labs(title = paste("Differential feature evidence:", title),
      subtitle = stringr::str_wrap("Look for sizeable effects with FDR support; positive values follow the positive model-contrast direction.", width = 100),
      caption = stringr::str_wrap(paste("Points are tested features; filled points mark BH FDR < 0.05 within this model and contrast. The x-axis uses the fitted model effect (log2 fold change for gene counts; activity-scale difference for activity models).",
        "Labels are selected from the top nominal-p features and shown only when FDR < 0.05. Open Targets colours, when present, describe external evidence and do not prove causality."), width = 110),
      x = "Fitted effect (model scale)", y = if (y_val == "log10FDR") "-log10(BH FDR)" else "-log10(nominal p-value)",
      shape = "BH FDR < 0.05", colour = "Open Targets evidence") +
    ggrepel::geom_text_repel(
      data = label_tibble,
      ggplot2::aes(label = feature_label),
      size = 2,
      min.segment.length = 0,
      max.overlaps = 40
    ) +
    ggplot2::scale_x_continuous(limits = symmetric_limits) +
    ggplot2::theme(legend.position = "bottom")
}


plot_pseudobulk_feature_volcanoes <- function(
  pseudobulk_differential_results_tibble,
  pseudobulk_differential_top_features_tibble = tibble::tibble(feature_id = character()),
  pseudobulk_differential_top_feature_open_targets_evidence_tibble = tibble::tibble(feature_id = character(), OT_GWAS_evidence = logical()),
  feature_labels = character()
) {
  pseudobulk_differential_results_tibble |>
    group_split_by("contrast") |>
    purrr::map(
      ~ plot_pseudobulk_feature_volcano(
        pseudobulk_differential_results_tibble = .x,
        pseudobulk_differential_top_features_tibble = pseudobulk_differential_top_features_tibble,
        pseudobulk_differential_top_feature_open_targets_evidence_tibble = pseudobulk_differential_top_feature_open_targets_evidence_tibble,
        feature_labels = feature_labels
      )
    )
}

#' Plot a formatted pseudobulk chromatin accessibility volcano
#'
#' Construct a chromatin accessibility volcano from minimal point and label tibbles.
#'
#' @param plot_tibble Point data containing `logFC`, the requested y-value
#'   column, `type`, and `sig`.
#' @param label_tibble Thirty-row label data containing `logFC`, the requested
#'   y-value column, `type`, and `geneName`.
#' @param y_val Column name mapped to the volcano plot y-axis, usually a P-value/FDR-derived statistic.
#' @param title Plot title combining model and contrast.
#' @return A ggplot ready for saving or composition.
#' @keywords internal

plot_pseudobulk_chromatin_accessibility_volcano_tibble <- function(
  plot_tibble,
  label_tibble,
  y_val,
  title
) {
  force(plot_tibble)
  force(label_tibble)
  force(y_val)
  force(title)
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(x = logFC, y = .data[[y_val]], color = type)
  ) +
    ggplot2::geom_point(ggplot2::aes(alpha = sig), size = 0.5) +
    ggplot2::labs(
      title = paste("Differential chromatin accessibility:", title),
      subtitle = "Positive log2 fold changes indicate higher accessibility in the positive contrast direction.\nPoint opacity marks FDR significance; the dashed line uses the threshold on the displayed p-value scale.",
      caption = stringr::str_wrap(paste0(
        "Each point is a peak. Opacity: BH FDR < 0.05 across tested peaks within this model and contrast. Dashed line: ",
        if (y_val == "log10FDR") "FDR = 0.05." else "nominal p = 0.05, not an FDR cutoff.",
        "\nColours group genomic annotations. Up to 30 peaks are selected for labels by -log10(p) x absolute log2 fold change; gene labels denote annotation, not proven targets."), width = 110),
      x = "Log2 fold change", y = if (y_val == "log10FDR") "-log10(BH FDR)" else "-log10(nominal p-value)",
      alpha = "FDR significance", colour = "Genomic annotation"
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), lty = 2) +
    ggrepel::geom_text_repel(
      data = label_tibble,
      ggplot2::aes(label = geneName),
      size = 3,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    ggplot2::scale_x_continuous(limits = symmetric_limits) +
    ggplot2::scale_alpha_manual(values = c("Sig" = 1, "NS" = 0.1)) +
    ggplot2::theme(legend.position = "top") +
    ignore_aes_in_color_legend()
}

#' Plot pseudobulk chromatin accessibility volcano
#'
#' Plot differential chromatin-accessibility peaks as genomic volcano panels.
#'
#' @param pseudobulk_differential_results_tibble_per_contrast Differential accessibility result
#'   tibble for one contrast, including peak IDs, log fold changes, and statistics.
#' @param ATAC_consensus_peak_GRanges GRanges object containing ATAC consensus peak GRanges coordinates and metadata.
#' @param y_val Column name mapped to the volcano plot y-axis, usually a P-value/FDR-derived statistic.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_pseudobulk_chromatin_accessibility_volcano <- function(
  pseudobulk_differential_results_tibble_per_contrast,
  ATAC_consensus_peak_GRanges,
  y_val = c("log10pvalue", "log10FDR")
) {
  y_val <- y_val[[1]]
  annot_tibble <- ATAC_consensus_peak_GRanges |>
    tibble::as_tibble() |>
    dplyr::transmute(
      region_vec,
      txType,
      geneName = dplyr::case_when(is.na(geneName) ~ region_vec, .default = geneName)
    )

  combined_coef_test_results_formatted <- pseudobulk_differential_results_tibble_per_contrast |>
    dplyr::select(feature_id, logFC, PValue, contrast, model) |>
    dplyr::left_join(annot_tibble, by = c("feature_id" = "region_vec")) |>
    dplyr::mutate(
      FDR = stats::p.adjust(PValue, method = "BH"),
      sig = dplyr::case_when(FDR < 0.05 ~ "Sig", .default = "NS"),
      log10pvalue = -log10(PValue),
      log10pvalue_effect = dplyr::case_when(logFC < 0 ~ -log10pvalue, .default = log10pvalue),
      log10FDR = -log10(FDR),
      log10FDR_effect = dplyr::case_when(logFC < 0 ~ -log10FDR, .default = log10FDR),
      # distance_category = case_match(txType, "TSS" ~ 0, "promoter" ~ 1, "fiveUTR" ~ 2, "threeUTR" ~ 3, "CDS" ~ 4, "intron" ~ 5, "proximal" ~ 6, "distal" ~ 7, ), # deprecated in favor of below, doesn't really make sense
      type = dplyr::case_when(
        txType %in% c("TSS", "promoter", "fiveUTR", "threeUTR", "CDS", "intron") ~ "genic",
        txType %in% c("proximal") ~ "near_regulatory",
        txType %in% c("distal") ~ "distal_regulatory",
        txType %in% c("intergenic") ~ "intergenic"
      ),
      importance = log10pvalue * abs(logFC)
    )

  top_peaks_tibble <- combined_coef_test_results_formatted |>
    dplyr::arrange(dplyr::desc(importance)) |>
    dplyr::slice_head(n = 30)

  plot_pseudobulk_chromatin_accessibility_volcano_tibble(
    plot_tibble = combined_coef_test_results_formatted |>
      dplyr::select(logFC, dplyr::all_of(y_val), type, sig),
    label_tibble = top_peaks_tibble |>
      dplyr::select(logFC, dplyr::all_of(y_val), type, geneName),
    y_val = y_val,
    title = stringr::str_c(
      unique(combined_coef_test_results_formatted$model),
      ": ",
      unique(combined_coef_test_results_formatted$contrast)
    )
  )
}


plot_pseudobulk_chromatin_accessibility_volcanoes <- function(pseudobulk_differential_results_tibble, ATAC_consensus_peak_GRanges) {
  pseudobulk_differential_results_tibble |>
    group_split_by("contrast") |>
    purrr::map(
      ~ plot_pseudobulk_chromatin_accessibility_volcano(
        pseudobulk_differential_results_tibble_per_contrast = .x,
        ATAC_consensus_peak_GRanges = ATAC_consensus_peak_GRanges
      )
    )
}


plot_gene_set_enrichment_results <- function(gene_set_enrichment_results_tibble) {
  if (nrow(gene_set_enrichment_results_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  gene_set_enrichment_results_tibble |>
    group_split_by("contrast") |>
    purrr::map(plot_gene_set_enrichment_contrast_results)
}


#' Plot gene-set enrichment contrast results
#'
#' Plot top competitive cameraPR terms for one pseudobulk contrast.
#'
#' @param gene_set_enrichment_results_tibble gene-set enrichment result tibble with contrast, pathway, enrichment score, and adjusted P-value columns.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_gene_set_enrichment_contrast_results <- function(gene_set_enrichment_results_tibble) {
  plotting_tibble <- gene_set_enrichment_results_tibble |>
    dplyr::mutate(
      ID = stringr::str_remove(ID, "^HALLMARK_|^WP_|^GOBP_|^GOCC_|^GOMF_|^HP_|^REACTOME_|^PID_|^MODULE_|^GAVISH_3CA_|^KEGG_"),
      ID = stringr::str_replace_all(ID, "_", " "),
      ID = stringr::str_wrap(ID, width = 48),
      signed_log10_FDR = dplyr::if_else(
        Direction == "Up",
        -log10(pmax(FDR, .Machine$double.xmin)),
        log10(pmax(FDR, .Machine$double.xmin))
      )
    ) |>
    dplyr::slice_min(PValue, n = 25, with_ties = FALSE) |>
    dplyr::mutate(ID = stats::reorder(ID, signed_log10_FDR)) |>
    dplyr::select(ID, signed_log10_FDR, FDR, model, contrast)

  build_gene_set_enrichment_plot(plotting_tibble)
}

# Keep the plot environment separate from the full enrichment results.
build_gene_set_enrichment_plot <- function(plotting_tibble) {
  force(plotting_tibble)

  ggplot2::ggplot(
    plotting_tibble,
    ggplot2::aes(x = signed_log10_FDR, y = ID, color = FDR < 0.05)
  ) +
    ggplot2::geom_vline(
      xintercept = c(log10(0.05), -log10(0.05)),
      linetype = "dashed",
      color = "grey55"
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey55")) +
    ggplot2::labs(
      title = stringr::str_c(
        "Competitive pathway enrichment: ",
        unique(plotting_tibble$model),
        ": ",
        unique(plotting_tibble$contrast)
      ) |>
        stringr::str_replace_all("_", " ") |>
        stringr::str_wrap(width = 70),
      subtitle = stringr::str_wrap("Read direction from the sign and support from the magnitude; enrichment is relative to genes outside each set.", width = 100),
      caption = stringr::str_wrap("cameraPR results: up to 25 terms with the smallest nominal p-values are displayed, including non-significant terms. Position is signed -log10(FDR), not an enrichment effect size. Dashed lines and colour mark FDR = 0.05. Overlapping gene sets are not independent findings.", width = 110),
      x = paste(
        "Signed -log10(FDR)",
        "Positive values indicate enrichment among positive statistics",
        sep = "\n"
      ),
      y = NULL,
      color = "FDR < 0.05"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "bottom"
    )
}

plot_pseudobulk_p_value_distribution <- function(combined_pseudobulk_differential_results_tibble) {
  combined_pseudobulk_differential_results_tibble |>
    dplyr::transmute(PValue, model = get_mixsorted_factor(model)) |>
    build_pseudobulk_p_value_distribution_plot()
}

# Force the compact input so its promise cannot retain the preparation scope.
build_pseudobulk_p_value_distribution_plot <- function(plotting_tibble) {
  force(plotting_tibble)
  plotting_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = PValue, color = model)) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~model, scales = "free_y") +
    ggplot2::labs(title = "Nominal p-value distributions by differential model",
      subtitle = stringr::str_wrap("Inspect excess small p-values and unusual shapes; a left-hand peak can reflect signal or model misspecification.", width = 100),
      caption = stringr::str_wrap("Kernel densities pool feature tests and contrasts within each model. These are unadjusted p-values; smoothing can extend beyond 0-1. Facets use separate density scales, and density height is not a count of discoveries.", width = 110),
      x = "Nominal p-value", y = "Density", colour = "Model") +
    ggplot2::theme(legend.position = "bottom")
}

#' Get one MSigDB gene-set collection
#'
#' Fetch one selected MSigDB collection or subcollection for competitive tests.
#' Pathway size is filtered later against each contrast's tested-gene universe.
#'
#' @param organism_chr Pipeline organism identifier.
#' @param collection_chr MSigDB collection identifier.
#' @param subcollection_chr Optional MSigDB subcollection identifier.
#' @return Named list of gene-symbol vectors, one per gene set.
#' @keywords internal

get_msigdb_gene_sets <- function(
  organism_chr,
  collection_chr,
  subcollection_chr = NA_character_
) {
  species_chr <- switch(
    organism_chr,
    Homo_sapiens = "human",
    Mus_musculus = "mouse",
    stop("Unsupported organism for MSigDB gene sets: ", organism_chr)
  )
  msigdbr_args <- list(species = species_chr, collection = collection_chr)
  if (!is.na(subcollection_chr)) {
    msigdbr_args$subcollection <- subcollection_chr
  }

  gene_sets_tibble <- rlang::exec(msigdbr::msigdbr, !!!msigdbr_args) |>
    dplyr::transmute(gs_name, gene_id = gene_symbol) |>
    dplyr::filter(!is.na(gene_id)) |>
    dplyr::distinct() |>
    dplyr::summarise(genes = list(gene_id), .by = gs_name)

  rlang::set_names(gene_sets_tibble$genes, gene_sets_tibble$gs_name)
}


#' Convert between ENSEMBL and symbol
#'
#' Convert gene symbols to Ensembl IDs, or Ensembl IDs to unique gene symbols.
#'
#' @param gene_ID_vec Character vector of gene symbols or Ensembl IDs to convert using `gene_features_df`.
#' @param gene_features_df Gene metadata data frame with Ensembl IDs, symbols, coordinates, or feature names used for ID conversion and annotation.
#' @param strict Logical; when TRUE, fail on ambiguous or missing gene ID mappings instead of returning NA.
#' @return Named character vector mapping each input ID to the opposite ID type.
#'   When `strict = TRUE`, missing mappings error instead of being omitted.
#' @keywords internal

convert_between_ENSEMBL_and_symbol <- function(gene_ID_vec, gene_features_df, strict = FALSE) {
  is_ENS_format <- all(stringr::str_detect(gene_ID_vec, "^ENSG"))

  named_output_vec <- if (is_ENS_format) {
    gene_features_df |>
      dplyr::filter(id %in% gene_ID_vec) |>
      dplyr::select(id, gene_name_unique) |>
      tibble::deframe()
  } else {
    gene_features_df |>
      dplyr::filter(gene_name_unique %in% gene_ID_vec) |>
      dplyr::select(gene_name_unique, id) |>
      tibble::deframe()
  }

  missing_IDs_idx <- which(!gene_ID_vec %in% names(named_output_vec))
  if (length(missing_IDs_idx) > 0 && strict) {
    stop(stringr::str_glue("No ENSEMBL IDs found for {gene_ID_vec[missing_IDs_idx]}"))
  }
  named_output_vec
}


#' Get OT GWAS gene evidence tibble
#'
#' Query Open Targets for GWAS evidence linked to marker genes and an EFO trait.
#'
#' @param gene_symbols Character vector of gene symbols to map and query for Open Targets evidence.
#' @param gene_features_df Gene metadata data frame with Ensembl IDs, symbols, coordinates, or feature names used for ID conversion and annotation.
#' @param efo_id Optional EFO trait ID used when querying Open Targets gene evidence; `NULL` skips evidence lookup.
#' @return A tibble with stable identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_OT_GWAS_gene_evidence_tibble <- function(gene_symbols, gene_features_df, efo_id = NULL) {
  gene_symbols <- unique(gene_symbols)

  if (length(gene_symbols) == 0) {
    return(tibble::tibble(gene = character(), ensembl_id = character(), OT_GWAS_evidence = logical()))
  }

  ensembl_ids_vec <- convert_between_ENSEMBL_and_symbol(gene_symbols, gene_features_df)
  gene_evidence_tibble <- tibble::tibble(
    gene = gene_symbols,
    ensembl_id = unname(ensembl_ids_vec[gene_symbols]),
    OT_GWAS_evidence = NA
  )

  if (is.null(efo_id) || identical(efo_id, "")) {
    return(gene_evidence_tibble)
  }

  graph_ql_client <- ghql::GraphqlClient$new(url = "https://api.platform.opentargets.org/api/v4/graphql")

  gene_evidence_tibble |>
    dplyr::mutate(
      OT_GWAS_evidence = purrr::map2_lgl(
        ensembl_id,
        gene,
        ~ {
          if (is.na(.x)) {
            return(NA)
          }

          query_result <- tryCatch(
            get_data_from_exec_query(
              "queries/GWAS_association_with_gene.txt",
              variables_list = list(ensemblId = .x, efoId = efo_id, size = 1),
              graph_ql_client = graph_ql_client
            ),
            error = function(error) {
              warning(stringr::str_glue("Open Targets query failed for {.y}: {conditionMessage(error)}"))
              NULL
            }
          )

          if (is.null(query_result)) {
            return(NA)
          }

          evidence_count <- query_result$disease$gwasCredibleSets$count %||% 0
          evidence_count > 0
        }
      )
    )
}

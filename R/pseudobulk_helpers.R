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

psbulk_rowSums <- function(x) {
  if (inherits(x, "IterableMatrix")) {
    BPCells::rowSums(x)
  } else {
    Matrix::rowSums(x)
  }
}

psbulk_colSums <- function(x) {
  if (inherits(x, "IterableMatrix")) {
    BPCells::colSums(x)
  } else {
    Matrix::colSums(x)
  }
}

get_pseudobulk_depth_tibble <- function(psbulk_data_matrix, count_col = "n_counts") {
  n_counts <- as.numeric(psbulk_colSums(psbulk_data_matrix))
  n_features <- as.numeric(psbulk_colSums(psbulk_data_matrix > 0))

  get_psbulk_sample_tibble(psbulk_data_matrix) |>
    dplyr::mutate(
      !!count_col := n_counts[match(sample_name, colnames(psbulk_data_matrix))],
      n_features = n_features[match(sample_name, colnames(psbulk_data_matrix))],
      counts_per_feature = .data[[count_col]] / n_features,
      log10_n_counts = log10(.data[[count_col]]),
      log10_n_features = log10(n_features),
      log10_counts_per_feature = log10(counts_per_feature)
    )
}

get_group_pseudobulk_depth_tibble <- function(psbulk_data_matrix, group_col = "cluster", count_col = "n_counts") {
  n_counts <- as.numeric(psbulk_colSums(psbulk_data_matrix))
  n_features <- as.numeric(psbulk_colSums(psbulk_data_matrix > 0))

  tibble::tibble(
    !!group_col := colnames(psbulk_data_matrix),
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
      title = "Pseudobulk depth per cluster-donor sample"
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
      ggplot2::labs(caption = stringr::str_glue("Dashed line marks DTFA minimum ATAC count threshold: {min_ATAC_sample_counts}."))
  }

  plot
}

normalize_psbulk_feature_models <- function(cfg_psbulk_feature_models) {
  if (is.null(cfg_psbulk_feature_models)) {
    return(list())
  }

  assert_with_info(
    is.list(cfg_psbulk_feature_models) && !is.null(names(cfg_psbulk_feature_models)) && all(names(cfg_psbulk_feature_models) != ""),
    glue_info = "cfg_psbulk_feature_models must be a named list."
  )

  cfg_psbulk_feature_models |>
    purrr::discard(is.null) |>
    purrr::map(\(model_list) {
      model_list$cell_type_subset <- normalize_psbulk_cell_type_subset(model_list$cell_type_subset)
      model_list
    })
}


normalize_psbulk_cell_type_subset <- function(cell_type_subset) {
  if (is.null(cell_type_subset)) {
    return(NULL)
  }

  stringr::str_replace_all(as.character(cell_type_subset), "_", "-")
}


get_psbulk_cell_type_subset_label <- function(cell_type_subset) {
  if (is.null(cell_type_subset)) {
    return("all")
  }

  stringr::str_c(cell_type_subset, collapse = "+")
}


get_psbulk_sample_tibble <- function(psbulk_data_matrix) {
  sample_tibble <- tibble::tibble(sample_name = colnames(psbulk_data_matrix)) |>
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


#' Filter psbulk data matrix
#'
#' Align and filter a pseudobulk feature matrix before model fitting.
#'
#' @param psbulk_data_matrix Feature-by-pseudobulk-sample count/activity matrix.
#' @param psbulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @param extended_donor_id_metadata_tibble Donor metadata after reaction/sample-level covariates have been added for model design.
#' @param sample_depth_tibble Optional sample depth tibble used to remove low
#'   count pseudobulk columns before modeling.
#' @param min_sample_counts Minimum sample counts required when
#'   `sample_depth_tibble` is supplied.
#' @return Filtered feature matrix with columns aligned to model samples and
#'   features passing expression/support filters.
#' @keywords internal

filter_psbulk_data_matrix <- function(
  psbulk_data_matrix,
  psbulk_feature_dynamic_tibble,
  extended_donor_id_metadata_tibble,
  sample_depth_tibble = NULL,
  min_sample_counts = NULL
) {
  model_list <- psbulk_feature_dynamic_tibble$model[[1]]
  cell_type_subset <- model_list$cell_type_subset

  # filter the full donor_cluster-x-Gene matrix, such that

  # 1) donors with missing measurements (or missing entirely) that are part of the model model_formula are removed.
  model_vars_vec <- all.vars(stats::as.formula(paste0("~", model_list$formula)))

  donor_ids_keep <- extended_donor_id_metadata_tibble |>
    dplyr::filter(dplyr::if_all(.cols = dplyr::any_of(model_vars_vec), .fns = ~ !is.na(.x))) |>
    dplyr::pull(donor_id) |>
    stringr::str_replace_all("_", "-")

  sample_tibble <- get_psbulk_sample_tibble(psbulk_data_matrix) |>
    dplyr::filter(.data$donor_id %in% donor_ids_keep)

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

  filtered_psbulk_data_matrix <- psbulk_data_matrix[, surviving_col_names, drop = FALSE]
  if (inherits(filtered_psbulk_data_matrix, "IterableMatrix")) {
    filtered_psbulk_data_matrix <- methods::as(filtered_psbulk_data_matrix, "dgCMatrix")
  }
  assertthat::assert_that(is.matrix(filtered_psbulk_data_matrix) || inherits(filtered_psbulk_data_matrix, "dgCMatrix"))
  filtered_psbulk_data_matrix
}

#' Get pseudobulk chromVAR activity matrix
#'
#' Compute pseudobulk chromVAR activity from ATAC counts and peak annotations.
#'
#' @param psbulk_ATAC_data_matrix Peak-by-pseudobulk-sample ATAC count matrix.
#' @param chromVAR_obj chromVAR SummarizedExperiment containing deviations, annotations, and background metadata.
#' @param annotation_matrix Peak-by-feature annotation/weight matrix used by
#'   chromVAR.
#' @param normalize Logical; when TRUE, normalize activity/count matrices before returning them.
#' @return A named matrix-like object with rows and columns aligned to the input feature/cell identifiers.
#' @keywords internal

get_pseudobulk_chromVAR_deviation_record <- function(
  psbulk_ATAC_data_matrix,
  chromVAR_obj,
  annotation_matrix
) {
  peaks_above_cut_off_names <- psbulk_ATAC_data_matrix |>
    psbulk_rowSums() |>
    magrittr::is_greater_than(0) |>
    which() |>
    names()

  rowranges <- SummarizedExperiment::rowRanges(chromVAR_obj)
  peak_range_idx <- match(peaks_above_cut_off_names, get_peak_names_from_GRanges(rowranges))
  if (anyNA(peak_range_idx)) {
    stop("Some pseudobulk ATAC peaks were not found in single-cell chromVAR row ranges.")
  }
  filtered_rowranges <- rowranges[peak_range_idx]

  peak_filtered_psbulk_ATAC_data_matrix <- psbulk_ATAC_data_matrix[peaks_above_cut_off_names, , drop = FALSE]
  if (inherits(peak_filtered_psbulk_ATAC_data_matrix, "IterableMatrix")) {
    peak_filtered_psbulk_ATAC_data_matrix <- methods::as(peak_filtered_psbulk_ATAC_data_matrix, "dgCMatrix")
  }
  annotation_matrix_peak_filtered <- annotation_matrix[rownames(peak_filtered_psbulk_ATAC_data_matrix), , drop = FALSE]

  psbulk_chromVAR_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_filtered_psbulk_ATAC_data_matrix),
    rowRanges = filtered_rowranges
  )
  SummarizedExperiment::rowData(psbulk_chromVAR_obj) <- SummarizedExperiment::rowData(chromVAR_obj)[peak_range_idx, , drop = FALSE]

  expectation_vec <- betterChromVAR::getExpectation(psbulk_chromVAR_obj)
  background_bins <- betterChromVAR::getBackgroundBins(
    x = expectation_vec,
    bias = SummarizedExperiment::rowData(psbulk_chromVAR_obj)$bias,
    flbias = SummarizedExperiment::rowData(psbulk_chromVAR_obj)$flbias,
    verbose = FALSE
  )
  background <- betterChromVAR::computeBackgrounds(
    object = peak_filtered_psbulk_ATAC_data_matrix,
    bins = background_bins,
    expectation = expectation_vec,
    verbose = FALSE
  )
  deviation_SE <- betterChromVAR::computeDeviationsAnalytic(
    object = psbulk_chromVAR_obj,
    background = background,
    annotations = annotation_matrix_peak_filtered,
    verbose = FALSE,
    retSE = TRUE,
    compute = c("deviations", "z")
  )

  list(
    deviation_SE = deviation_SE,
    background = background,
    peak_names = rownames(peak_filtered_psbulk_ATAC_data_matrix)
  )
}

get_pseudobulk_chromVAR_activity_matrix <- function(
  psbulk_ATAC_data_matrix,
  chromVAR_obj,
  annotation_matrix,
  normalize = TRUE
) {
  chromVAR_dev <- get_pseudobulk_chromVAR_deviation_record(
    psbulk_ATAC_data_matrix = psbulk_ATAC_data_matrix,
    chromVAR_obj = chromVAR_obj,
    annotation_matrix = annotation_matrix
  )$deviation_SE
  chromVAR_z_scores <- SummarizedExperiment::assay(chromVAR_dev, "z")

  if (isFALSE(normalize)) {
    return(chromVAR_z_scores)
  }

  chromVAR_z_scores |>
    sweep(2, colMeans(chromVAR_z_scores), FUN = "-") |>
    limma::normalizeBetweenArrays(method = "quantile")
}

scale_within_vector <- function(x) {
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / x_sd)
}

format_cell_type_GWAS_chromVAR_deviations <- function(chromVAR_dev, GWAS_inputs_tibble, cell_type_support_tibble = NULL) {
  deviation_tibble <- SummarizedExperiment::assay(chromVAR_dev, "deviations") |>
    tibble::as_tibble(rownames = "GWAS_ID") |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "cluster", values_to = "deviation")

  z_tibble <- SummarizedExperiment::assay(chromVAR_dev, "z") |>
    tibble::as_tibble(rownames = "GWAS_ID") |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "cluster", values_to = "z")

  out <- deviation_tibble |>
    dplyr::left_join(z_tibble, by = c("GWAS_ID", "cluster")) |>
    dplyr::group_by(.data$GWAS_ID) |>
    dplyr::mutate(
      relative_deviation = scale_within_vector(.data$deviation),
      median_score = .data$relative_deviation
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      z_p = stats::pnorm(.data$z, lower.tail = FALSE),
      z_q = stats::p.adjust(.data$z_p, method = "BH"),
      support_label = dplyr::case_when(
        .data$z >= 3 ~ "**",
        .data$z >= 2 ~ "*",
        .default = ""
      ),
      score = .data$deviation,
      score_is_sig = .data$z >= 2
    ) |>
    add_GWAS_heatmap_categories(GWAS_tibble = GWAS_inputs_tibble)

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
      median_score,
      deviation,
      relative_deviation,
      z,
      z_p,
      z_q,
      support_label,
      dplyr::any_of(c("n_cells", "n_counts", "n_features", "counts_per_feature"))
    )
}

get_pseudobulk_activity_matrix.TFA <- function(
  psbulk_ATAC_data_matrix,
  chromVAR_obj,
  chromVAR_motif_matrix
) {
  get_pseudobulk_chromVAR_activity_matrix(
    psbulk_ATAC_data_matrix = psbulk_ATAC_data_matrix,
    chromVAR_obj = chromVAR_obj,
    annotation_matrix = chromVAR_motif_matrix,
    normalize = TRUE
  )
}

#' Format pseudobulk GWAS chromVAR Z-scores
#'
#' Convert a GWAS-by-pseudobulk-sample chromVAR activity matrix to the long
#' score tibble used by GWAS score boxplots.
#'
#' @param psbulk_GWAS_chromVAR_activity_matrix GWAS-by-pseudobulk-sample
#'   chromVAR Z-score matrix.
#' @param GWAS_inputs_tibble GWAS configuration tibble with labels, Open
#'   Targets study IDs, finemapping methods, and weighting modes.
#' @return A long tibble with GWAS, pseudobulk sample, cluster, donor, and score
#'   columns.
#' @keywords internal

format_psbulk_GWAS_chromVAR_ZScores <- function(psbulk_GWAS_chromVAR_activity_matrix, GWAS_inputs_tibble) {
  if (nrow(psbulk_GWAS_chromVAR_activity_matrix) == 0 || ncol(psbulk_GWAS_chromVAR_activity_matrix) == 0) {
    return(tibble::tibble(
      GWAS_ID = character(),
      Category = character(),
      variant_weighting_mode = character(),
      sample_name = character(),
      cluster = character(),
      donor_id = character(),
      score = numeric(),
      score_is_sig = logical()
    ))
  }

  sample_tibble <- get_psbulk_sample_tibble(psbulk_GWAS_chromVAR_activity_matrix)

  psbulk_GWAS_chromVAR_activity_matrix |>
    tibble::as_tibble(rownames = "GWAS_ID") |>
    tidyr::pivot_longer(-GWAS_ID, names_to = "sample_name", values_to = "score") |>
    dplyr::left_join(sample_tibble, by = "sample_name") |>
    dplyr::mutate(score_is_sig = abs(score) > 1.95) |>
    add_GWAS_heatmap_categories(GWAS_tibble = GWAS_inputs_tibble) |>
    dplyr::relocate(GWAS_ID, Category, variant_weighting_mode, sample_name, cluster, donor_id, score, score_is_sig)
}

#' Fit psbulk feature matrix model
#'
#' Fit the configured pseudobulk model for one feature matrix branch.
#'
#' @param psbulk_feature_matrix Feature-by-sample matrix to fit with voom/edgeR.
#' @param extended_donor_id_metadata_tibble Donor/sample metadata containing the
#'   covariates referenced by the configured model.
#' @param psbulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @return A named matrix-like object with rows and columns aligned to the input feature/cell identifiers.
#' @keywords internal

fit_psbulk_feature_matrix_model <- function(psbulk_feature_matrix, extended_donor_id_metadata_tibble, psbulk_feature_dynamic_tibble) {
  formatted_extended_donor_id_metadata_tibble <- extended_donor_id_metadata_tibble |>
    dplyr::mutate(donor_id = stringr::str_replace_all(donor_id, "_", "-"))

  # set up sample tibble
  final_sample_tibble <- colnames(psbulk_feature_matrix) |>
    tibble::enframe(value = "ID", name = NULL) |>
    tidyr::separate(ID, into = c("cluster", "donor_id"), sep = "_", remove = FALSE) |>
    dplyr::left_join(formatted_extended_donor_id_metadata_tibble, by = "donor_id")

  # Set up design
  model_list <- psbulk_feature_dynamic_tibble$model[[1]]
  design_and_basis_matrices <- get_design_and_basis_matrices_from_model_list(sample_tibble = final_sample_tibble, model_list = model_list)
  design_matrix <- design_and_basis_matrices$design_matrix
  basis_matrix <- design_and_basis_matrices$basis_matrix

  n_samples <- nrow(design_matrix)
  n_params <- ncol(design_matrix)
  residual_df <- n_samples - n_params
  assert_with_info(
    residual_df > 0,
    glue_info = "Residual degrees of freedom in model fitting are less than 1. Please reduce the number of parameters in the model. Current design matrix: {design_matrix}"
  )

  # Determine analysis type based on data characteristics and model configuration
  is_count_data <- is_count_matrix(psbulk_feature_matrix)
  has_random_effect <- !(is.null(model_list$random_effect) || model_list$random_effect == "")
  has_sufficient_params <- n_params > 1

  analysis_type <- if (has_random_effect) {
    # Random effects require limma with duplicateCorrelation (works with any data type)
    "limma_duplicateCorrelation"
  } else if (is_count_data && has_sufficient_params) {
    # EdgeR requires integer count matrix and at least 2 parameters
    "edgeR"
  } else {
    # Basic limma as fallback for continuous data or insufficient parameters
    "limma_basic"
  }

  if (analysis_type == "edgeR" || (analysis_type == "limma_duplicateCorrelation" && is_count_data)) {
    # EdgeR or limma_duplicateCorrelation with count data requires DGEList
    DGE_list <- edgeR::DGEList(counts = psbulk_feature_matrix, samples = final_sample_tibble)

    sample_cols_expressed <- DGE_list$counts |>
      colSums() |>
      magrittr::is_greater_than(0)

    DGE_list <- DGE_list[, sample_cols_expressed, keep.lib.sizes = FALSE]
    design_matrix <- design_matrix[sample_cols_expressed, , drop = FALSE]
    feature_rows_expressed <- edgeR::filterByExpr(DGE_list, design = design_matrix)

    DGE_list <- edgeR::normLibSizes(DGE_list[feature_rows_expressed, , keep.lib.sizes = FALSE]) # consider using method="TMMwsp" that is designed for 0-inflated data. However, needs quantification of 0-inflatedness.
    residual_df <- nrow(design_matrix) - ncol(design_matrix)

    assert_with_info(
      assertthat::not_empty(DGE_list$counts) && residual_df > 0,
      glue_info = "No valid count model remains after expression/sample filtering. Please check count depth and model complexity."
    )
  } else {
    # limma framework for continuous data (limma_basic or limma_duplicateCorrelation without count data)
    DGE_list <- {
      features_tibble <- tibble::enframe(rownames(psbulk_feature_matrix), value = "feature_id", name = NULL)

      # sanity check that sample order matches
      assert_with_info(
        identical(colnames(psbulk_feature_matrix), final_sample_tibble$ID),
        glue_info = "Sample order mismatch between psbulk_feature_matrix and final_sample_tibble$ID."
      )

      # build EList for limma
      structure(
        list(
          E = psbulk_feature_matrix,
          samples = tibble::column_to_rownames(final_sample_tibble, var = "ID"),
          genes = features_tibble
        ),
        class = "EList"
      )
    }
  }

  # Fit model based on analysis_type
  random_effect <- model_list$random_effect
  DGE_list_fit <- if (analysis_type == "limma_duplicateCorrelation") {
    # Limma with duplicateCorrelation for random effects
    # Apply voom transformation only for count data; use EList directly for non-count data
    expression_object <- if (is_count_data) {
      limma::voom(DGE_list, design = design_matrix, plot = FALSE)
    } else {
      DGE_list
    }
    corfit <- limma::duplicateCorrelation(expression_object, design = design_matrix, block = DGE_list$samples[[random_effect]])
    expression_object |>
      limma::lmFit(design = design_matrix, correlation = corfit$consensus, block = DGE_list$samples[[random_effect]]) |>
      magrittr::inset2("samples", DGE_list$samples) |>
      magrittr::inset2("data_matrix", dplyr::coalesce(DGE_list$E, DGE_list$counts))
  } else if (analysis_type == "edgeR") {
    # EdgeR for count data without random effects
    DGE_list |>
      edgeR::estimateDisp.DGEList(design = design_matrix) |>
      edgeR::glmQLFit(design = design_matrix, robust = TRUE)
  } else {
    # limma_basic as fallback for continuous data or insufficient parameters
    DGE_list |>
      limma::lmFit(design = design_matrix) |>
      magrittr::inset2("samples", DGE_list$samples) |>
      magrittr::inset2("data_matrix", DGE_list$E %||% DGE_list$counts)
  }
  # add basis matrix possibly required for complex contrast functions
  psbulk_feature_matrix_fit <- DGE_list_fit |>
    magrittr::inset2("basis_matrix", basis_matrix) |>
    magrittr::inset2("analysis_type", analysis_type) |>
    magrittr::inset2("samples", DGE_list$samples) |>
    magrittr::inset2("data_matrix", DGE_list$E %||% DGE_list$counts)
  return(psbulk_feature_matrix_fit)
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
#' @param psbulk_feature_matrix_fit Fitted pseudobulk model object with design
#'   matrix and sample metadata.
#' @param model_list Model specification list containing literal contrast specs
#'   and/or custom contrast function names.
#' @return A named list whose elements preserve downstream branch or plot labels where applicable.
#' @keywords internal

get_contrast_vec_list <- function(psbulk_feature_matrix_fit, model_list) {
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
          limma::makeContrasts(contrasts = sanitized_contrast_specs, levels = psbulk_feature_matrix_fit$design) |>
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
            get(function_name)(psbulk_feature_matrix_fit) |>
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

#' Get psbulk feature model contrast support
#'
#' Count samples and donors supporting each side of a model contrast.
#'
#' @param psbulk_feature_matrix_fit Fitted pseudobulk model object with design
#'   matrix and sample metadata.
#' @param contrast_vec Named numeric contrast vector aligned to the model coefficient columns.
#' @return One-row tibble describing contrast support: sample counts, donor
#'   counts, positive/negative side counts, paired donors, and analysis type.
#' @keywords internal

get_psbulk_feature_model_contrast_support <- function(psbulk_feature_matrix_fit, contrast_vec) {
  design_matrix <- psbulk_feature_matrix_fit$design
  samples_tibble <- psbulk_feature_matrix_fit$samples |>
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
    analysis_type = psbulk_feature_matrix_fit$analysis_type %||% NA_character_
  )
}

#' Get psbulk feature model results
#'
#' Run configured contrasts and collect pseudobulk model test results.
#'
#' @param psbulk_feature_matrix_fit Fitted pseudobulk model object, either an
#'   edgeR GLM fit or limma/voom fit.
#' @param psbulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @return Long differential-result tibble with one row per feature/contrast,
#'   model metadata, and contrast-support diagnostics.
#' @keywords internal

get_psbulk_feature_model_results <- function(psbulk_feature_matrix_fit, psbulk_feature_dynamic_tibble) {
  model_list <- psbulk_feature_dynamic_tibble$model[[1]]
  contrast_vec_list <- get_contrast_vec_list(psbulk_feature_matrix_fit, model_list)

  test_table_list <- if (class(psbulk_feature_matrix_fit) == "DGEGLM") {
    purrr::imap(
      contrast_vec_list,
      ~ edgeR::glmQLFTest(psbulk_feature_matrix_fit, contrast = .x)$table |>
        dplyr::rename(AveExpr = logCPM) |>
        tibble::rownames_to_column("feature_id") |>
        dplyr::mutate(FDR = stats::p.adjust(PValue, method = "BH"), contrast = .y) |>
        dplyr::bind_cols(get_psbulk_feature_model_contrast_support(psbulk_feature_matrix_fit, .x))
    )
  } else if (class(psbulk_feature_matrix_fit) == "MArrayLM") {
    purrr::imap(
      contrast_vec_list,
      ~ psbulk_feature_matrix_fit |>
        limma::contrasts.fit(contrast = .x) |>
        limma::eBayes() |>
        limma::topTable(coef = 1, number = Inf) |>
        dplyr::rename(PValue = P.Value, FDR = adj.P.Val) |>
        tibble::rownames_to_column("feature_id") |>
        dplyr::mutate(contrast = .y) |>
        dplyr::bind_cols(get_psbulk_feature_model_contrast_support(psbulk_feature_matrix_fit, .x))
    ) # format to EdgeR style
  }

  formatted_out <- test_table_list |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      cell_type_subset = psbulk_feature_dynamic_tibble$cell_type_subset,
      model = psbulk_feature_dynamic_tibble$model_name
    ) |>
    tibble::as_tibble()

  return(formatted_out)
}

#' Get psbulk DX significant elements tibble
#'
#' Extract significant pseudobulk differential features at an FDR threshold.
#'
#' @param combined_psbulk_DX_results_tibble Combined differential-results tibble across models and contrasts.
#' @param FDR_threshold Adjusted-P-value cutoff used to classify features as significant.
#' @param FDR_col Name of the adjusted-P-value column to threshold, commonly `FDR`.
#' @return A tibble with stable identifiers and derived columns consumed by downstream targets.
#' @keywords internal

get_psbulk_DX_significant_elements_tibble <- function(combined_psbulk_DX_results_tibble, FDR_threshold = 0.05, FDR_col = "FDR") {
  psbulk_DX_results_tibble <- combined_psbulk_DX_results_tibble |>
    dplyr::bind_rows()

  if (nrow(psbulk_DX_results_tibble) == 0) {
    return(tibble::tibble(model = character(), contrast = character(), direction = character(), n_significant = integer(), n_signed = integer()))
  }

  assert_with_info(
    FDR_col %in% names(psbulk_DX_results_tibble),
    glue_info = "psbulk_DX_results_tibble must contain the requested FDR column: {FDR_col}."
  )

  significant_elements_tibble <- psbulk_DX_results_tibble |>
    dplyr::filter(.data[[FDR_col]] < FDR_threshold, logFC != 0) |>
    dplyr::mutate(direction = dplyr::if_else(logFC > 0, "up", "down")) |>
    dplyr::count(model, contrast, direction, name = "n_significant")

  psbulk_DX_results_tibble |>
    dplyr::distinct(model, contrast) |>
    tidyr::expand_grid(direction = c("down", "up")) |>
    dplyr::left_join(significant_elements_tibble, by = c("model", "contrast", "direction")) |>
    dplyr::mutate(
      n_significant = tidyr::replace_na(n_significant, 0L),
      n_signed = dplyr::if_else(direction == "down", -n_significant, n_significant)
    )
}

plot_psbulk_DX_significant_elements <- function(significant_elements_tibble) {
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
    ggplot2::labs(x = NULL, y = "Significant elements", fill = "Model") +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot psbulk DX significant elements modality distribution
#'
#' Plot the modality/category composition of significant differential features.
#'
#' @param significant_elements_modality_distribution_tibble Summary tibble counting significant features by contrast/model/modality.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_psbulk_DX_significant_elements_modality_distribution <- function(significant_elements_modality_distribution_tibble) {
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
          modality = get_mixsorted_factor(modality)
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
          title = model_name
        ) +
        ggplot2::theme(legend.position = "bottom")
    })
}


#' Get gene idxs and annot per gene list
#'
#' Convert gene sets into model-row index lists for limma camera.
#'
#' @param DGE_list edgeR-style DGE/model list with counts or log-expression, sample metadata, and fitted model components.
#' @param gene_set_subcollection Named gene-set collection or subcollection used for GSEA.
#' @param min_genes_per_set Minimum number of genes from a gene set that must be present in the model feature universe.
#' @return Named list of integer row indices for gene sets meeting
#'   `min_genes_per_set`.
#' @keywords internal

get_gene_idxs_and_annot_per_gene_list <- function(DGE_list, gene_set_subcollection, min_genes_per_set = 5) {
  # convert to indices in relation to DGE_list
  gene_idx_in_DGE_list <- limma::ids2indices(gene.sets = gene_set_subcollection, identifiers = rownames(DGE_list), remove.empty = FALSE)

  set_names_w_fewer_than_n_genes <- gene_idx_in_DGE_list |>
    purrr::map_vec(length) |>
    magrittr::is_less_than(min_genes_per_set) |>
    which() |>
    names()
  if (length(set_names_w_fewer_than_n_genes) > 0) {
    warning(stringr::str_glue(
      "Removed the following sets because they had fewer than {min_genes_per_set} genes: {str_c(set_names_w_fewer_than_n_genes, collapse = ', ')}"
    ))
  }

  # assert_that()
  out_gene_idx_in_DGE_list <- gene_idx_in_DGE_list[names(gene_idx_in_DGE_list) %!in% set_names_w_fewer_than_n_genes]

  assert_with_info(
    length(out_gene_idx_in_DGE_list) > 0,
    glue_info = "No genesets were left after removing those with fewer than {min_genes_per_set} detected genes in the DGE list (after expression filtering).\nConsider other gene sets or lowering the min_genes_per_set in get_gene_idxs_and_annot_per_gene_list()."
  )
  return(out_gene_idx_in_DGE_list)
}

#' Get GSEA results
#'
#' Run limma camera gene-set tests for configured pseudobulk contrasts.
#'
#' @param psbulk_feature_matrix_fit Fitted pseudobulk model object used as input
#'   to limma camera.
#' @param gene_indices_per_gene_list Named list of integer feature indices for gene sets present in the fitted model matrix.
#' @param psbulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model, contrast, and feature matrix being processed.
#' @return Camera/GSEA result tibble for each contrast and gene set; random-effect
#'   models currently return an empty tibble with a warning.
#' @keywords internal

get_GSEA_results <- function(psbulk_feature_matrix_fit, gene_indices_per_gene_list, psbulk_feature_dynamic_tibble) {
  model_list <- psbulk_feature_dynamic_tibble$model[[1]]

  if (!(is.null(model_list$random_effect) || model_list$random_effect == "")) {
    warning("Camera results not implemented for random effect, returning empty tibble")
    return(tibble::tibble())
  }

  non_empty_set_IDs <- gene_indices_per_gene_list |>
    purrr::discard(~ length(.x) == 0) |>
    names()
  empty_custom_set_IDs <- gene_indices_per_gene_list |>
    purrr::keep(~ length(.x) == 0) |>
    names() |>
    stringr::str_subset("CUSTOM:")
  if (length(empty_custom_set_IDs) > 0) {
    warning(stringr::str_glue("Empty custom sets: {str_c(empty_custom_set_IDs, collapse = ', ')}"))
  }

  assert_with_info(
    ncol(psbulk_feature_matrix_fit$design) > 1,
    glue_info = "GSEA only supported for models with more than one parameter."
  )
  contrast_vec_list <- get_contrast_vec_list(psbulk_feature_matrix_fit, model_list)

  # Note: Random effect errors out here
  # but this works for limma object on TF activity??
  # we manually add FDR since it is missing in the case of single-gene sets.
  fry_results <- contrast_vec_list |>
    purrr::imap(
      \(contrast_vec, contrast_name) {
        edgeR::fry.DGEGLM(
          psbulk_feature_matrix_fit,
          index = gene_indices_per_gene_list[non_empty_set_IDs],
          design = psbulk_feature_matrix_fit$design,
          contrast = contrast_vec
        ) |>
          dplyr::mutate(PValue = PValue.Mixed, FDR = stats::p.adjust(PValue, "BH"), contrast = contrast_name) |>
          tibble::rownames_to_column(var = "ID")
      }
    ) |>
    dplyr::bind_rows()

  camera_results <- contrast_vec_list |>
    purrr::imap(
      \(contrast_vec, contrast_name) {
        edgeR::camera.DGEGLM(
          psbulk_feature_matrix_fit,
          index = gene_indices_per_gene_list[non_empty_set_IDs],
          design = psbulk_feature_matrix_fit$design,
          contrast = contrast_vec
        ) |>
          dplyr::mutate(FDR = stats::p.adjust(PValue, "BH"), contrast = contrast_name) |>
          tibble::rownames_to_column(var = "ID")
      }
    ) |>
    dplyr::bind_rows()

  excluded_by_camera <- non_empty_set_IDs[non_empty_set_IDs %!in% camera_results$ID]
  excluded_by_fry <- non_empty_set_IDs[non_empty_set_IDs %!in% fry_results$ID]

  if (length(excluded_by_camera) > 0 || length(excluded_by_fry) > 0) {
    camera_custom_gene_sets <- stringr::str_subset(excluded_by_camera, "CUSTOM:")
    fry_custom_gene_sets <- stringr::str_subset(excluded_by_fry, "CUSTOM:")
    stop(stringr::str_glue(
      "CAMERA excluded {length(excluded_by_camera)} gene sets, including the custom sets: {str_c(camera_custom_gene_sets, collapse = ', ')}",
      "FRY excluded {length(excluded_by_fry)} gene sets, including the custom sets: {str_c(fry_custom_gene_sets, collapse = ', ')}",
      .sep = "\n"
    ))
  }

  combined_results <-
    list(fry = fry_results, camera = camera_results) |>
    dplyr::bind_rows(.id = "method") |>
    dplyr::mutate(
      `-log10(PValue)` = -log10(PValue),
      `-log10(FDR)` = -log10(FDR),
      cell_type_subset = psbulk_feature_dynamic_tibble$cell_type_subset,
      model = psbulk_feature_dynamic_tibble$model_name,
      color_category = dplyr::case_when(
        FDR < 0.05 & Direction == "Up" ~ "SigUP",
        FDR < 0.05 & Direction == "Down" ~ "SigDown",
        TRUE ~ "NS"
      )
    )

  return(combined_results)
}

#' Plot psbulk DGE volcano
#'
#' Plot a pseudobulk differential gene-expression volcano for one contrast.
#'
#' @param psbulk_DX_results_tibble Differential result tibble with model,
#'   contrast, feature, statistics, and significance columns.
#' @param x_val Column name mapped to the volcano plot x-axis.
#' @param y_val Column name mapped to the volcano plot y-axis, usually a P-value/FDR-derived statistic.
#' @param psbulk_DX_top_features_tibble Optional top-feature annotations used for
#'   volcano labels.
#' @param psbulk_DX_top_feature_OT_evidence_tibble Optional Open Targets evidence
#'   annotations for labelled genes.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_psbulk_DGE_volcano <- function(
  psbulk_DX_results_tibble,
  x_val = "logFC",
  y_val = c("log10pvalue", "log10FDR")[1],
  psbulk_DX_top_features_tibble = tibble::tibble(feature_id = character()),
  psbulk_DX_top_feature_OT_evidence_tibble = tibble::tibble(feature_id = character(), OT_GWAS_evidence = logical())
) {
  test_results_formatted_tibble <- psbulk_DX_results_tibble |>
    dplyr::mutate(
      log10FDR = -log10(FDR),
      log10pvalue = -log10(PValue),
      log10FDR_effect = dplyr::case_when(logFC < 0 ~ -log10FDR, .default = log10FDR),
      log10pvalue_effect = dplyr::case_when(logFC < 0 ~ -log10pvalue, .default = log10pvalue)
    ) |>
    dplyr::ungroup()

  top_features <- psbulk_DX_top_features_tibble |>
    dplyr::filter(.data$contrast %in% unique(psbulk_DX_results_tibble$contrast)) |>
    dplyr::pull(feature_id) |>
    unique()

  OT_evidence_tibble <- psbulk_DX_top_feature_OT_evidence_tibble |>
    dplyr::select(feature_id, OT_GWAS_evidence) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  has_OT_evidence <- any(!is.na(OT_evidence_tibble$OT_GWAS_evidence))

  test_results_w_GWAS_evidence_tibble <- test_results_formatted_tibble |>
    dplyr::left_join(OT_evidence_tibble, by = "feature_id") |>
    dplyr::mutate(OT_GWAS_evidence = tidyr::replace_na(as.character(OT_GWAS_evidence), "Not tested"))

  volcano_plot <- ggplot2::ggplot(test_results_w_GWAS_evidence_tibble, ggplot2::aes(x = .data[[x_val]], y = .data[[y_val]], shape = FDR < 0.05)) +
    {
      if (has_OT_evidence) {
        ggplot2::geom_point(ggplot2::aes(color = OT_GWAS_evidence), alpha = 0.5)
      } else {
        ggplot2::geom_point(alpha = 0.5)
      }
    } +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
    ggplot2::labs(title = stringr::str_c(unique(test_results_formatted_tibble$model), ": ", unique(test_results_formatted_tibble$contrast))) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(test_results_w_GWAS_evidence_tibble, feature_id %in% top_features, FDR < 0.05),
      ggplot2::aes(label = feature_id),
      size = 2,
      min.segment.length = 0,
      max.overlaps = 40
    ) +
    ggplot2::scale_x_continuous(limits = symmetric_limits) +
    ggplot2::theme(legend.position = "bottom")
  # scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +

  volcano_plot
}


plot_psbulk_DGE_volcanoes <- function(
  psbulk_DX_results_tibble,
  psbulk_DX_top_features_tibble = tibble::tibble(feature_id = character()),
  psbulk_DX_top_feature_OT_evidence_tibble = tibble::tibble(feature_id = character(), OT_GWAS_evidence = logical())
) {
  psbulk_DX_results_tibble |>
    group_split_by("contrast") |>
    purrr::map(
      ~ plot_psbulk_DGE_volcano(
        psbulk_DX_results_tibble = .x,
        psbulk_DX_top_features_tibble = psbulk_DX_top_features_tibble,
        psbulk_DX_top_feature_OT_evidence_tibble = psbulk_DX_top_feature_OT_evidence_tibble
      )
    )
}

#' Plot psbulk DCA volcano
#'
#' Plot differential chromatin-accessibility peaks as genomic volcano panels.
#'
#' @param psbulk_DX_results_tibble_per_contrast Differential accessibility result
#'   tibble for one contrast, including peak IDs, log fold changes, and statistics.
#' @param ATAC_consensus_peak_GRanges GRanges object containing ATAC consensus peak GRanges coordinates and metadata.
#' @param y_val Column name mapped to the volcano plot y-axis, usually a P-value/FDR-derived statistic.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_psbulk_DCA_volcano <- function(
  psbulk_DX_results_tibble_per_contrast,
  ATAC_consensus_peak_GRanges,
  y_val = c("log10pvalue", "log10FDR")
) {
  annot_tibble <- ATAC_consensus_peak_GRanges |>
    tibble::as_tibble() |>
    dplyr::mutate(geneName = dplyr::case_when(is.na(geneName) ~ region_vec, .default = geneName))

  combined_coef_test_results_formatted <- psbulk_DX_results_tibble_per_contrast |>
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

  volcano_plot <- ggplot2::ggplot(combined_coef_test_results_formatted, ggplot2::aes(x = logFC, y = .data[[y_val[1]]], color = type)) + #txType)) +
    ggplot2::geom_point(ggplot2::aes(alpha = sig), size = 0.5) +
    ggplot2::labs(title = stringr::str_c(unique(combined_coef_test_results_formatted$model), ": ", unique(combined_coef_test_results_formatted$contrast))) +
    ggplot2::geom_hline(yintercept = -log10(0.05), lty = 2) +
    ggrepel::geom_text_repel(data = top_peaks_tibble, ggplot2::aes(label = geneName), size = 3, min.segment.length = 0, max.overlaps = 50) +
    ggplot2::scale_x_continuous(limits = symmetric_limits) +
    ggplot2::scale_alpha_manual(values = c("Sig" = 1, "NS" = 0.1)) +
    ggplot2::theme(legend.position = "top") +
    ignore_aes_in_color_legend()

  suppress_warnings_matching(
    expr = {
      print(volcano_plot)
      volcano_plot
    },
    pattern = "too many overlaps"
  )
}


plot_psbulk_DCA_volcanoes <- function(psbulk_DX_results_tibble, ATAC_consensus_peak_GRanges) {
  psbulk_DX_results_tibble |>
    group_split_by("contrast") |>
    purrr::map(
      ~ plot_psbulk_DCA_volcano(
        psbulk_DX_results_tibble_per_contrast = .x,
        ATAC_consensus_peak_GRanges = ATAC_consensus_peak_GRanges
      )
    )
}


plot_GSEA_results <- function(GSEA_results_tibble) {
  if (nrow(GSEA_results_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  GSEA_results_tibble |>
    group_split_by("contrast") |>
    purrr::map(plot_GSEA_contrast_results)
}


#' Plot GSEA contrast results
#'
#' Plot top camera/GSEA terms for one pseudobulk contrast.
#'
#' @param GSEA_results_tibble GSEA result tibble with contrast, pathway, enrichment score, and adjusted P-value columns.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_GSEA_contrast_results <- function(GSEA_results_tibble) {
  # not used in interactive mode

  sanitized_ID_GSEA_results_tibble <- GSEA_results_tibble |>
    dplyr::mutate(
      ID = stringr::str_remove(ID, "_TARGET_GENES$|^HALLMARK_|^KEGG_MEDICUS_(REFERENCE_)?|^WP_|^GOBP_|^GOCC_|^GOMF_|^HP_|^REACTOME_|^PID_|^MODULE_|^GAVISH_3CA_|^KEGG_")
    )

  top_ID_tibble <- sanitized_ID_GSEA_results_tibble |>
    dplyr::slice_min(order_by = PValue, n = 25, by = c(model, method, contrast)) |>
    dplyr::distinct(ID, model, contrast)

  plotting_tibble <- sanitized_ID_GSEA_results_tibble |>
    tidyr::pivot_wider(names_from = method, values_from = c(`-log10(FDR)`, Direction), id_cols = c(ID, model, contrast)) |>
    dplyr::mutate(
      `-log10(FDR)_camera` = dplyr::case_when(
        Direction_camera == "Up" ~ `-log10(FDR)_camera`,
        Direction_camera == "Down" ~ -`-log10(FDR)_camera`,
        .default = 0
      ),
      color_cat = dplyr::case_when(
        `-log10(FDR)_camera` > -log10(0.05) ~ "CameraSigUp",
        `-log10(FDR)_camera` < log10(0.05) ~ "CameraSigDown",
        `-log10(FDR)_fry` > -log10(0.05) ~ "SigMixed",
        .default = "NS"
      )
    )

  label_tibble <- plotting_tibble |>
    dplyr::inner_join(top_ID_tibble)

  plot <- ggplot2::ggplot(plotting_tibble, ggplot2::aes(x = `-log10(FDR)_camera`, y = `-log10(FDR)_fry`, label = ID, text = ID, color = color_cat)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = stringr::str_c(unique(plotting_tibble$model), ": ", unique(plotting_tibble$contrast))) +
    ggplot2::geom_vline(xintercept = -log10(0.05), lty = 2, color = "red") +
    ggplot2::geom_vline(xintercept = log10(0.05), lty = 2, color = "blue") +
    ggplot2::geom_hline(yintercept = -log10(0.05), lty = 2) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_color_manual(values = c("CameraSigUp" = "red", "CameraSigDown" = "blue", "SigMixed" = "black", "NS" = "grey")) +
    ggrepel::geom_text_repel(data = label_tibble, ggplot2::aes(label = ID), size = 2, min.segment.length = 0, max.overlaps = 20) +
    ggplot2::scale_x_continuous(limits = symmetric_limits)

  return(plot)
}

plot_psbulk_DX_PValue_density <- function(combined_psbulk_DX_results_tibble) {
  combined_psbulk_DX_results_tibble |>
    dplyr::mutate(model = get_mixsorted_factor(model)) |>
    ggplot2::ggplot(ggplot2::aes(x = PValue, color = model)) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~model, scales = "free_y") +
    ggplot2::theme(legend.position = "bottom")
}

#' Format psbulk GWAS chromVAR results
#'
#' Add GWAS metadata and plot labels to pseudobulk GWAS chromVAR results.
#'
#' @param psbulk_GWAS_chromVAR_results_tibble Raw model-result tibble(s) for
#'   GWAS chromVAR features.
#' @param GWAS_inputs_tibble GWAS configuration tibble with labels, Open Targets study IDs, finemapping methods, and weighting modes.
#' @return Formatted GWAS chromVAR result tibble with GWAS metadata, labels, and
#'   standardized result columns.
#' @keywords internal

format_psbulk_GWAS_chromVAR_results <- function(psbulk_GWAS_chromVAR_results_tibble, GWAS_inputs_tibble) {
  combined_results_tibble <- psbulk_GWAS_chromVAR_results_tibble |>
    dplyr::bind_rows()

  if (nrow(combined_results_tibble) == 0) {
    return(tibble::tibble(
      GWAS_ID = character(),
      Category = character(),
      variant_weighting_mode = character(),
      model = character(),
      contrast = character(),
      logFC = numeric(),
      AveExpr = numeric(),
      t = numeric(),
      PValue = numeric(),
      FDR = numeric(),
      FDR_within_GWAS_category = numeric(),
      cell_type_subset = character(),
      n_samples = integer(),
      n_donors = integer(),
      n_paired_donors = integer(),
      min_group_n_samples = integer(),
      analysis_type = character()
    ))
  }

  combined_results_tibble |>
    dplyr::rename(GWAS_ID = feature_id) |>
    dplyr::left_join(
      GWAS_inputs_tibble |> dplyr::select(Category, GWAS_ID, variant_weighting_mode),
      by = "GWAS_ID"
    ) |>
    dplyr::mutate(
      FDR_within_GWAS_category = stats::p.adjust(PValue, method = "BH"),
      .by = c(Category, GWAS_ID)
    ) |>
    dplyr::relocate(GWAS_ID, Category, variant_weighting_mode, model, contrast)
}

#' Plot psbulk GWAS chromVAR coefficient ranges
#'
#' Plot pseudobulk GWAS chromVAR model coefficients and approximate moderated
#' standard-error intervals by contrast within each GWAS row.
#'
#' @param psbulk_GWAS_chromVAR_results_tibble Formatted GWAS chromVAR results
#'   containing GWAS, model, contrast, and statistic columns.
#' @param contrast_compartment_patterns Named patterns used to group contrast labels before plotting.
#' @return A named list of ggplots, one per model.
#' @keywords internal

plot_psbulk_GWAS_chromVAR_coefficient_ranges <- function(psbulk_GWAS_chromVAR_results_tibble, contrast_compartment_patterns = NULL) {
  if (nrow(psbulk_GWAS_chromVAR_results_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  GWAS_levels <- psbulk_GWAS_chromVAR_results_tibble |>
    dplyr::distinct(Category, GWAS_ID) |>
    dplyr::arrange(Category, GWAS_ID) |>
    dplyr::distinct(GWAS_ID, .keep_all = TRUE) |>
    dplyr::pull(GWAS_ID) |>
    rev()

  plot_tibble <- psbulk_GWAS_chromVAR_results_tibble |>
    dplyr::mutate(
      model = get_mixsorted_factor(model),
      GWAS_ID = factor(GWAS_ID, levels = GWAS_levels),
      coefficient_se = dplyr::if_else(is.finite(t) & t != 0, abs(logFC / t), NA_real_),
      coefficient_low = logFC - 1.96 * coefficient_se,
      coefficient_high = logFC + 1.96 * coefficient_se,
      coefficient_class = factor(
        dplyr::case_when(
          FDR_within_GWAS_category < 0.05 & logFC > 0 ~ "Positive FDR < 0.05",
          FDR_within_GWAS_category < 0.05 & logFC < 0 ~ "Negative FDR < 0.05",
          .default = "FDR >= 0.05"
        ),
        levels = c("Positive FDR < 0.05", "Negative FDR < 0.05", "FDR >= 0.05")
      )
    )

  plot_tibble |>
    group_split_by("model") |>
    purrr::imap(\(model_tibble, model_name) {
      contrast_metadata <- get_compartment_metadata(
        model_tibble$contrast,
        contrast_compartment_patterns,
        type_col = "contrast",
        default_compartment = "Contrast"
      )
      contrast_levels <- contrast_metadata |>
        dplyr::arrange(compartment, contrast) |>
        dplyr::pull(contrast)
      contrast_breaks <- contrast_metadata |>
        dplyr::arrange(match(contrast, contrast_levels)) |>
        dplyr::pull(compartment) |>
        get_plot_group_breaks()

      coefficient_limits <- model_tibble |>
        dplyr::summarise(
          contrast = contrast_levels[[1]],
          coefficient_limit = max(abs(c(logFC, coefficient_low, coefficient_high)), na.rm = TRUE),
          .by = GWAS_ID
        ) |>
        dplyr::mutate(coefficient_limit = pmax(coefficient_limit, 1e-6)) |>
        tidyr::uncount(2, .id = "limit_direction") |>
        dplyr::mutate(
          coefficient_limit = dplyr::if_else(limit_direction == 1L, -coefficient_limit, coefficient_limit),
          contrast = factor(contrast, levels = contrast_levels)
        )

      model_tibble |>
        dplyr::mutate(contrast = factor(contrast, levels = contrast_levels)) |>
        ggplot2::ggplot(ggplot2::aes(x = contrast, y = logFC, color = coefficient_class)) +
        ggplot2::geom_blank(
          data = coefficient_limits,
          ggplot2::aes(x = contrast, y = coefficient_limit),
          inherit.aes = FALSE
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
        ggplot2::geom_vline(xintercept = contrast_breaks, color = "grey70", linewidth = 0.25) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = coefficient_low, ymax = coefficient_high), linewidth = 0.45, na.rm = TRUE) +
        ggplot2::geom_point(size = 1.8) +
        ggplot2::facet_grid(GWAS_ID ~ ., scales = "free_y") +
        ggplot2::scale_y_continuous(breaks = scales::breaks_extended(n = 3)) +
        ggplot2::scale_color_manual(
          values = c("Positive FDR < 0.05" = "#3B4CC0", "Negative FDR < 0.05" = "#B40426", "FDR >= 0.05" = "grey55"),
          name = NULL
        ) +
        ggplot2::labs(x = NULL, y = "Model coefficient (logFC)") +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          panel.border = ggplot2::element_rect(color = "grey50", fill = NA, linewidth = 0.35),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "bottom",
          strip.text.y = ggplot2::element_blank()
        )
    })
}

#' Plot psbulk GWAS chromVAR QC
#'
#' Plot QC summaries for pseudobulk GWAS chromVAR result distributions.
#'
#' @param psbulk_GWAS_chromVAR_results_tibble Formatted GWAS chromVAR results
#'   used for p-value and feature-support QC plots.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_psbulk_GWAS_chromVAR_QC <- function(psbulk_GWAS_chromVAR_results_tibble) {
  if (nrow(psbulk_GWAS_chromVAR_results_tibble) == 0) {
    return(structure(list(), class = c("empty_plot_list", "list")))
  }

  support_tibble <- psbulk_GWAS_chromVAR_results_tibble |>
    dplyr::distinct(model, contrast, n_samples, n_donors, min_group_n_samples, n_paired_donors)

  list(
    sample_counts = support_tibble |>
      tidyr::pivot_longer(
        cols = c(n_samples, n_donors, min_group_n_samples, n_paired_donors),
        names_to = "metric",
        values_to = "value"
      ) |>
      dplyr::mutate(metric = get_mixsorted_factor(metric)) |>
      ggplot2::ggplot(ggplot2::aes(x = contrast, y = value, fill = metric)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
      ggplot2::facet_wrap(~model, scales = "free_x") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom") +
      ggplot2::labs(x = NULL, y = "Count", fill = NULL, title = "GWAS chromVAR pseudobulk model support"),
    PValue_density = psbulk_GWAS_chromVAR_results_tibble |>
      plot_psbulk_DX_PValue_density()
  )
}

#' Plot DCTC by phenotype per cluster
#'
#' Plot donor cell-type composition against a phenotype for each cluster.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param DCTC_plot_phenotype_vars Phenotype columns to plot against predicted or observed cell-type composition.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param categorical Logical; TRUE treats phenotype variables as categorical, FALSE as continuous.
#' @param DCTC_color_by_categorical_metadata_column Optional categorical metadata column used for point colors in DCTC plots.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for `save_plots_structured()` or composition.
#' @keywords internal

plot_DCTC_by_phenotype_per_cluster <- function(
  metadata_tibble,
  DCTC_plot_phenotype_vars,
  cluster_col,
  categorical = FALSE,
  DCTC_color_by_categorical_metadata_column = NULL
) {
  w_pheno_class_tibble <-
    if (is_non_binary_numeric_vec(metadata_tibble[[DCTC_plot_phenotype_vars]]) && categorical) {
      get_grouped_phenotype(metadata_tibble, DCTC_plot_phenotype_vars)
    } else {
      dplyr::mutate(metadata_tibble, pheno_class = .data[[DCTC_plot_phenotype_vars]])
    }

  plot_tibble <- w_pheno_class_tibble |>
    dplyr::mutate(
      DTCT_color_by_col = as.factor(if (!is.null(DCTC_color_by_categorical_metadata_column)) .data[[DCTC_color_by_categorical_metadata_column]] else 1)
    ) |>
    dplyr::group_by(.data[[cluster_col]], donor_id, pheno_class, DTCT_color_by_col) |>
    dplyr::tally() |>
    dplyr::ungroup() |>
    dplyr::group_by(donor_id, pheno_class, DTCT_color_by_col) |>
    dplyr::mutate(prop = n / sum(n))

  is_continuous <- is_non_binary_numeric_vec(plot_tibble$pheno_class)

  geom_list <- if (is_continuous) {
    list(
      ggplot2::geom_point(alpha = 0.5),
      ggplot2::geom_smooth(method = "lm", formula = y ~ stats::poly(x, 2), linetype = "dashed")
    )
  } else {
    list(
      ggplot2::geom_boxplot(outlier.shape = NA),
      ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.5),
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    )
  }

  ggplot2::ggplot(plot_tibble, ggplot2::aes(x = pheno_class, y = prop, color = DTCT_color_by_col)) +
    geom_list +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[[cluster_col]]),
      cols = ggplot2::vars(DTCT_color_by_col)
    ) + # migth make sense to use facet_wrap to allow free scales - then limits below should also be removed.
    ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
    ggplot2::labs(x = stringr::str_glue("Phenotype class ({DCTC_plot_phenotype_vars})"), y = "Proportion of nuclei per cluster", title = stringr::str_glue("{DCTC_plot_phenotype_vars}"))
}

#' Get MSigDB genesets split by subcollection list list
#'
#' Fetch MSigDB gene sets and split them by subcollection for GSEA targets.
#'
#' @param species_chr Species label passed to MSigDB/msigdbr gene-set lookup.
#' @param use_id_chr Gene identifier type requested from MSigDB, such as symbols or Ensembl IDs.
#' @param min_size_num Minimum gene-set size retained before downstream enrichment testing.
#' @param max_size_num Maximum gene-set size retained before downstream enrichment testing.
#' @return Named list of subcollections; each subcollection is a named list of
#'   gene vectors after size filtering.
#' @keywords internal

get_MSigDB_genesets_split_by_subcollection_list_list <- function(
  species_chr,
  use_id_chr = c("gene_symbol", "db_gene_symbol", "ensembl_gene", "source_gene")[1],
  min_size_num = 15,
  max_size_num = 500
) {
  msig_tibble <- msigdbr::msigdbr(species = species_chr)

  filtered_msig_tibble <- msig_tibble |>
    dplyr::mutate(gs_subcollection = stringr::str_replace_all(dplyr::case_when(gs_subcollection == "" ~ gs_collection, .default = gs_subcollection), ":", "_")) |>
    dplyr::transmute(gs_name, gene_id = .data[[use_id_chr]], gs_subcollection) |>
    dplyr::filter(!is.na(gene_id)) |>
    dplyr::distinct() |>
    dplyr::add_count(gs_name, name = "n_genes") |>
    dplyr::filter(n_genes >= min_size_num, n_genes <= max_size_num)

  # list of character vectors per set, then to index list over y_DGEList rows
  MSigDB_genesets_split_by_subcollection_list_list <- filtered_msig_tibble |>
    group_split_by("gs_subcollection") |>
    purrr::map(
      \(gene_set_tibble) {
        summarised_tibble <- gene_set_tibble |>
          dplyr::group_by(gs_name) |>
          dplyr::summarise(genes = list(gene_id), .groups = "drop")

        rlang::set_names(summarised_tibble$genes, summarised_tibble$gs_name)
      }
    )

  MSigDB_genesets_split_by_subcollection_list_list
}


#' Get DCTC model results
#'
#' Fit donor cell-type composition models from cell metadata and donor covariates.
#'
#' @param metadata_w_cell_types_tibble Cell metadata tibble after cell-type labels have been joined; used for composition or export helpers.
#' @param extended_donor_id_metadata_tibble Donor metadata after reaction/sample-level covariates have been added for model design.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param DCTC_formula_chr Model formula string used for donor/cell-type composition testing.
#' @return Donor cell-type composition model results for the requested formula
#'   and cluster column.
#' @keywords internal

get_DCTC_model_results <- function(
  metadata_w_cell_types_tibble,
  extended_donor_id_metadata_tibble,
  cluster_col,
  DCTC_formula_chr
) {
  assert_cfg_is_set(cluster_col)
  assert_cfg_is_set(DCTC_formula_chr, "cfg_DCTC_formula_chr")
  formula <- stats::as.formula(DCTC_formula_chr)

  count_tibble_w_phenotypes <- metadata_w_cell_types_tibble |>
    dplyr::mutate(cluster = .data[[cluster_col]]) |>
    dplyr::count(donor_id, cluster, name = "n_nuclei") |>
    dplyr::group_by(donor_id) |>
    dplyr::mutate(n_total_nuclei = sum(n_nuclei)) |>
    dplyr::ungroup() |>
    dplyr::mutate(n_other_nuclei = n_total_nuclei - n_nuclei) |>
    dplyr::left_join(extended_donor_id_metadata_tibble, by = "donor_id") |>
    dplyr::mutate(n_donors = dplyr::n_distinct(donor_id), .by = "cluster") |>
    dplyr::filter(n_donors > 1)

  results <- count_tibble_w_phenotypes |>
    group_split_by("cluster") |>
    purrr::imap(
      \(cluster_tibble, cluster_name) {
        glmmTMB::glmmTMB(
          formula = formula,
          family = glmmTMB::betabinomial(link = "logit"),
          data = cluster_tibble
        ) |>
          broom.mixed::tidy(effects = "fixed")
      }
    ) |>
    dplyr::bind_rows(.id = "cluster")

  baseline_frac_tibble <- count_tibble_w_phenotypes |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(baseline_frac = sum(n_nuclei) / sum(n_total_nuclei))

  results_w_percent_change <- results |>
    dplyr::left_join(baseline_frac_tibble) |>
    dplyr::mutate(
      new_frac = stats::plogis(stats::qlogis(baseline_frac) + estimate * 1),
      middle_frac = (new_frac + baseline_frac) / 2,
      frac_diff = new_frac - baseline_frac,
      relative_frac_diff = frac_diff / baseline_frac,
      model_name = "DCTC"
    )

  return(results_w_percent_change)
}


plot_DCTC_model_change_per_unit <- function(results_tibble) {
  results_tibble |>
    dplyr::filter(!stringr::str_detect(term, "Intercept")) |>
    ggplot2::ggplot(ggplot2::aes(x = baseline_frac, y = cluster, fill = cluster, color = cluster, alpha = p.value < 0.05)) +
    # geom_col(position = "dodge2", width = 1) +
    ggplot2::geom_segment(ggplot2::aes(xend = new_frac), linewidth = 2, lineend = 'round', linejoin = 'round', arrow = grid::arrow(length = grid::unit(2, "mm"), type = "open")) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = baseline_frac), linewidth = 0.5, lineend = 'round', linejoin = 'round', linetype = "dashed") +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top") +
    ggplot2::facet_wrap(~term, ncol = 1) +
    ggplot2::geom_vline(xintercept = 0, lty = 2) +
    ggplot2::geom_text(
      ggplot2::aes(x = pmax(baseline_frac, new_frac), label = stringr::str_glue(" {round(relative_frac_diff * 100)}% (p = {signif(p.value, 2)})")),
      hjust = "left",
      nudge_x = 0.02
    ) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(), limits = c(0, 0.60), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(y = "Variable", x = "Absolute cell proportion", title = "Absolute change in cell type proportion per unit increase in variable")
}


plot_DCTC_model_coefs_forest <- function(results_tibble) {
  results_tibble |>
    dplyr::filter(!stringr::str_detect(term, "Intercept")) |>
    ggplot2::ggplot(ggplot2::aes(x = estimate, y = cluster, color = p.value < 0.05)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = estimate - std.error, xmax = estimate + std.error)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::facet_wrap(~term, scales = "free_x") +
    ggplot2::scale_x_continuous(limits = symmetric_limits)
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

  missing_IDs_idx <- which(gene_ID_vec %!in% names(named_output_vec))
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

#' Read cellbender h5 matrix
#'
#' Read a Cell Ranger HDF5 count matrix into a sparse gene-by-cell matrix.
#'
#' @param cellbender_h5_file Path to a CellBender HDF5 output file with a 10x-style `matrix` group.
#' @param feature_type 10x feature type to retain from the HDF5 file, usually `Gene Expression`.
#' @return A named matrix-like object with rows and columns aligned to the input feature/cell identifiers.
#' @keywords internal

read_cellbender_h5_matrix <- function(cellbender_h5_file, feature_type = "Gene Expression") {
  cellbender_h5_file_con <- hdf5r::H5File$new(cellbender_h5_file, mode = "r")
  on.exit(cellbender_h5_file_con$close_all(), add = TRUE)

  matrix_group <- cellbender_h5_file_con[["matrix"]]
  matrix_shape <- as.integer(matrix_group[["shape"]][])
  matrix_data <- as.numeric(matrix_group[["data"]][])
  feature_index <- as.integer(matrix_group[["indices"]][]) + 1L
  column_pointer <- as.integer(matrix_group[["indptr"]][])
  column_lengths <- diff(column_pointer)
  barcode_index <- rep.int(seq_along(column_lengths), column_lengths)

  counts_matrix <- Matrix::sparseMatrix(
    i = feature_index,
    j = barcode_index,
    x = matrix_data,
    dims = matrix_shape
  )

  barcodes <- as.character(matrix_group[["barcodes"]][])
  feature_names <- as.character(matrix_group[["features/name"]][])
  feature_types <- as.character(matrix_group[["features/feature_type"]][])
  keep_features <- feature_types == feature_type

  if (!any(keep_features)) {
    stop("No features with feature_type '", feature_type, "' found in CellBender h5 file.")
  }

  counts_matrix <- counts_matrix[keep_features, , drop = FALSE]
  rownames(counts_matrix) <- make.unique(feature_names[keep_features])
  colnames(counts_matrix) <- barcodes
  methods::as(counts_matrix, "dgCMatrix")
}


#' Compute GEX residual PCA
#'
#' Compute Pearson-residual PCA for the metadata cells of a BPCells-backed GEX
#' count matrix, with residuals from Seurat SCTransform or BPCells.
#'
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param metadata_tibble Cell metadata with `barcode_w_prefix`; its cells present in the matrix are analyzed.
#' @param GEX_PCA_backend Backend used to construct the PCA input matrix (`Seurat_SCT` or `BPCells_native`).
#' @param SCT_regress_vars Optional metadata columns, or cell-cycle score columns, regressed from the residuals.
#' @param n_components Number of output dimensions/components to compute.
#' @param n_variable_features Number of high-residual-variance genes kept before PCA.
#' @param min_feature_count Minimum total counts required for a gene to be considered in GEX PCA.
#' @param threads Number of threads used for the PCA step.
#' @param return_normalized Whether to also return the normalized data for the
#'   Seurat export: the lazy BPCells residuals of the variable genes, or the
#'   Seurat `SCTAssay`, which then includes UMI-corrected counts, and the
#'   Seurat cell-cycle scores.
#' @return A list with cell embeddings, gene loadings, singular values, and
#'   variable-feature diagnostics from the residual PCA workflow. With
#'   `return_normalized`, a list of these `PCA_results` and `normalized`.
#' @keywords internal

run_GEX_PCA_BPCells <- function(
  GEX_counts_matrix,
  metadata_tibble,
  GEX_PCA_backend = c("Seurat_SCT", "BPCells_native"),
  SCT_regress_vars = NULL,
  n_components,
  n_variable_features = 3000,
  min_feature_count = 50,
  threads = 1,
  return_normalized = FALSE
) {
  GEX_PCA_backend <- match.arg(GEX_PCA_backend)
  metadata_tibble <- dplyr::distinct(metadata_tibble, .data$barcode_w_prefix, .keep_all = TRUE)
  barcodes <- intersect(metadata_tibble$barcode_w_prefix, colnames(GEX_counts_matrix))
  if (length(barcodes) == 0) {
    stop("No requested barcodes were found in the GEX count matrix.")
  }
  counts_matrix <- GEX_counts_matrix[, barcodes, drop = FALSE]
  feature_counts <- BPCells::rowSums(counts_matrix)
  keep_features <- names(feature_counts)[feature_counts > min_feature_count]
  if (length(keep_features) < 2) {
    stop("Too few expressed GEX features remain for PCA.")
  }
  counts_matrix <- counts_matrix[keep_features, , drop = FALSE]

  SCT_regress_vars <- setdiff(as.character(unlist(SCT_regress_vars)), c(NA, "", "NULL"))
  cell_attr <- as.data.frame(metadata_tibble[match(barcodes, metadata_tibble$barcode_w_prefix), ])
  rownames(cell_attr) <- barcodes
  if (any(SCT_regress_vars %in% c("S.Score", "G2M.Score", "Phase", "CC.Difference"))) {
    cell_attr <- add_cell_cycle_scores_to_cell_attr(counts_matrix, cell_attr)
  }
  missing_regress_vars <- setdiff(SCT_regress_vars, colnames(cell_attr))
  if (length(missing_regress_vars) > 0) {
    stop(
      "SCT regression variable(s) were requested but are not available in the aligned metadata: ",
      paste(missing_regress_vars, collapse = ", ")
    )
  }

  residuals <- if (identical(GEX_PCA_backend, "BPCells_native")) {
    get_BPCells_native_GEX_residuals(counts_matrix, cell_attr, SCT_regress_vars, n_variable_features, threads)
  } else {
    get_Seurat_SCT_GEX_residuals(counts_matrix, cell_attr, SCT_regress_vars, n_variable_features,
      keep_assay = return_normalized)
  }
  pearson_residuals <- residuals$pearson_residuals
  n_components <- min(as.integer(n_components), nrow(pearson_residuals) - 1L, ncol(pearson_residuals) - 1L)
  if (n_components < 1) {
    stop("Too few cells or features remain for PCA.")
  }

  if (identical(GEX_PCA_backend, "BPCells_native")) {
    svd <- BPCells::svds(pearson_residuals, k = n_components, threads = threads)
    singular_values <- svd$d
    feature_loadings <- svd$u
    cell_embeddings <- sweep(svd$v, 2, svd$d, FUN = "*")
  } else {
    # The dense residuals have few features, so an eigendecomposition of their
    # feature Gram matrix is exact and fast with threaded BLAS.
    previous_threads <- RhpcBLASctl::blas_get_num_procs()
    RhpcBLASctl::blas_set_num_threads(threads)
    on.exit(RhpcBLASctl::blas_set_num_threads(previous_threads), add = TRUE)
    eigen_out <- eigen(tcrossprod(pearson_residuals), symmetric = TRUE)
    keep_components <- seq_len(min(n_components, length(eigen_out$values)))
    singular_values <- sqrt(pmax(eigen_out$values[keep_components], 0))
    feature_loadings <- eigen_out$vectors[, keep_components, drop = FALSE]
    cell_embeddings <- crossprod(pearson_residuals, feature_loadings)
  }

  component_names <- paste0("PCA_", seq_along(singular_values))
  dimnames(cell_embeddings) <- list(colnames(pearson_residuals), component_names)
  dimnames(feature_loadings) <- list(rownames(pearson_residuals), component_names)
  variable_features <- rownames(pearson_residuals)
  PCA_results <- list(
    cell_embeddings = cell_embeddings,
    feature_loadings = feature_loadings,
    singular_values = singular_values,
    variable_features = variable_features,
    variable_feature_stats = tibble::tibble(
      gene = variable_features,
      residual_variance = unname(residuals$residual_variance[variable_features]),
      PCA_weighted_loading_strength = sqrt(rowSums(sweep(feature_loadings, 2, singular_values, "*")^2))
    )
  )
  if (!return_normalized) {
    return(PCA_results)
  }
  # Scored after the PCA so that it cannot change its value; regression may
  # already have added the scores.
  cell_cycle_cols <- c("S.Score", "G2M.Score", "Phase", "CC.Difference")
  if (!all(cell_cycle_cols %in% colnames(cell_attr))) {
    cell_attr <- tryCatch(add_cell_cycle_scores_to_cell_attr(counts_matrix, cell_attr), error = function(error) {
      warning("Cell-cycle scores are not exported: ", conditionMessage(error), call. = FALSE)
      cell_attr
    })
  }
  list(
    PCA_results = PCA_results,
    normalized = list(
      backend = GEX_PCA_backend,
      regressed_vars = SCT_regress_vars,
      scale_data = if (identical(GEX_PCA_backend, "BPCells_native")) pearson_residuals,
      SCT_assay = residuals$SCT_assay,
      cell_cycle_tibble = if (all(cell_cycle_cols %in% colnames(cell_attr))) {
        tibble::tibble(barcode_w_prefix = rownames(cell_attr), cell_attr[cell_cycle_cols])
      }
    )
  )
}

get_BPCells_native_GEX_residuals <- function(
  counts_matrix,
  cell_attr,
  SCT_regress_vars,
  n_variable_features,
  threads,
  min_var = 0,
  clip_range = c(-10, 10),
  min_theta = 1e-6,
  max_theta = 1e6
) {
  cell_read_counts <- BPCells::colSums(counts_matrix)
  row_stats <- BPCells::matrix_stats(counts_matrix, row_stats = "variance", threads = threads)$row_stats
  gene_mean <- row_stats["mean", ]
  gene_var <- row_stats["variance", ]
  gene_beta <- BPCells::rowSums(counts_matrix) / sum(cell_read_counts)
  gene_beta <- pmax(gene_beta, .Machine$double.eps)
  gene_theta <- ifelse(
    gene_var > gene_mean,
    gene_mean^2 / pmax(gene_var - gene_mean, .Machine$double.eps),
    max_theta
  )
  gene_theta <- pmin(pmax(gene_theta, min_theta), max_theta)

  pearson_residuals <- BPCells::sctransform_pearson(
    mat = counts_matrix,
    gene_theta = gene_theta,
    gene_beta = gene_beta,
    cell_read_counts = cell_read_counts,
    min_var = min_var,
    clip_range = clip_range
  )
  if (length(SCT_regress_vars) > 0) {
    pearson_residuals <- BPCells::regress_out(
      mat = pearson_residuals,
      latent_data = cell_attr[, SCT_regress_vars, drop = FALSE],
      prediction_axis = "row"
    )
  }

  residual_stats <- BPCells::matrix_stats(pearson_residuals, row_stats = "variance", threads = threads)$row_stats
  residual_variance <- residual_stats["variance", ]
  variable_features <- residual_variance |>
    sort(decreasing = TRUE) |>
    utils::head(n = min(n_variable_features, length(residual_variance))) |>
    names()
  list(
    pearson_residuals = pearson_residuals[variable_features, , drop = FALSE],
    residual_variance = residual_variance
  )
}

get_Seurat_SCT_GEX_residuals <- function(counts_matrix, cell_attr, SCT_regress_vars, n_variable_features,
                                         keep_assay = FALSE) {
  previous_options <- options(future.globals.maxSize = 40 * 1024^3)
  on.exit(options(previous_options), add = TRUE)

  sct_args <- list(
    object = SeuratObject::CreateAssay5Object(counts = counts_matrix),
    cell.attr = cell_attr,
    variable.features.n = n_variable_features,
    conserve.memory = TRUE,
    # Corrected counts and their log1p data layer are needed only for export.
    do.correct.umi = keep_assay
  )
  if (length(SCT_regress_vars) > 0) {
    sct_args$vars.to.regress <- SCT_regress_vars
    sct_args$latent.data <- cell_attr[, SCT_regress_vars, drop = FALSE]
  }
  sct_assay <- do.call(Seurat::SCTransform, sct_args)
  scale_data <- SeuratObject::GetAssayData(sct_assay, layer = "scale.data")
  variable_features <- intersect(SeuratObject::VariableFeatures(sct_assay), rownames(scale_data))
  if (!keep_assay) {
    rm(sct_assay)
    gc()
  }

  pearson_residuals <- scale_data[variable_features, , drop = FALSE]
  residual_variance <- matrixStats::rowVars(as.matrix(pearson_residuals))
  names(residual_variance) <- rownames(pearson_residuals)
  list(pearson_residuals = pearson_residuals, residual_variance = residual_variance,
    SCT_assay = if (keep_assay) sct_assay)
}

#' Add Seurat-compatible cell-cycle scores
#'
#' Score the Seurat 2019 S and G2M gene sets as `Seurat::CellCycleScoring()`
#' does, on BPCells log-normalized counts, and add `S.Score`, `G2M.Score`,
#' `Phase`, and `CC.Difference` to `cell_attr`. Gene symbols are matched case
#' insensitively, as for mouse.
#' @keywords internal
add_cell_cycle_scores_to_cell_attr <- function(counts_matrix, cell_attr, min_cell_cycle_features = 5) {
  features <- list(
    S.Score = Seurat::CaseMatch(Seurat::cc.genes.updated.2019$s.genes, rownames(counts_matrix)),
    G2M.Score = Seurat::CaseMatch(Seurat::cc.genes.updated.2019$g2m.genes, rownames(counts_matrix))
  )
  if (any(lengths(features) < min_cell_cycle_features)) {
    stop(
      "Cell-cycle regression was requested, but too few cell-cycle genes were found in the GEX matrix: ",
      "S-phase ", length(features$S.Score), " < ", min_cell_cycle_features, "; ",
      "G2M-phase ", length(features$G2M.Score), " < ", min_cell_cycle_features, "."
    )
  }

  cell_counts <- BPCells::colSums(counts_matrix)
  normalized_data <- counts_matrix |>
    BPCells::multiply_cols(ifelse(cell_counts > 0, 10000 / cell_counts, 0)) |>
    log1p()
  n_unique_vals <- min(vapply(
    features,
    \(feature_set) length(unique(BPCells::rowMeans(normalized_data[feature_set, , drop = FALSE]))),
    integer(1)
  ))
  if (n_unique_vals < 3) {
    stop("Cell-cycle scoring was requested, but the normalized GEX data has too few unique cell-cycle feature means.")
  }

  # Seurat::AddModuleScore(): control genes are drawn from the expression bin of each module gene.
  gene_means <- sort(BPCells::rowMeans(normalized_data))
  set.seed(1)
  gene_bins <- ggplot2::cut_number(
    gene_means + stats::rnorm(length(gene_means)) / 1e30,
    n = min(n_unique_vals - 2, 24),
    labels = FALSE,
    right = FALSE
  )
  names(gene_bins) <- names(gene_means)
  n_controls <- min(lengths(features))
  control_features <- lapply(features, function(feature_set) {
    unique(unlist(lapply(feature_set, function(feature) {
      bin_features <- names(gene_bins)[gene_bins == gene_bins[[feature]]]
      sample(bin_features, size = min(n_controls, length(bin_features)), replace = FALSE)
    }), use.names = FALSE))
  })
  scores <- vapply(
    names(features),
    \(score) BPCells::colMeans(normalized_data[features[[score]], , drop = FALSE]) -
      BPCells::colMeans(normalized_data[control_features[[score]], , drop = FALSE]),
    numeric(ncol(normalized_data))
  )

  phase <- apply(scores, 1, function(cell_scores) {
    if (all(cell_scores < 0)) {
      return("G1")
    }
    if (sum(cell_scores == max(cell_scores)) > 1) {
      return("Undecided")
    }
    c("S", "G2M")[which(cell_scores == max(cell_scores))]
  })
  cell_attr[colnames(normalized_data), "S.Score"] <- scores[, "S.Score"]
  cell_attr[colnames(normalized_data), "G2M.Score"] <- scores[, "G2M.Score"]
  cell_attr[colnames(normalized_data), "Phase"] <- phase
  cell_attr$CC.Difference <- cell_attr$S.Score - cell_attr$G2M.Score
  cell_attr
}

#' Prepare GEX metadata tibble
#'
#' Join PCA embeddings, UMAP coordinates, and clusters back onto GEX metadata.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param barcode_vec Optional barcode subset; cells not present in the matrix are ignored, and an empty intersection is an error.
#' @param donor_id_metadata_tibble Donor-level metadata tibble keyed by donor identifier.
#' @param GEM_well_metadata_tibble GEM well-level metadata tibble keyed by `GEM_well_ID`.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

prepare_GEX_metadata_tibble <- function(metadata_tibble, barcode_vec, donor_id_metadata_tibble = NULL, GEM_well_metadata_tibble = NULL) {
  metadata_out <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% barcode_vec) |>
    dplyr::arrange(match(.data$barcode_w_prefix, barcode_vec))

  if (!is.null(donor_id_metadata_tibble)) {
    metadata_out <- dplyr::left_join(metadata_out, donor_id_metadata_tibble, by = "donor_id")
  }
  if (!is.null(GEM_well_metadata_tibble)) {
    metadata_out <- dplyr::left_join(metadata_out, GEM_well_metadata_tibble, by = "GEM_well_ID")
  }

  metadata_out
}

#' Prepare scDblFinder GEM well tibble
#'
#' Build one scDblFinder branch record per 10x Genomics GEM well.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param cluster_collapse_list Optional named mapping used to collapse detailed
#'   cluster labels before passing them to scDblFinder.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param sample_col Metadata column identifying the GEM well/sample used as the
#'   dynamic-branch key.
#' @return A tibble with one row per GEM well and list-columns for barcodes and
#'   named collapsed cluster labels.
#' @keywords internal

prepare_scDblFinder_GEM_well_tibble <- function(metadata_tibble,
                                                cluster_collapse_list = NULL,
                                                cluster_col,
                                                sample_col = "GEM_well_ID") {
  metadata <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c("barcode_w_prefix", sample_col, cluster_col))) |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(!is.na(.data[[sample_col]]))

  metadata$scDblFinder_cluster <- collapse_cell_type_labels(
    metadata[[cluster_col]],
    cluster_collapse_list
  )

  metadata |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      barcode_vec = list(.data$barcode_w_prefix),
      cluster_vec = list(stats::setNames(.data$scDblFinder_cluster, .data$barcode_w_prefix)),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data[[sample_col]])
}

#' Run scDblFinder BPCells GEM well
#'
#' Run scDblFinder on one GEM well-specific slice of a BPCells feature matrix.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param scDblFinder_GEM_well_tibble One-row branch tibble with GEM well ID,
#'   barcode vector, and named cluster vector.
#' @param output_suffix Suffix appended to `scDblFinder.class_` and
#'   `scDblFinder.score_` output columns.
#' @param dbr.sd Doublet-rate uncertainty passed to `scDblFinder::scDblFinder()`.
#' @param sample_col Column containing the GEM well/sample identifier in the
#'   branch tibble.
#' @param ... Additional arguments passed to `scDblFinder::scDblFinder()`.
#' @return A tibble keyed by `barcode_w_prefix` with GEM well ID, doublet class,
#'   and doublet score columns for the requested output suffix.
#' @keywords internal

run_scDblFinder_BPCells_GEM_well <- function(feature_matrix,
                                             scDblFinder_GEM_well_tibble,
                                             output_suffix,
                                             dbr.sd = 1.0,
                                             sample_col = "GEM_well_ID",
                                             ...) {
  GEM_well_ID <- scDblFinder_GEM_well_tibble[[sample_col]][[1]]
  barcode_vec <- scDblFinder_GEM_well_tibble$barcode_vec[[1]]
  clusters <- scDblFinder_GEM_well_tibble$cluster_vec[[1]]
  class_col <- paste0("scDblFinder.class_", output_suffix)
  score_col <- paste0("scDblFinder.score_", output_suffix)

  barcode_vec <- intersect(barcode_vec, colnames(feature_matrix))
  if (length(barcode_vec) == 0) {
    stop("No barcodes found for scDblFinder GEM well ", GEM_well_ID)
  }

  # Keep scDblFinder memory bounded by materializing only one GEM well branch,
  # and avoid MulticoreParam forked copies inside the branch.
  counts_matrix <- feature_matrix[, barcode_vec, drop = FALSE] |>
    methods::as("dgCMatrix")
  clusters <- clusters[barcode_vec]

  scDblFinder::scDblFinder(
    counts_matrix,
    clusters = clusters,
    dbr.sd = dbr.sd,
    returnType = "scores",
    BPPARAM = BiocParallel::SerialParam(),
    ...
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("barcode_w_prefix") |>
    dplyr::transmute(
      barcode_w_prefix = .data$barcode_w_prefix,
      !!class_col := .data$class,
      !!score_col := .data$score,
      GEM_well_ID = GEM_well_ID
    )
}

#' Apply scDblFinder cell- and cluster-level filters
#'
#' Add scDblFinder annotations to cell metadata and optionally remove cells
#' called as doublets or all cells in clusters whose called-doublet fraction
#' exceeds a configured threshold. Cluster fractions are always calculated
#' from the original scDblFinder calls, independently of cell-level removal.
#'
#' @param metadata_tibble Cell metadata containing `barcode_w_prefix` and the
#'   configured cluster column.
#' @param scDblFinder_results_df scDblFinder results containing
#'   `barcode_w_prefix`, or row names that contain the prefixed barcodes.
#' @param class_col Name of the scDblFinder class column.
#' @param cluster_col Name of the cluster column in `metadata_tibble`.
#' @param remove_called_doublets Whether to remove individual cells called as
#'   doublets.
#' @param max_doublet_fraction_per_cluster Maximum allowed fraction of called
#'   doublets in a cluster. `NULL` disables whole-cluster removal.
#' @return Filtered cell metadata with scDblFinder columns joined by barcode.
#' @keywords internal

filter_metadata_by_scDblFinder <- function(
  metadata_tibble,
  scDblFinder_results_df,
  class_col,
  cluster_col,
  remove_called_doublets,
  max_doublet_fraction_per_cluster
) {
  if (!"barcode_w_prefix" %in% names(scDblFinder_results_df)) {
    scDblFinder_results_df <- scDblFinder_results_df |>
      tibble::rownames_to_column("barcode_w_prefix")
  }

  called_doublet_barcodes <- scDblFinder_results_df |>
    dplyr::filter(.data[[class_col]] == "doublet") |>
    dplyr::pull("barcode_w_prefix")

  excluded_called_doublet_barcodes <- if (remove_called_doublets) {
    called_doublet_barcodes
  } else {
    character()
  }

  excluded_high_doublet_cluster_barcodes <- if (
    is.null(max_doublet_fraction_per_cluster)
  ) {
    character()
  } else {
    high_doublet_clusters <- metadata_tibble |>
      dplyr::mutate(
        is_called_doublet = .data$barcode_w_prefix %in%
          called_doublet_barcodes
      ) |>
      dplyr::group_by(.data[[cluster_col]]) |>
      dplyr::summarise(
        doublet_fraction = mean(.data$is_called_doublet),
        .groups = "drop"
      ) |>
      dplyr::filter(
        .data$doublet_fraction > max_doublet_fraction_per_cluster
      )
    high_doublet_clusters <- high_doublet_clusters[[cluster_col]]

    metadata_tibble |>
      dplyr::filter(.data[[cluster_col]] %in% high_doublet_clusters) |>
      dplyr::pull("barcode_w_prefix")
  }

  excluded_barcodes <- base::union(
    excluded_called_doublet_barcodes,
    excluded_high_doublet_cluster_barcodes
  )

  metadata_tibble |>
    dplyr::filter(!.data$barcode_w_prefix %in% excluded_barcodes) |>
    dplyr::left_join(scDblFinder_results_df, by = "barcode_w_prefix")
}

#' Get feature groups from LSI loadings
#'
#' Cluster LSI loading profiles into feature groups for scDblFinder aggregation.
#'
#' @param LSI_loadings_tibble Feature loading tibble with `peak` and `LSI_<dim>`
#'   columns.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param n_groups Number of k-means groups to create; must not exceed the
#'   number of feature rows.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @return Named character vector mapping each peak to a
#'   `scDblFinder_feature_<group>` label.
#' @keywords internal

get_feature_groups_from_LSI_loadings <- function(LSI_loadings_tibble, dims, n_groups, seed = 1) {
  if (length(dims) == 0) {
    stop("At least one LSI dimension is required for feature grouping.")
  }
  if (nrow(LSI_loadings_tibble) < n_groups) {
    stop("n_groups cannot exceed the number of features in LSI_loadings_tibble.")
  }

  dim_cols <- paste0("LSI_", dims)
  missing_cols <- setdiff(c("peak", dim_cols), colnames(LSI_loadings_tibble))
  if (length(missing_cols) > 0) {
    stop("Missing columns in LSI_loadings_tibble: ", paste(missing_cols, collapse = ", "))
  }

  set.seed(seed)
  loadings_matrix <- LSI_loadings_tibble |>
    dplyr::select(dplyr::all_of(dim_cols)) |>
    as.matrix()

  groups <- stats::kmeans(loadings_matrix, centers = n_groups, iter.max = 50, nstart = 1)$cluster
  stats::setNames(paste0("scDblFinder_feature_", groups), LSI_loadings_tibble$peak)
}

aggregate_BPCells_rows_by_group <- function(feature_matrix, feature_groups, threads = 1) {
  feature_groups <- feature_groups[rownames(feature_matrix)]
  if (anyNA(feature_groups)) {
    stop("All feature_matrix row names must have a feature group.")
  }

  BPCells::pseudobulk_matrix(
    mat = t(feature_matrix),
    cell_groups = factor(feature_groups),
    method = "sum",
    threads = threads
  ) |>
    t()
}

calculate_BPCells_UCell_scores_from_matrix <- function(counts_matrix,
                                                       features,
                                                       max_rank = 1500,
                                                       chunk_size = 1000,
                                                       workers = 1,
                                                       w_neg = 1,
                                                       ties_method = "average",
                                                       missing_genes = c("impute", "skip")) {
  missing_genes <- match.arg(missing_genes)
  if (is.null(w_neg)) {
    w_neg <- 1
  }
  if (!is.numeric(w_neg) || w_neg < 0) {
    stop("Weight on negative signatures (w_neg) must be >= 0.")
  }
  if (!is.numeric(max_rank)) {
    stop("Rank cutoff (max_rank) must be numeric.")
  }
  if (!is.numeric(workers) || workers < 1) {
    stop("Number of workers must be >= 1.")
  }

  features <- normalize_UCell_signature_names(features)
  max_rank <- min(as.integer(max_rank), nrow(counts_matrix))
  workers <- as.integer(workers)
  if (any(lengths(features) > max_rank)) {
    stop("One or more signatures contain more genes than max_rank. Increase max_rank or use shorter signatures.")
  }

  feature_indices <- prepare_UCell_signature_indices(
    features = features,
    feature_names = rownames(counts_matrix),
    missing_genes = missing_genes
  )

  # BPCells-native reimplementation of UCell 2.14.0 scoring, kept separate from
  # the Seurat AddModuleScore-compatible helper above. We intentionally mirror
  # UCell's per-cell descending ranks, max-rank capping, signed positive/negative
  # signatures, and lower-bound clipping. Validation compares this helper against
  # UCell::ScoreSignatures_UCell() on generated matrices; numerical output should
  # match exactly apart from ordinary floating-point representation.
  chunks <- split(seq_len(ncol(counts_matrix)), ceiling(seq_len(ncol(counts_matrix)) / chunk_size))
  score_chunk <- function(chunk_idx) {
    counts_chunk <- as.matrix(counts_matrix[, chunk_idx, drop = FALSE])
    rank_chunk <- rank_UCell_count_chunk(counts_chunk, ties_method = ties_method)
    calculate_UCell_scores_from_rank_chunk(
      rank_chunk = rank_chunk,
      feature_indices = feature_indices,
      max_rank = max_rank,
      w_neg = w_neg
    )
  }

  chunk_scores <- if (workers == 1) {
    lapply(chunks, score_chunk)
  } else {
    parallel::mclapply(chunks, score_chunk, mc.cores = workers)
  }
  score_matrix <- do.call(rbind, chunk_scores)
  score_matrix <- score_matrix[colnames(counts_matrix), , drop = FALSE]

  as.data.frame(score_matrix, check.names = FALSE)
}

normalize_UCell_signature_names <- function(features) {
  default_names <- paste0("signature_", seq_along(features))
  if (is.null(names(features))) {
    names(features) <- default_names
    return(features)
  }

  invalid_names <- names(features) == "" | duplicated(names(features))
  names(features)[invalid_names] <- default_names[invalid_names]
  features
}

prepare_UCell_signature_indices <- function(features, feature_names, missing_genes) {
  lapply(features, function(signature) {
    negative_features <- grep("-$", unlist(signature), perl = TRUE, value = TRUE)
    positive_features <- setdiff(unlist(signature), negative_features)
    positive_features <- gsub("\\+$", "", positive_features, perl = TRUE)
    negative_features <- gsub("-$", "", negative_features, perl = TRUE)

    list(
      positive = get_UCell_feature_indices(feature_names, positive_features, missing_genes = missing_genes),
      negative = get_UCell_feature_indices(feature_names, negative_features, missing_genes = missing_genes)
    )
  })
}

get_UCell_feature_indices <- function(feature_names, signature, missing_genes) {
  idx <- match(signature, feature_names)
  if (identical(missing_genes, "skip")) {
    idx <- idx[!is.na(idx)]
  } else {
    idx[is.na(idx)] <- -1L
  }
  idx
}

rank_UCell_count_chunk <- function(counts_chunk, ties_method = "average") {
  rank_chunk <- matrixStats::colRanks(
    -counts_chunk,
    ties.method = ties_method,
    preserveShape = TRUE
  )
  dimnames(rank_chunk) <- dimnames(counts_chunk)
  rank_chunk
}

calculate_UCell_scores_from_rank_chunk <- function(rank_chunk, feature_indices, max_rank, w_neg) {
  scores <- vapply(feature_indices, function(signature_indices) {
    positive_score <- calculate_UCell_score_from_indices(rank_chunk, signature_indices$positive, max_rank)
    negative_score <- calculate_UCell_score_from_indices(rank_chunk, signature_indices$negative, max_rank)
    score <- positive_score - w_neg * negative_score
    score[score < 0] <- 0
    score
  }, numeric(ncol(rank_chunk)))

  if (is.vector(scores)) {
    scores <- matrix(scores, nrow = ncol(rank_chunk))
  }
  rownames(scores) <- colnames(rank_chunk)
  colnames(scores) <- names(feature_indices)
  scores
}

calculate_UCell_score_from_indices <- function(rank_chunk, feature_idx, max_rank) {
  signature_length <- length(feature_idx)
  if (signature_length == 0) {
    return(rep(0, ncol(rank_chunk)))
  }

  present_idx <- feature_idx[feature_idx > 0]
  missing_idx <- feature_idx[feature_idx < 0]
  rank_sum <- rep(length(missing_idx) * max_rank, ncol(rank_chunk))
  if (length(present_idx) > 0) {
    signature_ranks <- rank_chunk[present_idx, , drop = FALSE]
    signature_ranks[signature_ranks >= max_rank] <- max_rank
    rank_sum <- rank_sum + colSums(signature_ranks)
  }

  minimum_rank_sum <- signature_length * (signature_length + 1) / 2
  1 - (rank_sum - minimum_rank_sum) / (signature_length * max_rank - minimum_rank_sum)
}

#' Prepare one marker-comparison record per multicluster GEX cell type.
make_cluster_marker_groups <- function(metadata_tibble) {
  metadata_tibble |>
    dplyr::transmute(
      barcode_w_prefix,
      cluster = as.character(PCA_harmony_SNN_cluster),
      cell_type = as.character(PCA_harmony_SNN_cluster_cell_type)
    ) |>
    dplyr::filter(!is.na(cluster), !is.na(cell_type)) |>
    tidyr::nest(metadata = c(barcode_w_prefix, cluster)) |>
    dplyr::mutate(
      clusters = purrr::map(metadata, \(cells) stringr::str_sort(unique(cells$cluster), numeric = TRUE))
    ) |>
    dplyr::filter(lengths(clusters) > 1L) |>
    dplyr::arrange(cell_type)
}

#' Test clusters within one cell type; adjust across genes per contrast.
get_cluster_markers_from_matrix <- function(feature_matrix, group_record) {
  clusters <- group_record$clusters[[1]]
  markers <- get_BPCells_markers_from_matrix(
    feature_matrix = feature_matrix,
    metadata_tibble = group_record$metadata[[1]],
    group_col = "cluster"
  )
  if (length(clusters) == 2L) {
    markers <- dplyr::filter(markers, cluster == clusters[[1]])
  }
  markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(p_val_adj = stats::p.adjust(p_val_raw, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(cluster = factor(cluster, levels = clusters))
  list(cell_type = group_record$cell_type[[1]], clusters = clusters, markers = markers)
}

#' Plot embedding loadings from tibble
#'
#' Plot top positive and negative feature loadings for selected embedding dimensions.
#'
#' @param loadings_tibble Loading tibble containing feature names and one loading
#'   column per requested dimension.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param feature_col Column containing feature labels for the y axis.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param nfeatures Number of strongest positive and strongest negative features
#'   to label per dimension.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_embedding_loadings_from_tibble <- function(loadings_tibble, dims, feature_col = "gene", dim_prefix = "PCA_", nfeatures = 50) {
  dims |>
    purrr::set_names(paste0(dim_prefix, dims)) |>
    purrr::map(\(dim_idx) {
      dim_col <- paste0(dim_prefix, dim_idx)
      loadings_tibble |>
        dplyr::select(dplyr::all_of(c(feature_col, dim_col))) |>
        dplyr::slice_max(order_by = abs(.data[[dim_col]]), n = nfeatures) |>
        dplyr::mutate(!!feature_col := forcats::fct_reorder(.data[[feature_col]], .data[[dim_col]])) %>%
        ggplot2::ggplot(ggplot2::aes(x = .data[[dim_col]], y = .data[[feature_col]])) +
        ggplot2::geom_col() +
        ggplot2::labs(x = "Loading", y = NULL, title = dim_col)
    })
}

collapse_cell_type_labels <- function(cell_type_labels_vec, collapse_list) {
  if (is.null(collapse_list)) {
    return(cell_type_labels_vec)
  }
  # Build old_label -> new_label lookup. The previous imap+unlist approach produced
  # compound names ("NewName.OldName") that never matched, so lookups always returned NA.
  lookup_vec <- stats::setNames(
    rep(names(collapse_list), lengths(collapse_list)),
    unlist(collapse_list)
  )
  labels_chr <- as.character(cell_type_labels_vec)
  collapsed <- lookup_vec[labels_chr]
  ifelse(is.na(collapsed), labels_chr, collapsed) |>
    stats::setNames(names(cell_type_labels_vec))
}

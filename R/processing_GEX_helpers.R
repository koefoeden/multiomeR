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


#' Normalize GEX counts with Seurat SCT and run PCA
#'
#' Compute SCTransform residual PCA for a BPCells-backed GEX count matrix.
#'
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param barcode_vec Optional barcode subset; cells not present in the matrix are ignored, and an empty intersection is an error.
#' @param metadata_tibble Optional cell metadata. If supplied, it must contain `barcode_w_prefix` and is aligned to the count matrix cells.
#' @param organism_chr Organism identifier used for built-in cell-cycle gene sets (`Homo_sapiens` or `Mus_musculus`).
#' @param GEX_PCA_backend Backend used to construct the PCA input matrix (`Seurat_SCT` or `BPCells_native`).
#' @param SCT_regress_vars Optional metadata columns passed to `Seurat::SCTransform(vars.to.regress = ...)`.
#' @param n_components Number of output dimensions/components to compute.
#' @param n_variable_features Number of high-residual-variance genes kept before PCA.
#' @param min_feature_count Minimum total counts required for a gene to be considered in GEX PCA.
#' @param threads Number of BLAS threads used for the dense residual PCA step.
#' @return A list with cell embeddings, gene loadings, singular values, and
#'   variable-feature diagnostics from the residual PCA workflow.
#' @keywords internal

run_GEX_PCA_BPCells <- function(
  GEX_counts_matrix,
  barcode_vec = NULL,
  metadata_tibble = NULL,
  organism_chr = NULL,
  GEX_PCA_backend = "Seurat_SCT",
  SCT_regress_vars = NULL,
  n_components,
  n_variable_features = 3000,
  min_feature_count = 50,
  threads = 1
) {
  counts_matrix <- GEX_counts_matrix
  if (is.null(barcode_vec) && !is.null(metadata_tibble)) {
    if (!"barcode_w_prefix" %in% colnames(metadata_tibble)) {
      stop("metadata_tibble must contain a 'barcode_w_prefix' column.")
    }
    barcode_vec <- metadata_tibble$barcode_w_prefix
  }

  if (!is.null(barcode_vec)) {
    barcode_vec <- intersect(barcode_vec, colnames(counts_matrix))
    if (length(barcode_vec) == 0) {
      stop("No requested barcodes were found in the GEX count matrix.")
    }
    counts_matrix <- counts_matrix[, barcode_vec, drop = FALSE]
  }

  feature_counts <- BPCells::rowSums(counts_matrix)
  keep_features <- names(feature_counts)[feature_counts > min_feature_count]
  if (length(keep_features) < 2) {
    stop("Too few expressed GEX features remain for subgroup PCA.")
  }
  counts_matrix <- counts_matrix[keep_features, , drop = FALSE]

  SCT_regress_vars <- normalize_SCT_regress_vars(SCT_regress_vars)
  cell_attr <- prepare_SCT_cell_attr(metadata_tibble, colnames(counts_matrix))
  requested_cell_cycle_cols <- intersect(SCT_regress_vars, cell_cycle_score_cols())
  if (length(requested_cell_cycle_cols) > 0) {
    cell_attr <- add_cell_cycle_scores_to_cell_attr(
      counts_matrix = counts_matrix,
      cell_attr = cell_attr,
      organism_chr = organism_chr
    )
  }

  missing_regress_vars <- setdiff(SCT_regress_vars, colnames(cell_attr))
  if (length(missing_regress_vars) > 0) {
    stop(
      "SCT regression variable(s) were requested but are not available in the aligned metadata: ",
      paste(missing_regress_vars, collapse = ", ")
    )
  }

  GEX_PCA_backend <- match_GEX_PCA_backend(GEX_PCA_backend)
  if (identical(GEX_PCA_backend, "BPCells_native")) {
    return(run_BPCells_native_GEX_PCA(
      counts_matrix = counts_matrix,
      cell_attr = cell_attr,
      SCT_regress_vars = SCT_regress_vars,
      n_components = n_components,
      n_variable_features = n_variable_features,
      threads = threads
    ))
  }

  final_sct <- run_Seurat_SCT_for_PCA(
    counts_matrix = counts_matrix,
    cell_attr = cell_attr,
    SCT_regress_vars = SCT_regress_vars,
    n_variable_features = n_variable_features
  )
  gc()

  pearson_residuals <- final_sct$y
  if (is.null(pearson_residuals) || nrow(pearson_residuals) < 2 || ncol(pearson_residuals) < 2) {
    stop("SCTransform returned too few residuals for GEX PCA.")
  }
  variable_features <- final_sct$variable_features
  if (is.null(variable_features)) {
    variable_features <- rownames(pearson_residuals)
  }
  variable_features <- intersect(variable_features, rownames(pearson_residuals))
  pearson_residuals <- pearson_residuals[variable_features, , drop = FALSE]

  n_components <- min(as.integer(n_components), nrow(pearson_residuals) - 1L, ncol(pearson_residuals) - 1L)
  if (n_components < 1) {
    stop("Too few cells or features remain for subgroup PCA.")
  }

  pca_out <- run_dense_feature_gram_PCA(
    feature_by_cell_matrix = pearson_residuals,
    n_components = n_components,
    threads = threads
  )
  cell_embeddings <- pca_out$cell_embeddings
  rownames(cell_embeddings) <- colnames(pearson_residuals)
  colnames(cell_embeddings) <- paste0("PCA_", seq_len(ncol(cell_embeddings)))

  feature_loadings <- pca_out$feature_loadings
  rownames(feature_loadings) <- rownames(pearson_residuals)
  colnames(feature_loadings) <- paste0("PCA_", seq_len(ncol(feature_loadings)))
  variable_feature_variance <- get_SCT_variable_feature_variance(final_sct, pearson_residuals)
  variable_feature_stats <- tibble::tibble(
    gene = variable_features,
    residual_variance = unname(variable_feature_variance[variable_features]),
    PCA_weighted_loading_strength = sqrt(rowSums(sweep(feature_loadings, 2, pca_out$singular_values, "*")^2))
  )

  list(
    cell_embeddings = cell_embeddings,
    feature_loadings = feature_loadings,
    singular_values = pca_out$singular_values,
    variable_features = variable_features,
    variable_feature_stats = variable_feature_stats,
    GEX_PCA_backend = GEX_PCA_backend,
    SCT_regress_vars = SCT_regress_vars,
    cell_cycle_score_cols = intersect(cell_cycle_score_cols(), colnames(cell_attr))
  )
}

match_GEX_PCA_backend <- function(GEX_PCA_backend) {
  if (is.null(GEX_PCA_backend) || length(GEX_PCA_backend) == 0 || is.na(GEX_PCA_backend[[1]]) || !nzchar(GEX_PCA_backend[[1]])) {
    return("Seurat_SCT")
  }

  backend_chr <- as.character(GEX_PCA_backend[[1]])
  backend_key <- tolower(gsub("[^A-Za-z0-9]+", "_", backend_chr))
  backend_lookup <- c(
    seurat_sct = "Seurat_SCT",
    sct = "Seurat_SCT",
    bpcells_native = "BPCells_native",
    bpcells = "BPCells_native",
    native = "BPCells_native"
  )
  backend <- unname(backend_lookup[[backend_key]])
  if (is.null(backend)) {
    stop("Unsupported GEX_PCA_backend '", backend_chr, "'. Use 'Seurat_SCT' or 'BPCells_native'.")
  }
  backend
}

run_BPCells_native_GEX_PCA <- function(
  counts_matrix,
  cell_attr,
  SCT_regress_vars = NULL,
  n_components,
  n_variable_features = 3000,
  min_var = 0,
  clip_range = c(-10, 10),
  min_theta = 1e-6,
  max_theta = 1e6,
  threads = 1
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
  variable_feature_variance <- residual_stats["variance", ]
  variable_features <- variable_feature_variance |>
    sort(decreasing = TRUE) |>
    utils::head(n = min(n_variable_features, length(variable_feature_variance))) |>
    names()
  pearson_residuals <- pearson_residuals[variable_features, , drop = FALSE]

  n_components <- min(as.integer(n_components), nrow(pearson_residuals) - 1L, ncol(pearson_residuals) - 1L)
  if (n_components < 1) {
    stop("Too few cells or features remain for subgroup PCA.")
  }

  svd <- BPCells::svds(pearson_residuals, k = n_components, threads = threads)
  cell_embeddings <- sweep(svd$v, 2, svd$d, FUN = "*")
  rownames(cell_embeddings) <- colnames(pearson_residuals)
  colnames(cell_embeddings) <- paste0("PCA_", seq_len(ncol(cell_embeddings)))

  feature_loadings <- svd$u
  rownames(feature_loadings) <- rownames(pearson_residuals)
  colnames(feature_loadings) <- paste0("PCA_", seq_len(ncol(feature_loadings)))
  variable_feature_stats <- tibble::tibble(
    gene = variable_features,
    residual_variance = unname(variable_feature_variance[variable_features]),
    PCA_weighted_loading_strength = sqrt(rowSums(sweep(feature_loadings, 2, svd$d, "*")^2))
  )

  list(
    cell_embeddings = cell_embeddings,
    feature_loadings = feature_loadings,
    singular_values = svd$d,
    variable_features = variable_features,
    variable_feature_stats = variable_feature_stats,
    GEX_PCA_backend = "BPCells_native",
    SCT_regress_vars = SCT_regress_vars,
    cell_cycle_score_cols = intersect(cell_cycle_score_cols(), colnames(cell_attr))
  )
}

run_dense_feature_gram_PCA <- function(feature_by_cell_matrix, n_components, threads = 1) {
  with_blas_threads(threads, {
    gram_matrix <- tcrossprod(feature_by_cell_matrix)
    eigen_out <- eigen(gram_matrix, symmetric = TRUE)
    keep_components <- seq_len(min(n_components, length(eigen_out$values)))
    singular_values <- sqrt(pmax(eigen_out$values[keep_components], 0))
    feature_loadings <- eigen_out$vectors[, keep_components, drop = FALSE]
    cell_embeddings <- crossprod(feature_by_cell_matrix, feature_loadings)
  })

  list(
    cell_embeddings = cell_embeddings,
    feature_loadings = feature_loadings,
    singular_values = singular_values
  )
}

with_blas_threads <- function(threads, code) {
  threads <- as.integer(threads)
  if (!requireNamespace("RhpcBLASctl", quietly = TRUE) || is.na(threads) || threads < 1) {
    return(force(code))
  }

  old_threads <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(threads)
  on.exit(RhpcBLASctl::blas_set_num_threads(old_threads), add = TRUE)
  force(code)
}

normalize_SCT_regress_vars <- function(SCT_regress_vars) {
  if (is.null(SCT_regress_vars)) {
    return(character())
  }
  SCT_regress_vars <- unlist(SCT_regress_vars, use.names = FALSE)
  SCT_regress_vars <- as.character(SCT_regress_vars)
  SCT_regress_vars <- SCT_regress_vars[!is.na(SCT_regress_vars) & nzchar(SCT_regress_vars) & SCT_regress_vars != "NULL"]
  unique(SCT_regress_vars)
}

prepare_SCT_cell_attr <- function(metadata_tibble, barcode_vec) {
  if (is.null(metadata_tibble)) {
    return(data.frame(row.names = barcode_vec))
  }
  if (!"barcode_w_prefix" %in% colnames(metadata_tibble)) {
    stop("metadata_tibble must contain a 'barcode_w_prefix' column.")
  }

  metadata_out <- metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% barcode_vec) |>
    dplyr::arrange(match(.data$barcode_w_prefix, barcode_vec))

  if (!identical(metadata_out$barcode_w_prefix, barcode_vec)) {
    missing_barcodes <- setdiff(barcode_vec, metadata_out$barcode_w_prefix)
    stop("metadata_tibble is missing ", length(missing_barcodes), " requested barcode(s).")
  }

  metadata_out <- as.data.frame(metadata_out, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(metadata_out) <- metadata_out$barcode_w_prefix
  metadata_out
}

run_Seurat_SCT_for_PCA <- function(counts_matrix, cell_attr, SCT_regress_vars = NULL, n_variable_features = 3000) {
  sct_args <- list(
    object = SeuratObject::CreateAssay5Object(counts = counts_matrix),
    cell.attr = cell_attr,
    variable.features.n = n_variable_features,
    conserve.memory = TRUE,
    do.correct.umi = FALSE
  )
  if (length(SCT_regress_vars) > 0) {
    sct_args$vars.to.regress <- SCT_regress_vars
    sct_args$latent.data <- cell_attr[, SCT_regress_vars, drop = FALSE]
  }
  sct_assay <- do.call(Seurat::SCTransform, sct_args)

  list(
    y = SeuratObject::GetAssayData(sct_assay, layer = "scale.data"),
    variable_features = SeuratObject::VariableFeatures(sct_assay)
  )
}

cell_cycle_score_cols <- function() {
  c("S.Score", "G2M.Score", "Phase", "CC.Difference")
}

cell_cycle_gene_sets <- function(organism_chr) {
  s_genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes

  if (is.null(organism_chr) || organism_chr %in% c("Homo_sapiens", "hsapiens", "human", "Mus_musculus", "mmusculus", "mouse")) {
    return(list(s_genes = s_genes, g2m_genes = g2m_genes))
  }

  stop("Cell-cycle scoring was requested, but organism '", organism_chr, "' is not supported.")
}

add_cell_cycle_scores_to_cell_attr <- function(counts_matrix, cell_attr, organism_chr, min_cell_cycle_features = 5) {
  gene_sets <- cell_cycle_gene_sets(organism_chr)
  s_features <- Seurat::CaseMatch(gene_sets$s_genes, rownames(counts_matrix))
  g2m_features <- Seurat::CaseMatch(gene_sets$g2m_genes, rownames(counts_matrix))
  if (length(s_features) < min_cell_cycle_features || length(g2m_features) < min_cell_cycle_features) {
    stop(
      "Cell-cycle regression was requested, but too few cell-cycle genes were found in the GEX matrix: ",
      "S-phase ", length(s_features), " < ", min_cell_cycle_features, "; ",
      "G2M-phase ", length(g2m_features), " < ", min_cell_cycle_features, "."
    )
  }

  normalized_assay <- SeuratObject::CreateAssay5Object(counts = counts_matrix) |>
    Seurat::NormalizeData(verbose = FALSE)
  normalized_data <- SeuratObject::GetAssayData(normalized_assay, layer = "data")
  n_unique_vals <- min(
    length(unique(Matrix::rowMeans(normalized_data[s_features, , drop = FALSE]))),
    length(unique(Matrix::rowMeans(normalized_data[g2m_features, , drop = FALSE])))
  )
  if (n_unique_vals < 3) {
    stop("Cell-cycle scoring was requested, but the normalized GEX data has too few unique cell-cycle feature means.")
  }

  cc_scores <- calculate_Seurat_cell_cycle_scores_from_matrix(
    assay_data = normalized_assay,
    s.features = s_features,
    g2m.features = g2m_features,
    nbin = min(n_unique_vals - 2, 24)
  )

  score_cols <- c("S.Score", "G2M.Score", "Phase")
  cell_attr[rownames(cc_scores), score_cols] <- cc_scores[, score_cols, drop = FALSE]
  cell_attr$CC.Difference <- cell_attr$S.Score - cell_attr$G2M.Score
  cell_attr
}

calculate_Seurat_cell_cycle_scores_from_matrix <- function(assay_data, s.features, g2m.features, nbin, seed = 1) {
  features <- list(S.Score = s.features, G2M.Score = g2m.features)
  set.seed(seed)
  feature_scores <- Seurat::AddModuleScore(
    object = assay_data,
    features = features,
    ctrl = min(lengths(features)),
    nbin = nbin,
    name = "Cell.Cycle",
    seed = seed,
    slot = "data"
  )
  colnames(feature_scores) <- c("S.Score", "G2M.Score")

  phase <- apply(feature_scores, 1, function(scores) {
    if (all(scores < 0)) {
      return("G1")
    }
    if (sum(scores == max(scores)) > 1) {
      return("Undecided")
    }
    c("S", "G2M")[which(scores == max(scores))]
  })

  data.frame(
    S.Score = feature_scores[, "S.Score"],
    G2M.Score = feature_scores[, "G2M.Score"],
    Phase = phase,
    row.names = rownames(feature_scores)
  )
}

get_SCT_variable_feature_variance <- function(SCT_out, pearson_residuals) {
  if (!is.null(SCT_out$gene_attr) && "residual_variance" %in% colnames(SCT_out$gene_attr)) {
    residual_variance <- SCT_out$gene_attr$residual_variance
    names(residual_variance) <- rownames(SCT_out$gene_attr)
    return(residual_variance)
  }
  if (inherits(SCT_out, "SCTAssay")) {
    feature_attributes <- Seurat::SCTResults(SCT_out, slot = "feature.attributes")
    if (!is.null(feature_attributes) && "residual_variance" %in% colnames(feature_attributes)) {
      residual_variance <- feature_attributes$residual_variance
      names(residual_variance) <- rownames(feature_attributes)
      return(residual_variance)
    }
  }

  residual_variance <- matrixStats::rowVars(as.matrix(pearson_residuals))
  names(residual_variance) <- rownames(pearson_residuals)
  residual_variance
}

#' Prepare GEX metadata tibble
#'
#' Join PCA embeddings, UMAP coordinates, and clusters back onto GEX metadata.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param barcode_vec Optional barcode subset; cells not present in the matrix are ignored, and an empty intersection is an error.
#' @param donor_id_metadata_tibble Donor-level metadata tibble keyed by donor identifier.
#' @param reaction_ID_metadata_tibble Reaction-level metadata tibble keyed by 10x reaction ID.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

prepare_GEX_metadata_tibble <- function(metadata_tibble, barcode_vec, donor_id_metadata_tibble = NULL, reaction_ID_metadata_tibble = NULL) {
  metadata_out <- metadata_tibble |>
    dplyr::distinct(barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% barcode_vec) |>
    dplyr::arrange(match(.data$barcode_w_prefix, barcode_vec))

  if (!is.null(donor_id_metadata_tibble)) {
    metadata_out <- dplyr::left_join(metadata_out, donor_id_metadata_tibble, by = "donor_id")
  }
  if (!is.null(reaction_ID_metadata_tibble)) {
    metadata_out <- dplyr::left_join(metadata_out, reaction_ID_metadata_tibble, by = "TENX_reaction_ID")
  }

  metadata_out
}

#' Prepare scDblFinder reaction tibble
#'
#' Build one scDblFinder branch record per TENX reaction.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param cluster_collapse_list Optional named mapping used to collapse detailed
#'   cluster labels before passing them to scDblFinder.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param sample_col Metadata column identifying the reaction/sample used as the
#'   dynamic-branch key.
#' @return A tibble with one row per reaction and list-columns for barcodes and
#'   named collapsed cluster labels.
#' @keywords internal

prepare_scDblFinder_reaction_tibble <- function(metadata_tibble,
                                                cluster_collapse_list = NULL,
                                                cluster_col,
                                                sample_col = "TENX_reaction_ID") {
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

#' Run scDblFinder BPCells reaction
#'
#' Run scDblFinder on one reaction-specific slice of a BPCells feature matrix.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param scDblFinder_reaction_tibble One-row branch tibble with reaction ID,
#'   barcode vector, and named cluster vector.
#' @param output_suffix Suffix appended to `scDblFinder.class_` and
#'   `scDblFinder.score_` output columns.
#' @param dbr.sd Doublet-rate uncertainty passed to `scDblFinder::scDblFinder()`.
#' @param sample_col Column containing the reaction/sample identifier in the
#'   branch tibble.
#' @param ... Additional arguments passed to `scDblFinder::scDblFinder()`.
#' @return A tibble keyed by `barcode_w_prefix` with reaction ID, doublet class,
#'   and doublet score columns for the requested output suffix.
#' @keywords internal

run_scDblFinder_BPCells_reaction <- function(feature_matrix,
                                             scDblFinder_reaction_tibble,
                                             output_suffix,
                                             dbr.sd = 1.0,
                                             sample_col = "TENX_reaction_ID",
                                             ...) {
  reaction_id <- scDblFinder_reaction_tibble[[sample_col]][[1]]
  barcode_vec <- scDblFinder_reaction_tibble$barcode_vec[[1]]
  clusters <- scDblFinder_reaction_tibble$cluster_vec[[1]]
  class_col <- paste0("scDblFinder.class_", output_suffix)
  score_col <- paste0("scDblFinder.score_", output_suffix)

  barcode_vec <- intersect(barcode_vec, colnames(feature_matrix))
  if (length(barcode_vec) == 0) {
    stop("No barcodes found for scDblFinder reaction ", reaction_id)
  }

  # Keep scDblFinder memory bounded by materializing only one reaction branch,
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
      TENX_reaction_ID = reaction_id
    )
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

#' Add GEX module scores to metadata
#'
#' Add simple marker module scores from log-normalized GEX counts to metadata.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param named_marker_genes_list Named list of marker gene vectors. A trailing
#'   `-` marks genes whose average expression is subtracted from the module score.
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param scale_factor Scale factor used when normalizing counts, coverage, or marker scores.
#' @return Metadata tibble with one added numeric module-score column per marker
#'   list name.
#' @keywords internal

add_GEX_module_scores_to_metadata <- function(metadata_tibble, named_marker_genes_list, GEX_counts_matrix, scale_factor = 10000) {
  counts_matrix <- GEX_counts_matrix
  keep_barcodes <- intersect(metadata_tibble$barcode_w_prefix, colnames(counts_matrix))
  if (length(keep_barcodes) == 0) {
    stop("No metadata barcodes were found in the GEX count matrix.")
  }
  counts_matrix <- counts_matrix[, keep_barcodes, drop = FALSE]

  cell_counts <- BPCells::colSums(counts_matrix)
  col_scaling <- ifelse(cell_counts > 0, scale_factor / cell_counts, 0)
  log_norm_matrix <- counts_matrix |>
    BPCells::multiply_cols(col_scaling) |>
    BPCells::log1p_slow()

  score_tibble <- tibble::tibble(barcode_w_prefix = colnames(log_norm_matrix))
  for (module_name in names(named_marker_genes_list)) {
    marker_vec <- named_marker_genes_list[[module_name]]
    positive_features <- marker_vec[!stringr::str_detect(marker_vec, "-$")] |>
      stringr::str_remove_all("[+-]$")
    negative_features <- marker_vec[stringr::str_detect(marker_vec, "-$")] |>
      stringr::str_remove_all("[+-]$")

    positive_features <- intersect(unique(positive_features), rownames(log_norm_matrix))
    negative_features <- intersect(unique(negative_features), rownames(log_norm_matrix))

    positive_score <- if (length(positive_features) > 0) {
      BPCells::colMeans(log_norm_matrix[positive_features, , drop = FALSE])
    } else {
      rep(0, ncol(log_norm_matrix))
    }
    negative_score <- if (length(negative_features) > 0) {
      BPCells::colMeans(log_norm_matrix[negative_features, , drop = FALSE])
    } else {
      rep(0, ncol(log_norm_matrix))
    }

    score_tibble[[module_name]] <- positive_score - negative_score
  }

  metadata_tibble |>
    dplyr::left_join(score_tibble, by = "barcode_w_prefix")
}

#' Get BPCells markers within parent groups from matrix
#'
#' Run BPCells Wilcoxon marker testing within parent cluster groups.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @param parent_group_col Single column name used for parent group col; the column must exist in the relevant metadata tibble.
#' @param features Character vector of feature names to extract from the matrix row names; missing features are handled by the called helper.
#' @return A named matrix-like object with rows and columns aligned to the input feature/cell identifiers.
#' @keywords internal

get_BPCells_markers_within_parent_groups_from_matrix <- function(
  feature_matrix,
  metadata_tibble,
  group_col,
  parent_group_col,
  features = NULL
) {
  group_to_parent_tibble <- metadata_tibble |>
    dplyr::select(dplyr::all_of(c(group_col, parent_group_col))) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data[[parent_group_col]]) & !is.na(.data[[group_col]])) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(group_col, parent_group_col)), as.character))

  marker_metadata <- metadata_tibble |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(group_col, parent_group_col)), as.character)) |>
    dplyr::semi_join(group_to_parent_tibble, by = c(group_col, parent_group_col))

  reference_tibble <- group_to_parent_tibble |>
    dplyr::left_join(
      group_to_parent_tibble,
      by = parent_group_col,
      suffix = c("", "_reference"),
      relationship = "many-to-many"
    ) |>
    dplyr::filter(.data[[group_col]] != .data[[paste0(group_col, "_reference")]]) |>
    dplyr::group_by(.data[[group_col]], .data[[parent_group_col]]) |>
    dplyr::summarise(reference_children = stringr::str_c(.data[[paste0(group_col, "_reference")]], collapse = ", "), .groups = "drop")

  within_parent_groups_vec <- group_to_parent_tibble |>
    dplyr::count(.data[[parent_group_col]], name = "n_groups") |>
    dplyr::filter(.data$n_groups > 1) |>
    dplyr::pull(.data[[parent_group_col]])

  within_parent_marker_tibble <- within_parent_groups_vec |>
    purrr::map(\(parent_group_name) {
      get_BPCells_markers_from_matrix(
        feature_matrix = feature_matrix,
        metadata_tibble = dplyr::filter(marker_metadata, .data[[parent_group_col]] == parent_group_name),
        group_col = group_col,
        features = features
      )
    }) |>
    dplyr::bind_rows()

  if (nrow(within_parent_marker_tibble) > 0) {
    within_parent_marker_tibble <- within_parent_marker_tibble |>
      dplyr::left_join(reference_tibble, by = c("cluster" = group_col)) |>
      dplyr::mutate(
        child_cluster = .data$cluster,
        parent_cluster = .data[[parent_group_col]]
      )
  }

  singleton_groups_vec <- group_to_parent_tibble |>
    dplyr::add_count(.data[[parent_group_col]], name = "n_groups") |>
    dplyr::filter(.data$n_groups == 1) |>
    dplyr::pull(.data[[group_col]])

  singleton_marker_tibble <- singleton_groups_vec |>
    purrr::map(\(group_name) {
      get_BPCells_markers_from_matrix(
        feature_matrix = feature_matrix,
        metadata_tibble = dplyr::mutate(
          marker_metadata,
          marker_group = dplyr::if_else(.data[[group_col]] == group_name, group_name, "background")
        ),
        group_col = "marker_group",
        features = features
      ) |>
        dplyr::filter(.data$cluster == group_name) |>
        dplyr::mutate(
          child_cluster = group_name,
          parent_cluster = group_to_parent_tibble[[parent_group_col]][group_to_parent_tibble[[group_col]] == group_name],
          reference_children = "all"
        )
    }) |>
    dplyr::bind_rows()

  dplyr::bind_rows(within_parent_marker_tibble, singleton_marker_tibble)
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

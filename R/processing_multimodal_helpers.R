# Preprocessing, dimension reduction and clustering ----
normalize_embedding_rows <- function(embedding_matrix) {
  row_norms <- sqrt(rowSums(embedding_matrix^2))
  row_norms[!is.finite(row_norms) | row_norms == 0] <- 1
  sweep(embedding_matrix, 1, row_norms, "/")
}

#' Average each cell's neighbours, given as index columns after the cell itself.
#' colMeans() over the neighbour dimension of a (neighbours x cells x dimensions)
#' array sums the same rows in the same order as colMeans() of one cell's rows.
impute_embedding_from_neighbor_indices <- function(embedding_matrix, neighbor_indices) {
  idx <- t(neighbor_indices[, -1, drop = FALSE])
  colMeans(array(embedding_matrix[idx, , drop = FALSE], c(nrow(idx), ncol(idx), ncol(embedding_matrix))), dims = 1)
}

row_euclidean_distance <- function(x, y) {
  sqrt(rowSums((x - y)^2))
}

select_embedding_dimensions <- function(embedding_matrix, dims, dim_prefix = NULL) {
  if (!is.null(dim_prefix)) {
    dim_cols <- paste0(dim_prefix, dims)
    if (all(dim_cols %in% colnames(embedding_matrix))) {
      return(embedding_matrix[, dim_cols, drop = FALSE])
    }
  }

  if (all(dims <= ncol(embedding_matrix))) {
    return(embedding_matrix[, dims, drop = FALSE])
  }

  stop("Requested embedding dimensions are not available.")
}

#' Find nearest neighbours with HNSW
#'
#' Search each cell's `k` nearest neighbours, the cell itself usually first. The
#' HNSW index is built on one thread, so the neighbours, and the clusters and
#' layouts built on them, are the same in every run; the search is threaded and
#' gives the same result for any thread count.
#'
#' @param data Numeric matrix with cells in rows; row names are kept.
#' @param k Number of neighbours per cell, including the cell itself.
#' @param metric `"cosine"` or `"euclidean"`.
#' @param threads Number of search threads.
#' @param ef HNSW search breadth; larger values improve recall at higher run time.
#' @return List with `idx` (one-based) and `dist` matrices.
#' @keywords internal
find_hnsw_neighbors <- function(data, k, metric = "cosine", threads = 1, ef = 500) {
  index <- RcppHNSW::hnsw_build(data, distance = metric, n_threads = 1, verbose = FALSE)
  knn <- RcppHNSW::hnsw_search(data, index, k = min(as.integer(k), nrow(data)), ef = ef, n_threads = threads, verbose = FALSE)
  rownames(knn$idx) <- rownames(knn$dist) <- rownames(data)
  knn
}

embedding_matrix_to_tibble <- function(embedding_matrix, cols = colnames(embedding_matrix)) {
  embedding_matrix[, cols, drop = FALSE] |>
    tibble::as_tibble(rownames = "barcode_w_prefix")
}

#' Get WNN embedding matrices
#'
#' Align and subset RNA and ATAC embeddings for weighted-nearest-neighbor analysis.
#'
#' @param GEX_embedding_matrix RNA/GEX cell embedding matrix with barcodes in
#'   row names.
#' @param ATAC_embedding_matrix ATAC cell embedding matrix with barcodes in row
#'   names.
#' @param barcode_vec Optional barcode subset; cells not present in the matrix are ignored, and an empty intersection is an error.
#' @param GEX_dims GEX embedding dimensions used when reconstructing reductions and neighbors in the export object.
#' @param ATAC_dims ATAC embedding dimensions used when reconstructing reductions and neighbors in the export object.
#' @return A named list with aligned `RNA` and `ATAC` embedding matrices,
#'   subset to shared barcodes and requested dimensions.
#' @keywords internal

get_WNN_embedding_matrices <- function(
  GEX_embedding_matrix,
  ATAC_embedding_matrix,
  barcode_vec,
  GEX_dims,
  ATAC_dims
) {
  shared_barcodes <- barcode_vec |>
    intersect(rownames(GEX_embedding_matrix)) |>
    intersect(rownames(ATAC_embedding_matrix))

  if (length(shared_barcodes) == 0) {
    stop("No shared barcodes were found across GEX embeddings, ATAC embeddings, and metadata.")
  }
  if (length(shared_barcodes) < length(barcode_vec)) {
    warning("Dropping ", length(barcode_vec) - length(shared_barcodes), " barcodes not present in both WNN embedding matrices.")
  }

  list(
    RNA = GEX_embedding_matrix[shared_barcodes, , drop = FALSE] |>
      select_embedding_dimensions(dims = GEX_dims, dim_prefix = "PCA_"),
    ATAC = ATAC_embedding_matrix[shared_barcodes, , drop = FALSE] |>
      select_embedding_dimensions(dims = ATAC_dims, dim_prefix = "LSI_")
  )
}

#' Calculate a small-SNN WNN kernel bandwidth
#'
#' Estimate each cell's kernel width from the farthest of its `k` lowest
#' nonzero shared-nearest-neighbour similarities, following Seurat's WNN
#' bandwidth strategy.
#'
#' @param embedding_matrix L2-normalized cell-by-dimension embedding matrix.
#' @param knn_idx Cell-by-neighbour matrix of one-based cell indices, including
#'   each query cell in its first column.
#' @param k Number of nearest neighbours used to construct the SNN and select
#'   low-similarity cells.
#' @param nearest_dist Distance to each cell's nearest non-self neighbour.
#' @param native_source_file Path to the tracked native C++ implementation.
#' @param threads Number of OpenMP threads; results do not depend on it.
#' @return Numeric bandwidth vector with one value per cell.
#' @keywords internal

calculate_small_SNN_bandwidth <- function(
  embedding_matrix,
  knn_idx,
  k,
  nearest_dist,
  native_source_file = file.path(
    get_project_root(),
    "src",
    "wnn_snn_bandwidth.cpp"
  ),
  threads = 1
) {
  embedding_matrix <- as.matrix(embedding_matrix)
  storage.mode(embedding_matrix) <- "double"
  knn_idx <- as.matrix(knn_idx)
  storage.mode(knn_idx) <- "integer"

  .Call(
    "multiomeR_wnn_small_snn_bandwidth",
    embedding_matrix,
    knn_idx,
    as.integer(k),
    as.numeric(nearest_dist),
    as.integer(threads),
    PACKAGE = load_WNN_native_library(native_source_file)
  )
}

# Floating-point contraction stays off so that the kernels reproduce R's rounding exactly.
load_WNN_native_library <- function(native_source_file) {
  load_native_library(native_source_file, "multiomeR_wnn", env = c(
    PKG_CXXFLAGS = "$(SHLIB_OPENMP_CXXFLAGS) -ffp-contract=off",
    PKG_LIBS = "$(SHLIB_OPENMP_CXXFLAGS)"
  ))
}

#' Find each modality's neighbours with the cell itself first
#'
#' HNSW occasionally misses a query cell among its own neighbours; it is then
#' put first and the farthest neighbour dropped. The small-SNN bandwidth uses
#' the first `k` columns.
find_WNN_modality_neighbors <- function(embeddings_norm, knn_k, k, threads, ef, native_source_file) {
  purrr::map(embeddings_norm, \(embedding) {
    nn <- find_hnsw_neighbors(embedding, k = knn_k, metric = "euclidean", threads = threads, ef = ef)
    for (cell in which(nn$idx[, 1] != seq_len(nrow(nn$idx)))) {
      keep <- nn$idx[cell, ] != cell
      nn$idx[cell, ] <- c(cell, nn$idx[cell, keep])[seq_len(knn_k)]
      nn$dist[cell, ] <- c(0, nn$dist[cell, keep])[seq_len(knn_k)]
    }
    nn$sigma <- pmax(calculate_small_SNN_bandwidth(embedding, nn$idx, k, nn$dist[, 2], native_source_file, threads),
      .Machine$double.eps)
    nn
  })
}

#' Seurat's WNN modality weights from within- versus cross-modality prediction
get_WNN_modality_weights <- function(embeddings_norm, modality_nn, k) {
  scores <- purrr::imap(embeddings_norm, \(embedding, modality) {
    own <- modality_nn[[modality]]
    kernel <- \(nn) {
      imputed <- impute_embedding_from_neighbor_indices(embedding, nn$idx[, seq_len(k + 1L), drop = FALSE])
      exp(-pmax(row_euclidean_distance(embedding, imputed) - own$dist[, 2], 0) / own$sigma)
    }
    within_kernel <- kernel(own)
    cross_scores <- purrr::map(modality_nn[names(modality_nn) != modality],
      \(nn) pmin(pmax(within_kernel / (kernel(nn) + 1e-4), 0), 200))
    rowSums(exp(do.call(cbind, cross_scores)))
  })
  total <- Reduce(`+`, scores)
  purrr::map(scores, \(score) score / total)
}

#' Weighted nearest neighbors BPCells
#'
#' Compute Seurat-style weighted nearest neighbors from aligned, L2-normalized
#' modality embeddings with Euclidean distances. In our tests, Seurat's modality
#' weights approached 0 or 1 as the number of cells grew, because the kernel
#' widths shrink with density, so they are estimated within random subsets of at
#' most `weight_cells` cells; neighbours are selected among all cells.
#' Adapted from Seurat WNN (Hao et al. 2021, doi:10.1016/j.cell.2021.04.048), without a Seurat object.
#'
#' @param embeddings_list Named list of modality embedding matrices with
#'   identical cell row names in the same order.
#' @param k Number of nearest neighbors to use for KNN/SNN construction.
#' @param candidate_k Number of candidate neighbors collected across modality
#'   KNN graphs before choosing the final weighted `k`.
#' @param threads Number of threads passed to BPCells, HNSW, or matrix-stat routines.
#' @param ef HNSW search breadth parameter; larger values improve recall at higher runtime/memory cost.
#' @param weight_cells Largest subset in which modality weights are estimated;
#'   `Inf` estimates them among all cells, as Seurat does.
#' @param seed Seed of the random assignment of cells to weight subsets.
#' @param extra_neighbors Neighbours returned beyond `k`, in score order, so that
#'   [restrict_WNN_results()] can later drop cells without a new search.
#' @param native_source_file Path to the tracked C++ implementation of the
#'   small-SNN kernel bandwidth.
#' @return A list containing modality weights per barcode and the weighted
#'   neighbor index and distance matrices.
#' @keywords internal

weighted_nearest_neighbors_BPCells <- function(
  embeddings_list,
  k = 20,
  candidate_k = 200,
  threads = 1,
  ef = 500,
  weight_cells = 1e4,
  seed = 1,
  extra_neighbors = 0,
  native_source_file = file.path(get_project_root(), "src", "wnn_snn_bandwidth.cpp")
) {
  cell_names <- rownames(embeddings_list[[1]])
  n_cells <- length(cell_names)
  k <- min(as.integer(k), n_cells - 1L)
  candidate_k <- min(max(as.integer(candidate_k), k), n_cells - 1L)
  knn_k <- min(max(k + 1L, candidate_k + 1L), n_cells)

  embeddings_norm <- purrr::map(embeddings_list, normalize_embedding_rows)
  modality_knn <- find_WNN_modality_neighbors(embeddings_norm, knn_k, k, threads, ef, native_source_file)
  nearest_dist <- purrr::map(modality_knn, \(nn) nn$dist[, 2])
  sigma_list <- purrr::map(modality_knn, "sigma")

  modality_weights <- if (n_cells <= weight_cells) {
    get_WNN_modality_weights(embeddings_norm, modality_knn, k)
  } else {
    subset <- withr::with_seed(seed, sample(rep_len(seq_len(ceiling(n_cells / weight_cells)), n_cells)))
    subset_weights <- lapply(split(seq_len(n_cells), subset), \(cells) {
      subset_embeddings <- purrr::map(embeddings_norm, \(embedding) embedding[cells, , drop = FALSE])
      get_WNN_modality_weights(subset_embeddings,
        find_WNN_modality_neighbors(subset_embeddings, k + 1L, k, threads, ef, native_source_file), k)
    })
    purrr::map(purrr::set_names(names(embeddings_norm)),
      \(modality) stats::setNames(unsplit(purrr::map(subset_weights, modality), subset), cell_names))
  }

  # Candidates are each modality's neighbours 2..candidate_k + 1 (first occurrence), scored by the
  # weighted kernel of every modality; the top k are kept in stable decreasing score order.
  selected <- .Call(
    "multiomeR_wnn_select",
    unname(purrr::map(embeddings_norm, base::t)),
    unname(purrr::map(modality_knn, "idx")),
    unname(nearest_dist),
    unname(sigma_list),
    unname(purrr::map(modality_weights, unname)),
    as.integer(candidate_k),
    as.integer(k + extra_neighbors),
    as.integer(threads),
    PACKAGE = load_WNN_native_library(native_source_file)
  )
  nn_idx <- selected[[1]]
  nn_dist <- selected[[2]]
  # The column names reproduce those of the former per-cell implementation.
  dimnames(nn_idx) <- list(cell_names, NULL)
  dimnames(nn_dist) <- list(cell_names, cell_names[nn_idx[1, ]])

  modality_weights_tibble <- tibble::tibble(barcode_w_prefix = cell_names)
  for (modality in names(modality_weights)) {
    modality_weights_tibble[[modality]] <- modality_weights[[modality]]
  }

  list(
    modality_weights = modality_weights_tibble,
    nn_idx = nn_idx,
    nn_dist = nn_dist
  )
}

#' Restrict WNN results to a subset of cells
#'
#' Keep each retained cell's first `k` neighbours among the retained cells, in
#' the original score order: the neighbours a selection among only the retained
#' candidates would give, with the same modality weights and kernel widths.
#'
#' @param WNN_results List returned by [weighted_nearest_neighbors_BPCells()],
#'   with at least `k` neighbours per cell.
#' @param barcodes Cells to keep.
#' @param k Number of neighbours to keep per cell.
#' @return WNN results for the retained cells with `k` neighbours each.
#' @keywords internal
restrict_WNN_results <- function(WNN_results, barcodes = rownames(WNN_results$nn_idx), k) {
  keep <- rownames(WNN_results$nn_idx) %in% barcodes
  new_index <- ifelse(keep, cumsum(keep), NA_integer_)
  idx <- matrix(new_index[WNN_results$nn_idx[keep, , drop = FALSE]], nrow = sum(keep))
  columns <- t(apply(!is.na(idx), 1, \(retained) which(retained)[seq_len(k)]))
  if (anyNA(columns)) {
    stop(sum(rowSums(is.na(columns)) > 0), " cells have fewer than ", k, " retained neighbours.")
  }
  positions <- cbind(rep(seq_len(nrow(idx)), k), as.vector(columns))
  list(
    modality_weights = WNN_results$modality_weights[keep, ],
    nn_idx = matrix(idx[positions], ncol = k, dimnames = list(rownames(WNN_results$nn_idx)[keep], NULL)),
    nn_dist = matrix(WNN_results$nn_dist[keep, , drop = FALSE][positions], ncol = k,
      dimnames = list(rownames(WNN_results$nn_idx)[keep], NULL))
  )
}

#' Run WNN UMAP
#'
#' Run UMAP from a precomputed weighted-nearest-neighbor graph.
#'
#' @param WNN_results List returned by WNN helpers, including KNN/SNN structures and modality weights aligned by barcode.
#' @param n_neighbors UMAP neighbor count; clipped below the number of input cells where the helper does that internally.
#' @param min_dist UMAP minimum-distance parameter controlling how tightly local neighborhoods are packed.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @param col_prefix Prefix assigned to generated coordinate columns, for example `LSI_UMAP` gives `LSI_UMAP_1` and `LSI_UMAP_2`.
#' @param n_components Number of UMAP output dimensions to compute.
#' @param threads Number of threads for the fuzzy graph and the optimization.
#' @return A tibble with `barcode_w_prefix` and generated UMAP coordinate
#'   columns named from `col_prefix`.
#' @keywords internal

run_WNN_UMAP <- function(WNN_results, n_neighbors, min_dist, seed = 1, col_prefix = "WNN_UMAP", n_components = 2, threads = 1) {
  run_UMAP_from_neighbors(
    knn = list(idx = WNN_results$nn_idx, dist = WNN_results$nn_dist),
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    seed = seed,
    col_prefix = col_prefix,
    n_components = n_components,
    threads = threads
  ) |>
    embedding_matrix_to_tibble()
}

#' Cluster WNN graph
#'
#' Cluster the WNN neighbor graph with Leiden and return metadata-ready labels.
#'
#' @param WNN_results List returned by WNN helpers, including KNN/SNN structures and modality weights aligned by barcode.
#' @param resolution Leiden clustering resolution; higher values generally split clusters more finely.
#' @param cluster_col Single metadata column name used as the cluster/grouping variable.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @param min_barcodes Minimum number of barcodes required for a cluster/group to be retained; smaller groups are dropped.
#' @return A tibble with `barcode_w_prefix` and one cluster column named by
#'   `cluster_col`; clusters below `min_barcodes` are dropped.
#' @keywords internal

cluster_WNN_graph <- function(WNN_results, resolution, cluster_col = "WNN_harmony_SNN_cluster", seed = 1, min_barcodes = 100) {
  clusters <- cluster_neighbors(
    knn = list(idx = WNN_results$nn_idx, dist = WNN_results$nn_dist),
    k = ncol(WNN_results$nn_idx),
    resolution = resolution,
    seed = seed,
    min_barcodes = min_barcodes
  )

  tibble::tibble(
    barcode_w_prefix = names(clusters),
    !!cluster_col := clusters
  )
}

get_SNN_matrix_from_knn <- function(knn, cell_names = rownames(knn$idx)) {
  snn <- BPCells::knn_to_snn_graph(knn, return_type = "list")
  snn_matrix <- Matrix::sparseMatrix(
    i = snn$i + 1L,
    j = snn$j + 1L,
    x = snn$weight,
    dims = c(snn$dim, snn$dim),
    dimnames = list(cell_names, cell_names)
  )
  Matrix::drop0(snn_matrix + Matrix::t(snn_matrix))
}

get_SNN_matrix_from_WNN_results <- function(WNN_results) {
  get_SNN_matrix_from_knn(
    knn = list(idx = WNN_results$nn_idx, dist = WNN_results$nn_dist),
    cell_names = rownames(WNN_results$nn_idx)
  )
}

#' Build WNN metadata tibble
#'
#' Join WNN weights, UMAP coordinates, and clusters back onto cell metadata.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param WNN_results List returned by WNN helpers, including KNN/SNN structures and modality weights aligned by barcode.
#' @param UMAP_embeddings_tibble WNN UMAP coordinate tibble keyed by
#'   `barcode_w_prefix`.
#' @param WNN_clusters_tibble WNN cluster-label tibble keyed by `barcode_w_prefix`.
#' @param GEX_UMAP_embeddings_tibble Optional RNA-only UMAP coordinate tibble to
#'   append for comparison plots.
#' @return Metadata restricted to cells present in WNN outputs, with old generated
#'   WNN/GEX columns removed before fresh outputs are joined.
#' @keywords internal

build_WNN_metadata_tibble <- function(metadata_tibble, WNN_results, UMAP_embeddings_tibble, WNN_clusters_tibble, GEX_UMAP_embeddings_tibble = NULL) {
  modality_weights_tibble <- WNN_results$modality_weights |>
    dplyr::rename_with(\(modality) paste0(modality, ".weight"), -barcode_w_prefix)

  generated_tibbles <- list(modality_weights_tibble, UMAP_embeddings_tibble, WNN_clusters_tibble)
  if (!is.null(GEX_UMAP_embeddings_tibble)) {
    generated_tibbles <- append(generated_tibbles, list(GEX_UMAP_embeddings_tibble))
  }

  # Rebuilding WNN metadata can start from metadata that already has older
  # generated WNN/GEX columns. Drop them before joining to avoid .x/.y suffixes.
  generated_cols <- generated_tibbles |>
    purrr::map(colnames) |>
    unlist(use.names = FALSE) |>
    setdiff("barcode_w_prefix")

  metadata_out <- metadata_tibble |>
    dplyr::semi_join(modality_weights_tibble, by = "barcode_w_prefix") |>
    dplyr::semi_join(WNN_clusters_tibble, by = "barcode_w_prefix") |>
    dplyr::select(-dplyr::any_of(generated_cols)) |>
    dplyr::left_join(modality_weights_tibble, by = "barcode_w_prefix") |>
    dplyr::left_join(UMAP_embeddings_tibble, by = "barcode_w_prefix") |>
    dplyr::left_join(WNN_clusters_tibble, by = "barcode_w_prefix")

  if (!is.null(GEX_UMAP_embeddings_tibble)) {
    metadata_out <- metadata_out |>
      dplyr::left_join(GEX_UMAP_embeddings_tibble, by = "barcode_w_prefix")
  }

  metadata_out
}

#' Scale each column of a count matrix to counts per 10,000.
normalize_BPCells_counts <- function(counts_matrix) {
  BPCells::multiply_cols(counts_matrix, 1e4 / pmax(BPCells::colSums(counts_matrix), 1))
}

#' Get BPCells markers from matrix
#'
#' Run BPCells Wilcoxon marker testing for selected features and metadata groups.
#'
#' @param feature_matrix Feature-by-cell matrix-like object with row names as feature IDs and column names as cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param group_col Single metadata column name used to group cells, samples, or features.
#' @param features Character vector of feature names to extract from the matrix row names; missing features are handled by the called helper.
#' @return A marker tibble from `BPCells::marker_features()` with Seurat-like
#'   `p_val`, `avg_log2FC`, `p_val_adj`, `cluster`, and `gene` columns added.
#' @keywords internal

get_BPCells_markers_from_matrix <- function(feature_matrix, metadata_tibble, group_col, features = NULL) {
  keep_barcodes <- intersect(colnames(feature_matrix), metadata_tibble$barcode_w_prefix)
  if (length(keep_barcodes) == 0) {
    stop("No overlapping barcodes between feature matrix and metadata.")
  }

  if (is.null(features)) {
    features <- rownames(feature_matrix)
  }
  features <- intersect(features, rownames(feature_matrix))
  if (length(features) == 0) {
    stop("No requested marker features were found in the feature matrix.")
  }

  if (!methods::is(feature_matrix, "IterableMatrix")) {
    feature_matrix <- feature_matrix |>
      Matrix::Matrix(sparse = TRUE) |>
      methods::as("dgCMatrix") |>
      BPCells::write_matrix_memory()
  }
  # Counts per 10,000 of each nucleus's total, so means and ranks compare
  # expression rather than sequencing depth.
  marker_mat <- normalize_BPCells_counts(feature_matrix[, keep_barcodes, drop = FALSE])[features, , drop = FALSE]

  # factor() drops levels without nuclei, such as marker sets that label no
  # cluster, on which BPCells::marker_features() fails.
  groups <- factor(stats::setNames(metadata_tibble[[group_col]], metadata_tibble$barcode_w_prefix)[keep_barcodes])
  BPCells::marker_features(marker_mat, groups, method = "wilcoxon") |>
    dplyr::mutate(
      p_val = p_val_raw,
      avg_log2FC = log2((foreground_mean + 1) / (background_mean + 1)),
      p_val_adj = stats::p.adjust(p_val_raw, method = "BH"),
      cluster = foreground,
      gene = feature
    )
}

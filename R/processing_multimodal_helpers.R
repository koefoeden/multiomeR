# Preprocessing, dimension reduction and clustering ----
normalize_embedding_rows <- function(embedding_matrix) {
  row_norms <- sqrt(rowSums(embedding_matrix^2))
  row_norms[!is.finite(row_norms) | row_norms == 0] <- 1
  sweep(embedding_matrix, 1, row_norms, "/")
}

impute_embedding_from_neighbor_indices <- function(embedding_matrix, neighbor_indices) {
  if (ncol(neighbor_indices) > 1 && all(neighbor_indices[, 1] == seq_len(nrow(neighbor_indices)))) {
    neighbor_indices <- neighbor_indices[, -1, drop = FALSE]
  }

  t(apply(
    neighbor_indices,
    1,
    \(idx) colMeans(embedding_matrix[idx, , drop = FALSE])
  ))
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
  )
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
    PACKAGE = load_native_library(native_source_file, "multiomeR_wnn")
  )
}

#' Weighted nearest neighbors BPCells
#'
#' Compute Seurat-style weighted nearest neighbors from aligned, L2-normalized
#' modality embeddings with Euclidean distances.
#'
#' @param embeddings_list Named list of modality embedding matrices with
#'   identical cell row names in the same order.
#' @param k Number of nearest neighbors to use for KNN/SNN construction.
#' @param candidate_k Number of candidate neighbors collected across modality
#'   KNN graphs before choosing the final weighted `k`.
#' @param threads Number of threads passed to BPCells, HNSW, or matrix-stat routines.
#' @param ef HNSW search breadth parameter; larger values improve recall at higher runtime/memory cost.
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
  native_source_file = file.path(get_project_root(), "src", "wnn_snn_bandwidth.cpp")
) {
  cell_names <- rownames(embeddings_list[[1]])
  n_cells <- length(cell_names)
  k <- min(as.integer(k), n_cells - 1L)
  candidate_k <- min(max(as.integer(candidate_k), k), n_cells - 1L)
  knn_k <- min(max(k + 1L, candidate_k + 1L), n_cells)

  embeddings_norm <- purrr::map(embeddings_list, normalize_embedding_rows)

  modality_knn <- purrr::map(embeddings_norm, \(embedding) {
    BPCells::knn_hnsw(
      data = embedding,
      k = knn_k,
      metric = "euclidean",
      verbose = FALSE,
      threads = threads,
      ef = ef
    )
  })

  nearest_dist <- purrr::map(modality_knn, \(nn) nn$dist[, 2])
  sigma_list <- purrr::map2(
    embeddings_norm,
    modality_knn,
    \(embedding, nn) {
      bandwidth <- calculate_small_SNN_bandwidth(
        embedding_matrix = embedding,
        knn_idx = nn$idx,
        k = k,
        nearest_dist = nn$dist[, 2],
        native_source_file = native_source_file
      )
      pmax(bandwidth, .Machine$double.eps)
    }
  )

  modality_scores <- purrr::map(names(embeddings_norm), \(modality) {
    own_embedding <- embeddings_norm[[modality]]
    own_nn <- modality_knn[[modality]]
    within_imputed <- impute_embedding_from_neighbor_indices(own_embedding, own_nn$idx[, seq_len(k + 1L), drop = FALSE])
    within_dist <- pmax(row_euclidean_distance(own_embedding, within_imputed) - nearest_dist[[modality]], 0)
    within_kernel <- exp(-1 * (within_dist / sigma_list[[modality]]))

    cross_scores <- purrr::map(setdiff(names(embeddings_norm), modality), \(other_modality) {
      other_nn <- modality_knn[[other_modality]]
      cross_imputed <- impute_embedding_from_neighbor_indices(own_embedding, other_nn$idx[, seq_len(k + 1L), drop = FALSE])
      cross_dist <- pmax(row_euclidean_distance(own_embedding, cross_imputed) - nearest_dist[[modality]], 0)
      cross_kernel <- exp(-1 * (cross_dist / sigma_list[[modality]]))
      pmin(pmax(within_kernel / (cross_kernel + 1e-4), 0), 200)
    })

    rowSums(exp(do.call(cbind, cross_scores)))
  }) |>
    purrr::set_names(names(embeddings_norm))

  total_modality_score <- Reduce(`+`, modality_scores)
  modality_weights <- purrr::map(modality_scores, \(score) score / total_modality_score)

  candidate_indices <- lapply(seq_len(n_cells), \(cell_idx) {
    Reduce(
      union,
      purrr::map(modality_knn, \(nn) nn$idx[cell_idx, seq_len(candidate_k + 1L)][-1])
    )
  })

  candidate_distances <- purrr::imap(embeddings_norm, \(embedding, modality) {
    lapply(seq_len(n_cells), \(cell_idx) {
      neighbors <- candidate_indices[[cell_idx]]
      query <- matrix(embedding[cell_idx, ], nrow = length(neighbors), ncol = ncol(embedding), byrow = TRUE)
      pmax(row_euclidean_distance(query, embedding[neighbors, , drop = FALSE]) - nearest_dist[[modality]][[cell_idx]], 0)
    })
  })

  weighted_scores <- lapply(seq_len(n_cells), \(cell_idx) {
    Reduce(
      `+`,
      purrr::map(names(embeddings_norm), \(modality) {
        exp(-1 * (candidate_distances[[modality]][[cell_idx]] / sigma_list[[modality]][[cell_idx]])) *
          modality_weights[[modality]][[cell_idx]]
      })
    )
  })

  selected_order <- purrr::map(weighted_scores, \(score) order(score, decreasing = TRUE))
  nn_idx <- t(vapply(
    seq_len(n_cells),
    \(cell_idx) {
      candidate_indices[[cell_idx]][selected_order[[cell_idx]][seq_len(k)]]
    },
    integer(k)
  ))
  nn_dist <- t(vapply(
    seq_len(n_cells),
    \(cell_idx) {
      score <- weighted_scores[[cell_idx]][selected_order[[cell_idx]][seq_len(k)]]
      sqrt(pmin(pmax((1 - score) / 2, 0), 1))
    },
    numeric(k)
  ))

  rownames(nn_idx) <- rownames(nn_dist) <- cell_names

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
#' @param threads Number of threads used for nearest-neighbor processing.
#' @return A tibble with `barcode_w_prefix` and generated UMAP coordinate
#'   columns named from `col_prefix`.
#' @keywords internal

run_WNN_UMAP <- function(WNN_results, n_neighbors, min_dist, seed = 1, col_prefix = "WNN_UMAP", n_components = 2, threads = 1) {
  n_neighbors <- min(as.integer(n_neighbors), ncol(WNN_results$nn_idx))
  set.seed(seed)
  umap <- uwot::umap(
    X = NULL,
    nn_method = list(
      idx = WNN_results$nn_idx[, seq_len(n_neighbors), drop = FALSE],
      dist = WNN_results$nn_dist[, seq_len(n_neighbors), drop = FALSE]
    ),
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    n_components = n_components,
    n_threads = threads,
    n_sgd_threads = 0,
    verbose = TRUE
  )

  rownames(umap) <- rownames(WNN_results$nn_idx)
  colnames(umap) <- paste0(col_prefix, "_", seq_len(ncol(umap)))
  embedding_matrix_to_tibble(umap)
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
  clusters <- cluster_knn_snn_leiden(
    knn = list(idx = WNN_results$nn_idx, dist = WNN_results$nn_dist),
    resolution = resolution,
    seed = seed
  )
  names(clusters) <- rownames(WNN_results$nn_idx)
  clusters <- filter_clusters_by_min_barcodes(clusters, min_barcodes = min_barcodes)

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

  marker_mat <- feature_matrix[features, keep_barcodes, drop = FALSE]
  if (!methods::is(marker_mat, "IterableMatrix")) {
    marker_mat <- marker_mat |>
      Matrix::Matrix(sparse = TRUE) |>
      methods::as("dgCMatrix") |>
      BPCells::write_matrix_memory()
  }

  groups <- stats::setNames(metadata_tibble[[group_col]], metadata_tibble$barcode_w_prefix)[keep_barcodes]
  BPCells::marker_features(marker_mat, groups, method = "wilcoxon") |>
    dplyr::mutate(
      p_val = p_val_raw,
      avg_log2FC = log2((foreground_mean + 1e-9) / (background_mean + 1e-9)),
      p_val_adj = stats::p.adjust(p_val_raw, method = "BH"),
      cluster = foreground,
      gene = feature
    )
}

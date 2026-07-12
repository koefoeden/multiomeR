get_seurat_export_cells <- function(metadata_tibble) {
  cells <- metadata_tibble$barcode_w_prefix
  if (is.null(cells) || anyNA(cells) || any(cells == "")) {
    stop("metadata_tibble must contain non-empty barcode_w_prefix values.")
  }
  if (anyDuplicated(cells)) {
    stop("metadata_tibble contains duplicated barcode_w_prefix values.")
  }
  cells
}

check_seurat_export_cells <- function(cells, named_cell_sets) {
  purrr::iwalk(named_cell_sets, \(available_cells, source_name) {
    missing_cells <- base::setdiff(cells, available_cells)
    if (length(missing_cells) > 0) {
      stop(
        "Missing ", length(missing_cells), " final WNN cell(s) from ",
        source_name, ". First missing: ",
        paste(utils::head(missing_cells, 20), collapse = ", ")
      )
    }
  })
  invisible(cells)
}

#' Prepare seurat export metadata
#'
#' Build the legacy Seurat/Signac export object from BPCells-native matrices and metadata.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param WNN_results List returned by WNN helpers, including KNN/SNN structures and modality weights aligned by barcode.
#' @param cells Character vector of cell barcodes defining the order and subset to export or align.
#' @return A Seurat/Signac object or component with metadata, assays, reductions, and graphs aligned by cell barcode.
#' @keywords internal

prepare_seurat_export_metadata <- function(metadata_tibble, cells, WNN_results = NULL) {
  metadata_df <- metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% cells) |>
    dplyr::arrange(match(.data$barcode_w_prefix, cells)) |>
    dplyr::select(-dplyr::matches("(_UMAP_[12]$)")) |>
    dplyr::select(-dplyr::any_of(c("RNA.weight", "SCTregr.weight", "ATAC.weight"))) |>
    base::as.data.frame()

  if (!is.null(WNN_results)) {
    modality_weights <- WNN_results$modality_weights |>
      dplyr::filter(.data$barcode_w_prefix %in% cells) |>
      dplyr::arrange(match(.data$barcode_w_prefix, cells))

    gex_weight_col <- base::setdiff(colnames(modality_weights), c("barcode_w_prefix", "ATAC"))[[1]]
    modality_weights <- modality_weights |>
      dplyr::transmute(
        barcode_w_prefix = .data$barcode_w_prefix,
        SCTregr.weight = .data[[gex_weight_col]],
        ATAC.weight = .data$ATAC
      ) |>
      base::as.data.frame()

    metadata_df <- dplyr::left_join(metadata_df, modality_weights, by = "barcode_w_prefix")
  }

  rownames(metadata_df) <- metadata_df$barcode_w_prefix
  metadata_df
}

align_matrix_columns <- function(matrix, cells, matrix_name) {
  missing_cells <- base::setdiff(cells, colnames(matrix))
  if (length(missing_cells) > 0) {
    stop(matrix_name, " is missing ", length(missing_cells), " requested cell(s).")
  }
  matrix[, cells, drop = FALSE]
}

align_matrix_rows <- function(matrix, features, matrix_name) {
  missing_features <- base::setdiff(features, rownames(matrix))
  if (length(missing_features) > 0) {
    stop(matrix_name, " is missing ", length(missing_features), " requested feature(s).")
  }
  matrix[features, , drop = FALSE]
}

add_feature_metadata <- function(assay, feature_metadata_df) {
  if (is.null(feature_metadata_df) || nrow(feature_metadata_df) == 0) {
    return(assay)
  }

  common_features <- base::intersect(rownames(assay), rownames(feature_metadata_df))
  if (length(common_features) == 0) {
    return(assay)
  }

  metadata_out <- feature_metadata_df[common_features, , drop = FALSE]
  SeuratObject::AddMetaData(assay, metadata = metadata_out)
}

get_RNA_feature_metadata <- function(gene_features_df, RNA_features) {
  if (is.null(gene_features_df) || nrow(gene_features_df) == 0) {
    return(NULL)
  }

  feature_metadata <- gene_features_df |>
    dplyr::mutate(gene_name_Seurat = rownames(gene_features_df)) |>
    dplyr::select(dplyr::any_of(c("gene_name_Seurat", "seqnames", "start", "end", "id", "name"))) |>
    base::as.data.frame()
  rownames(feature_metadata) <- rownames(gene_features_df)
  feature_metadata[base::intersect(RNA_features, rownames(feature_metadata)), , drop = FALSE]
}

align_peak_GRanges_to_matrix <- function(peak_GRanges, peak_names) {
  available_peak_names <- get_peak_names_from_GRanges(peak_GRanges)
  peak_idx <- match(peak_names, available_peak_names)
  if (anyNA(peak_idx)) {
    stop("Some ATAC peak matrix rows were not found in ATAC_peak_GRanges.")
  }

  peak_GRanges <- peak_GRanges[peak_idx]
  names(peak_GRanges) <- peak_names
  peak_GRanges
}

get_peak_feature_metadata <- function(annotated_peak_GRanges, peak_names) {
  if (is.null(annotated_peak_GRanges)) {
    return(NULL)
  }

  annotated_peak_GRanges <- align_peak_GRanges_to_matrix(annotated_peak_GRanges, peak_names)
  metadata_df <- GenomicRanges::as.data.frame(annotated_peak_GRanges) |>
    dplyr::select(-dplyr::any_of(c("seqnames", "start", "end", "width", "strand"))) |>
    base::as.data.frame()
  rownames(metadata_df) <- peak_names
  metadata_df
}

#' Ensure signac annotation GRanges
#'
#' Build the legacy Seurat/Signac export object from BPCells-native matrices and metadata.
#'
#' @param annotation_GRanges Gene annotation GRanges; missing Signac-required metadata columns are added with safe defaults.
#' @return A GRanges object with coordinates and metadata columns expected by downstream ATAC helpers.
#' @keywords internal

ensure_signac_annotation_GRanges <- function(annotation_GRanges) {
  if (is.null(annotation_GRanges)) {
    annotation_GRanges <- GenomicRanges::GRanges()
  }

  if (!"gene_name" %in% colnames(S4Vectors::mcols(annotation_GRanges))) {
    annotation_GRanges$gene_name <- names(annotation_GRanges) %||% character(length(annotation_GRanges))
  }
  if (!"gene_id" %in% colnames(S4Vectors::mcols(annotation_GRanges))) {
    annotation_GRanges$gene_id <- names(annotation_GRanges) %||% character(length(annotation_GRanges))
  }
  if (!"gene_biotype" %in% colnames(S4Vectors::mcols(annotation_GRanges))) {
    annotation_GRanges$gene_biotype <- rep("unknown", length(annotation_GRanges))
  }
  if (!"tx_id" %in% colnames(S4Vectors::mcols(annotation_GRanges)) && !"transcript_id" %in% colnames(S4Vectors::mcols(annotation_GRanges))) {
    annotation_GRanges$tx_id <- names(annotation_GRanges) %||% character(length(annotation_GRanges))
  }
  if (!"type" %in% colnames(S4Vectors::mcols(annotation_GRanges))) {
    annotation_GRanges$type <- rep("gene", length(annotation_GRanges))
  }
  annotation_GRanges
}

build_signac_annotation_GRanges <- function(reference_Ensembl_annotations_GRanges_list) {
  ensure_signac_annotation_GRanges(reference_Ensembl_annotations_GRanges_list$genes)
}

#' Make signac fragment records
#'
#' Build the legacy Seurat/Signac export object from BPCells-native matrices and metadata.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param GEM_well_ID_vec Character vector of `GEM_well_ID` values to include in the
#'   exported fragment records.
#' @param cellranger_summary_files Character vector/list of Cell Ranger summary
#'   file paths aligned to `GEM_well_ID_vec`.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

make_signac_fragment_records <- function(metadata_tibble, GEM_well_ID_vec, cellranger_summary_files) {
  cellranger_summary_files <- unlist(cellranger_summary_files, use.names = FALSE)

  if (length(cellranger_summary_files) != length(GEM_well_ID_vec)) {
    stop("cellranger_summary_files must have the same length as GEM_well_ID_vec.")
  }

  purrr::map2_dfr(
    GEM_well_ID_vec,
    seq_along(GEM_well_ID_vec),
    \(GEM_well_ID, idx) {
      fragment_file <- file.path(dirname(cellranger_summary_files[[idx]]), "atac_fragments.tsv.gz")
      assay_cell_names <- metadata_tibble$barcode_w_prefix[metadata_tibble$GEM_well_ID == GEM_well_ID]
      fragment_cell_names <- sub(paste0("^", GEM_well_ID, "_"), "", assay_cell_names)

      tibble::tibble(
        GEM_well_ID = GEM_well_ID,
        fragment_file = fragment_file,
        fragment_index_file = paste0(fragment_file, ".tbi"),
        assay_cell_names = list(assay_cell_names),
        fragment_cell_names = list(fragment_cell_names)
      )
    }
  )
}

create_signac_fragment_objects <- function(fragment_records_tibble) {
  if (is.null(fragment_records_tibble) || nrow(fragment_records_tibble) == 0) {
    return(list())
  }

  purrr::pmap(
    fragment_records_tibble,
    \(GEM_well_ID, fragment_file, fragment_index_file, assay_cell_names, fragment_cell_names) {
      cells <- fragment_cell_names
      names(cells) <- assay_cell_names
      Signac::CreateFragmentObject(
        path = fragment_file,
        index = fragment_index_file,
        cells = cells,
        validate.fragments = FALSE
      )
    }
  )
}

#' Add dimreduc from matrix
#'
#' Add an embedding matrix as a Seurat dimensional reduction.
#'
#' @param object Seurat object being modified in place and returned.
#' @param name Name of the Seurat dimensional reduction to create.
#' @param embeddings Cell embedding matrix with row names matching Seurat cells.
#' @param assay Seurat assay name associated with the reduction, graph, or metadata being added.
#' @param key Dimensional-reduction key prefix used by Seurat, conventionally ending in an underscore.
#' @param loadings Optional feature-loading matrix aligned to the embedding columns.
#' @param cells Character vector of cell barcodes defining the order and subset to export or align.
#' @return The Seurat object with the requested reduction added, or unchanged
#'   when `embeddings` is `NULL`.
#' @keywords internal

add_dimreduc_from_matrix <- function(object, name, embeddings, assay, key, loadings = NULL, cells = colnames(object)) {
  if (is.null(embeddings)) {
    return(object)
  }

  missing_cells <- base::setdiff(cells, rownames(embeddings))
  if (length(missing_cells) > 0) {
    stop("Reduction ", name, " is missing ", length(missing_cells), " cell(s).")
  }

  embeddings <- embeddings[cells, , drop = FALSE]
  if (is.null(loadings)) {
    loadings <- matrix(numeric(), nrow = 0, ncol = 0)
  }

  object[[name]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    loadings = loadings,
    assay = assay,
    key = key
  )
  object
}

embedding_tibble_to_matrix <- function(embedding_tibble) {
  if (is.null(embedding_tibble)) {
    return(NULL)
  }

  embedding_tibble |>
    base::as.data.frame() |>
    tibble::column_to_rownames("barcode_w_prefix") |>
    as.matrix()
}

first_or_null <- function(x) {
  if (length(x) == 0) {
    return(NULL)
  }
  x[[1]]
}

set_seurat_export_defaults <- function(object, preferred_cluster_col = NULL, preferred_reduction = NULL, preferred_graph = NULL, preferred_neighbor = NULL) {
  cluster_col <- preferred_cluster_col %||%
    first_or_null(base::intersect(
      c("WNN_harmony_SNN_cluster", "PCA_harmony_SNN_cluster_cell_type", "PCA_harmony_SNN_cluster"),
      colnames(object@meta.data)
    ))

  reduction <- preferred_reduction %||%
    first_or_null(base::intersect(
      c("WNN_harmony_NN_UMAP", "PCA_harmony_UMAP", "PCA_UMAP", "PCA_harmony", "PCA"),
      names(object@reductions)
    ))

  graph <- preferred_graph %||%
    first_or_null(base::intersect(c("WNN_harmony_SNN", "PCA_harmony_SNN"), names(object@graphs))) %||%
    NA_character_

  neighbor <- preferred_neighbor %||%
    first_or_null(base::intersect(c("WNN_harmony_NN"), names(object@neighbors))) %||%
    NA_character_

  object@misc$def_dim_reduc_full <- reduction
  object@misc$def_dim_reduc <- reduction
  object@misc$def_graph <- graph
  object@misc$def_cluster_col <- cluster_col
  object@misc$def_NN_object <- neighbor

  SeuratObject::Idents(object) <- object@meta.data[[cluster_col]]
  SeuratObject::DefaultAssay(object) <- if ("SCTregr" %in% names(object@assays)) "SCTregr" else "RNA"
  object
}

#' Embedding matrix to knn
#'
#' Build a cosine HNSW KNN graph from selected embedding dimensions.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param cells Character vector of cell barcodes defining the order and subset to export or align.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param k Number of nearest neighbors to use for KNN/SNN construction.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param threads Number of threads passed to BPCells, HNSW, or matrix-stat routines.
#' @return BPCells KNN result with neighbor index and distance matrices aligned
#'   to `cells`.
#' @keywords internal

embedding_matrix_to_knn <- function(embedding_matrix, cells, dims, k, dim_prefix, threads = 1) {
  missing_cells <- base::setdiff(cells, rownames(embedding_matrix))
  if (length(missing_cells) > 0) {
    stop("Embedding matrix is missing ", length(missing_cells), " requested cell(s).")
  }

  graph_input <- select_embedding_dimensions(
    embedding_matrix = embedding_matrix[cells, , drop = FALSE],
    dims = dims,
    dim_prefix = dim_prefix
  )
  k <- min(as.integer(k), nrow(graph_input) - 1L)
  BPCells::knn_hnsw(graph_input, k = k, metric = "cosine", threads = threads, ef = 500)
}

knn_to_sparse_knn_matrix <- function(knn, cell_names) {
  n_cells <- length(cell_names)
  knn_idx <- knn$idx
  knn_idx <- matrix(as.integer(knn_idx), nrow = nrow(knn_idx), dimnames = dimnames(knn_idx))
  if (!identical(rownames(knn_idx), cell_names)) {
    rownames(knn_idx) <- cell_names
  }

  Matrix::sparseMatrix(
    i = rep(seq_len(n_cells), each = ncol(knn_idx)),
    j = as.vector(t(knn_idx)),
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(cell_names, cell_names)
  ) |>
    Matrix::drop0()
}

#' Add graph from sparse matrix
#'
#' Store a sparse adjacency matrix as a Seurat graph for one assay.
#'
#' @param object Seurat object being modified in place and returned.
#' @param name Name of the Seurat graph slot to create.
#' @param graph_matrix Sparse cell-by-cell graph matrix whose row and column names define the Seurat graph cells.
#' @param assay Seurat assay name associated with the reduction, graph, or metadata being added.
#' @return The Seurat object with `graph_matrix` stored under `name`.
#' @keywords internal

add_graph_from_sparse_matrix <- function(object, name, graph_matrix, assay) {
  graph <- SeuratObject::as.Graph(graph_matrix)
  methods::slot(graph, "assay.used") <- assay
  object[[name]] <- graph
  object
}

#' Add neighbor from KNN
#'
#' Store a BPCells KNN result as a Seurat Neighbor object.
#'
#' @param object Seurat object being modified in place and returned.
#' @param name Name of the Seurat neighbor slot to create.
#' @param knn Nearest-neighbor result with index and distance components, usually from BPCells HNSW helpers.
#' @param cells Character vector of cell barcodes defining the order and subset to export or align.
#' @return The Seurat object with `knn$idx` and `knn$dist` stored under
#'   `object@neighbors[[name]]`.
#' @keywords internal

add_neighbor_from_knn <- function(object, name, knn, cells) {
  object@neighbors[[name]] <- methods::new(
    "Neighbor",
    nn.idx = as.matrix(knn$idx),
    nn.dist = as.matrix(knn$dist),
    alg.idx = NULL,
    alg.info = list(),
    cell.names = cells
  )
  object
}

align_WNN_knn_to_cells <- function(WNN_results, cells) {
  wnn_cells <- rownames(WNN_results$nn_idx)
  if (!base::setequal(cells, wnn_cells)) {
    stop("Final Seurat cells must exactly match WNN_results cells.")
  }

  old_to_new <- match(seq_along(wnn_cells), match(cells, wnn_cells))
  nn_idx <- WNN_results$nn_idx[cells, , drop = FALSE]
  nn_idx[] <- old_to_new[nn_idx]
  nn_dist <- WNN_results$nn_dist[cells, , drop = FALSE]
  rownames(nn_idx) <- rownames(nn_dist) <- cells
  list(idx = nn_idx, dist = nn_dist)
}

#' Build seurat signac convenience object
#'
#' Build the legacy Seurat/Signac export object from BPCells-native matrices and metadata.
#'
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param PCA_results List returned by the GEX PCA helper, including cell embeddings, feature loadings, and variable features.
#' @param GEX_harmony_embeddings Harmony-corrected GEX embedding matrix aligned by barcode.
#' @param GEX_UMAP_embeddings_tibble Tibble of GEX UMAP coordinates keyed by `barcode_w_prefix`.
#' @param GEX_non_harmony_UMAP_embeddings_tibble Tibble of non-Harmony GEX UMAP coordinates keyed by `barcode_w_prefix`.
#' @param gene_features_df Gene metadata data frame with Ensembl IDs, symbols, coordinates, or feature names used for ID conversion and annotation.
#' @param GEX_dims GEX dimensions used when reconstructing reductions and neighbors in the export object.
#' @param data_nNNs Nearest-neighbor count used for graph/neighbor reconstruction in the export object.
#' @param ATAC_peak_matrix Optional peak-by-cell ATAC matrix used for the Signac assay.
#' @param ATAC_peak_GRanges Optional GRanges for ATAC peaks, aligned by peak name to the corresponding ATAC matrix rows.
#' @param ATAC_annotated_peak_GRanges Optional ATAC peak GRanges with annotation metadata used for Signac feature metadata or locus plots.
#' @param ATAC_LSI_results Optional list returned by ATAC LSI with cell embeddings, feature loadings, and singular values.
#' @param ATAC_harmony_embeddings Optional Harmony-corrected ATAC embedding matrix aligned by barcode.
#' @param WNN_results Optional list returned by WNN helpers, including KNN/SNN structures and modality weights aligned by barcode.
#' @param ATAC_UMAP_embeddings_tibble Optional tibble of ATAC UMAP coordinates keyed by `barcode_w_prefix`.
#' @param WNN_UMAP_embeddings_tibble Optional tibble of WNN UMAP coordinates keyed by `barcode_w_prefix`.
#' @param TF_activity_matrix Optional TF-by-cell activity matrix, usually chromVAR deviation scores or a derived activity assay.
#' @param signac_annotation_GRanges Optional GRanges object containing signac annotation GRanges coordinates and metadata.
#' @param fragment_records_tibble Optional tibble describing fragment files and cell sets used to create Signac Fragment objects.
#' @param ATAC_dims Optional ATAC dimensions used when reconstructing reductions and neighbors in the export object.
#' @param graph_threads Thread count used for KNN graph reconstruction during export.
#' @return A Seurat/Signac object or component with metadata, assays, reductions, and graphs aligned by cell barcode.
#' @keywords internal

build_seurat_signac_convenience_object <- function(
  GEX_counts_matrix,
  metadata_tibble,
  PCA_results,
  GEX_harmony_embeddings,
  GEX_UMAP_embeddings_tibble,
  GEX_non_harmony_UMAP_embeddings_tibble,
  gene_features_df,
  GEX_dims,
  data_nNNs,
  ATAC_peak_matrix = NULL,
  ATAC_peak_GRanges = NULL,
  ATAC_annotated_peak_GRanges = NULL,
  ATAC_LSI_results = NULL,
  ATAC_harmony_embeddings = NULL,
  WNN_results = NULL,
  ATAC_UMAP_embeddings_tibble = NULL,
  WNN_UMAP_embeddings_tibble = NULL,
  TF_activity_matrix = NULL,
  signac_annotation_GRanges = NULL,
  fragment_records_tibble = NULL,
  ATAC_dims = NULL,
  graph_threads = 1
) {
  cells <- get_seurat_export_cells(metadata_tibble)
  check_seurat_export_cells(
    cells = cells,
    named_cell_sets = purrr::compact(list(
      GEX_counts_matrix = colnames(GEX_counts_matrix),
      PCA_results = rownames(PCA_results$cell_embeddings),
      GEX_harmony_embeddings = rownames(GEX_harmony_embeddings),
      ATAC_peak_matrix = colnames(ATAC_peak_matrix),
      ATAC_LSI_results = rownames(ATAC_LSI_results$cell_embeddings),
      ATAC_harmony_embeddings = rownames(ATAC_harmony_embeddings),
      WNN_results = rownames(WNN_results$nn_idx),
      TF_activity_matrix = colnames(TF_activity_matrix)
    ))
  )

  metadata_df <- prepare_seurat_export_metadata(metadata_tibble, cells, WNN_results)
  RNA_counts <- align_matrix_columns(GEX_counts_matrix, cells, "GEX_counts_matrix")

  object <- SeuratObject::CreateSeuratObject(
    counts = RNA_counts,
    assay = "RNA",
    meta.data = metadata_df,
    project = "SeuratProject"
  )
  object[["RNA"]] <- add_feature_metadata(
    assay = object[["RNA"]],
    feature_metadata_df = get_RNA_feature_metadata(gene_features_df, rownames(object[["RNA"]]))
  )

  if (!is.null(ATAC_peak_matrix)) {
    peak_counts <- align_matrix_columns(ATAC_peak_matrix, cells, "ATAC_peak_matrix")
    peak_ranges <- align_peak_GRanges_to_matrix(ATAC_peak_GRanges, rownames(peak_counts))
    object[["ATAC"]] <- Signac::CreateGRangesAssay(counts = peak_counts, ranges = peak_ranges)
    Signac::Annotation(object[["ATAC"]]) <- ensure_signac_annotation_GRanges(signac_annotation_GRanges)
    Signac::Fragments(object[["ATAC"]]) <- create_signac_fragment_objects(fragment_records_tibble)
    SeuratObject::VariableFeatures(object[["ATAC"]]) <- rownames(object[["ATAC"]])
    object[["ATAC"]] <- add_feature_metadata(
      assay = object[["ATAC"]],
      feature_metadata_df = get_peak_feature_metadata(ATAC_annotated_peak_GRanges, rownames(peak_counts))
    )
  }

  sct_features <- rownames(RNA_counts)[BPCells::rowSums(RNA_counts) > 50]
  if (length(sct_features) == 0) {
    sct_features <- rownames(RNA_counts)
  }
  SCTregr_counts <- align_matrix_rows(RNA_counts, sct_features, "RNA_counts")
  object[["SCTregr"]] <- SeuratObject::CreateAssay5Object(counts = SCTregr_counts)
  SeuratObject::VariableFeatures(object[["SCTregr"]]) <- base::intersect(PCA_results$variable_features, rownames(object[["SCTregr"]]))

  if (!is.null(TF_activity_matrix)) {
    TF_activity_data <- align_matrix_columns(TF_activity_matrix, cells, "TF_activity_matrix")
    object[["TF_activity"]] <- SeuratObject::CreateAssay5Object(data = TF_activity_data)
  }

  object <- add_dimreduc_from_matrix(
    object = object,
    name = "PCA",
    embeddings = PCA_results$cell_embeddings,
    loadings = PCA_results$feature_loadings,
    assay = "SCTregr",
    key = "PC_",
    cells = cells
  )
  object <- add_dimreduc_from_matrix(
    object = object,
    name = "PCA_UMAP",
    embeddings = embedding_tibble_to_matrix(GEX_non_harmony_UMAP_embeddings_tibble),
    assay = "SCTregr",
    key = "PCAUMAP_",
    cells = cells
  )
  object <- add_dimreduc_from_matrix(
    object = object,
    name = "PCA_harmony",
    embeddings = GEX_harmony_embeddings,
    assay = "SCTregr",
    key = "PCAharmony_",
    cells = cells
  )
  object <- add_dimreduc_from_matrix(
    object = object,
    name = "PCA_harmony_UMAP",
    embeddings = embedding_tibble_to_matrix(GEX_UMAP_embeddings_tibble),
    assay = "SCTregr",
    key = "PCAharmonyUMAP_",
    cells = cells
  )
  if ("ATAC" %in% names(object@assays)) {
    object <- add_dimreduc_from_matrix(
      object = object,
      name = "LSI",
      embeddings = ATAC_LSI_results$cell_embeddings,
      loadings = ATAC_LSI_results$feature_loadings,
      assay = "ATAC",
      key = "LSI_",
      cells = cells
    )
    object <- add_dimreduc_from_matrix(
      object = object,
      name = "LSI_harmony",
      embeddings = ATAC_harmony_embeddings,
      assay = "ATAC",
      key = "LSIharmony_",
      cells = cells
    )
    object <- add_dimreduc_from_matrix(
      object = object,
      name = "LSI_harmony_UMAP",
      embeddings = embedding_tibble_to_matrix(ATAC_UMAP_embeddings_tibble),
      assay = "ATAC",
      key = "LSIharmonyUMAP_",
      cells = cells
    )
  }
  object <- add_dimreduc_from_matrix(
    object = object,
    name = "WNN_harmony_NN_UMAP",
    embeddings = embedding_tibble_to_matrix(WNN_UMAP_embeddings_tibble),
    assay = if ("ATAC" %in% names(object@assays)) "ATAC" else "SCTregr",
    key = "WNNharmonyNNUMAP_",
    cells = cells
  )

  GEX_knn <- embedding_matrix_to_knn(GEX_harmony_embeddings, cells, GEX_dims, data_nNNs, "PCA_", graph_threads)
  object <- add_graph_from_sparse_matrix(object, "PCA_harmony_NN", knn_to_sparse_knn_matrix(GEX_knn, cells), "SCTregr")
  object <- add_graph_from_sparse_matrix(object, "PCA_harmony_SNN", get_SNN_matrix_from_knn(GEX_knn, cells), "SCTregr")

  if ("ATAC" %in% names(object@assays) && !is.null(ATAC_harmony_embeddings)) {
    ATAC_knn <- embedding_matrix_to_knn(ATAC_harmony_embeddings, cells, ATAC_dims, data_nNNs, "LSI_", graph_threads)
    object <- add_graph_from_sparse_matrix(object, "LSI_harmony_NN", knn_to_sparse_knn_matrix(ATAC_knn, cells), "ATAC")
    object <- add_graph_from_sparse_matrix(object, "LSI_harmony_SNN", get_SNN_matrix_from_knn(ATAC_knn, cells), "ATAC")
  }

  if (!is.null(WNN_results)) {
    WNN_knn <- align_WNN_knn_to_cells(WNN_results, cells)
    wnn_assay <- if ("ATAC" %in% names(object@assays)) "ATAC" else "SCTregr"
    object <- add_graph_from_sparse_matrix(object, "WNN_harmony_KNN", knn_to_sparse_knn_matrix(WNN_knn, cells), wnn_assay)
    object <- add_graph_from_sparse_matrix(object, "WNN_harmony_SNN", get_SNN_matrix_from_knn(WNN_knn, cells), wnn_assay)
    object <- add_neighbor_from_knn(object, "WNN_harmony_NN", WNN_knn, cells)
  }

  set_seurat_export_defaults(object)
}

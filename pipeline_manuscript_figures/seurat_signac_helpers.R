comparison_aggregation_names <- function() {
  paste0("comparison_", c(1, 2, 5, 10, 20), "x")
}

# Evaluates expr with one forked future worker per allocated core, then restores
# the previous plan: crew workers run several targets, and a plan left in place
# would oversubscribe the threads of Seurat's own kernels in later targets.
comparison_with_workers <- function(expr) {
  withr::local_options(future.globals.maxSize = Inf)
  with(future::plan(future::multicore, workers = parallelly::availableCores()), local = TRUE)
  expr
}

comparison_build_config <- function(aggregation) {
  if (!aggregation %in% comparison_aggregation_names()) {
    stop("Unsupported comparison aggregation: ", aggregation, call. = FALSE)
  }

  aggregation_config <- benchmark_aggregation_tibble(aggregation)
  sample_tibble <- tibble::tibble(
    sample_id = aggregation_config$aggregation_GEM_well_IDs[[1]],
    cellranger_count_dir = aggregation_config$benchmark_cellranger_count_dirs[[1]]
  ) |>
    dplyr::mutate(
      matrix_h5 = file.path(.data$cellranger_count_dir, "outs", "filtered_feature_bc_matrix.h5"),
      fragment_file = file.path(.data$cellranger_count_dir, "outs", "atac_fragments.tsv.gz"),
      fragment_index_file = paste0(.data$fragment_file, ".tbi")
    )

  list(
    aggregation = aggregation,
    sample_tibble = sample_tibble,
    RNA_dims = aggregation_config$aggregation_GEX_data_PCs[[1]],
    ATAC_dims = aggregation_config$aggregation_ATAC_data_PCs[[1]],
    GEX_cluster_resolution = aggregation_config$aggregation_GEX_cluster_res[[1]],
    WNN_cluster_resolution = aggregation_config$aggregation_WNN_cluster_res[[1]],
    WNN_k = aggregation_config$aggregation_data_nNNs[[1]],
    UMAP_n_neighbors = aggregation_config$aggregation_UMAP_nNNs[[1]],
    UMAP_min_dist = aggregation_config$aggregation_UMAP_min_dist[[1]],
    genome = "GRCh38",
    random_seed = 1234L
  )
}

comparison_input_files <- function(config) {
  unlist(
    config$sample_tibble[c("matrix_h5", "fragment_file", "fragment_index_file")],
    use.names = FALSE
  )
}

# One BPCells directory of GEM-well-prefixed GEX counts per GEM well.
comparison_write_RNA_dirs <- function(config, out_dir) {
  purrr::map2_chr(config$sample_tibble$sample_id, config$sample_tibble$matrix_h5, \(sample_id, matrix_h5) {
    counts <- BPCells::open_matrix_10x_hdf5(matrix_h5, feature_type = "Gene Expression")
    colnames(counts) <- paste0(sample_id, "_", colnames(counts))
    dir <- file.path(out_dir, sample_id)
    BPCells::write_matrix_dir(counts, dir, overwrite = TRUE)
    dir
  })
}

comparison_process_RNA <- function(RNA_dirs, config) {
  # SCTransform passes its gene models to future even in a sequential plan.
  withr::local_options(future.globals.maxSize = Inf)
  counts <- purrr::map(RNA_dirs, open_BPCells_dir)
  object <- SeuratObject::CreateSeuratObject(purrr::reduce(counts, cbind), assay = "RNA")
  object$GEM_well_ID <- rep(config$sample_tibble$sample_id, purrr::map_int(counts, ncol))
  object <- Seurat::SCTransform(object, conserve.memory = TRUE, seed.use = config$random_seed)
  object <- Seurat::RunPCA(object, npcs = max(config$RNA_dims), seed.use = config$random_seed)
  object <- Seurat::FindNeighbors(object, reduction = "pca", dims = config$RNA_dims, k.param = config$WNN_k)
  object <- Seurat::FindClusters(object, resolution = config$GEX_cluster_resolution, algorithm = 4L,
    random.seed = config$random_seed)
  object$RNA_clusters <- object$seurat_clusters
  Seurat::RunUMAP(object, reduction = "pca", dims = config$RNA_dims, n.neighbors = config$UMAP_n_neighbors,
    min.dist = config$UMAP_min_dist, reduction.name = "umap.rna", reduction.key = "rnaUMAP_",
    seed.use = config$random_seed)
}

comparison_call_consensus_peaks <- function(
  fragment_files,
  cluster_names,
  config,
  blacklist_GRanges
) {
  fragment_files <- unlist(fragment_files, use.names = FALSE)
  cluster_names <- as.character(cluster_names)
  if (length(fragment_files) != length(cluster_names)) {
    stop("Peak-calling fragment files and cluster names must have equal lengths.", call. = FALSE)
  }

  within_cluster_peaks <- comparison_with_workers(future.apply::future_lapply(
    seq_along(fragment_files),
    \(cluster_index) {
      cluster_name <- cluster_names[[cluster_index]]
      narrow_peak_file <- call_peaks_w_MACS3(
        ATAC_fragments_per_cluster = fragment_files[[cluster_index]],
        ATAC_peak_calling_cluster_names = cluster_name,
        genome = config$genome,
        output_suffix = cluster_name,
        allow_no_peaks = TRUE
      )
      on.exit(unlink(narrow_peak_file), add = TRUE)

      get_peak_GRanges_w_fixed_width(
        narrow_peak_file,
        genome = config$genome,
        blacklist_GRanges = blacklist_GRanges
      ) |>
        list() |>
        combine_collapse_GRanges_ArchR(by = "neg_log10pvalue_summit") |>
        format_peak_GRanges(cluster_name = cluster_name)
    },
    future.seed = config$random_seed
  ))

  combine_collapse_GRanges_ArchR(
    within_cluster_peaks,
    by = "fold_change"
  )
}

# Object cell names map to the barcodes of the GEM well's fragment file.
comparison_fragments <- function(fragment_file, sample_id, cells, ...) {
  Signac::CreateFragmentObject(fragment_file, cells = stats::setNames(sub(paste0("^", sample_id, "_"), "", cells), cells), ...)
}

# Quantifies the consensus peaks per GEM well, as Signac's merging vignette does,
# with FeatureMatrix()'s default fragtk backend, into one BPCells directory each.
comparison_write_ATAC_dirs <- function(RNA_dirs, config, peaks, out_dir) {
  unlist(comparison_with_workers(future.apply::future_Map(\(RNA_dir, sample_id, fragment_file) {
    dir <- file.path(out_dir, sample_id)
    unlink(dir, recursive = TRUE)
    fragments <- comparison_fragments(fragment_file, sample_id, colnames(BPCells::open_matrix_dir(RNA_dir)))
    Signac::FeatureMatrix(fragments, features = peaks, bpcells = TRUE, bpcells.dir = dir)
    dir
  }, RNA_dirs, config$sample_tibble$sample_id, config$sample_tibble$fragment_file, future.seed = TRUE)), use.names = FALSE)
}

comparison_build_object <- function(RNA_object, ATAC_dirs, config) {
  cells <- split(colnames(RNA_object), RNA_object$GEM_well_ID)[config$sample_tibble$sample_id]
  fragments <- purrr::pmap(list(config$sample_tibble$fragment_file, config$sample_tibble$sample_id, cells),
    comparison_fragments, validate.fragments = FALSE)
  counts <- purrr::reduce(purrr::map(ATAC_dirs, open_BPCells_dir), cbind)
  object <- RNA_object
  object[["ATAC"]] <- Signac::CreateChromatinAssay5(counts = counts[, colnames(object)], fragments = fragments,
    validate.fragments = FALSE)
  SeuratObject::DefaultAssay(object) <- "ATAC"
  object <- Signac::RunTFIDF(object)
  object <- Signac::FindTopFeatures(object, min.cutoff = "q0")
  object <- Signac::RunSVD(object, n = max(config$ATAC_dims))
  object <- Seurat::FindMultiModalNeighbors(object, reduction.list = list("pca", "lsi"),
    dims.list = list(config$RNA_dims, config$ATAC_dims), k.nn = config$WNN_k)
  object <- Seurat::RunUMAP(object, nn.name = "weighted.nn", min.dist = config$UMAP_min_dist,
    reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = config$random_seed)
  Seurat::FindClusters(object, graph.name = "wsnn", algorithm = 3L, resolution = config$WNN_cluster_resolution,
    random.seed = config$random_seed)
}

comparison_summarize <- function(object, config, peaks) {
  tibble::tibble(
    aggregation = config$aggregation,
    GEM_well_count = nrow(config$sample_tibble),
    cells = ncol(object),
    peaks = length(peaks),
    RNA_features = nrow(object[["RNA"]]),
    WNN_clusters = dplyr::n_distinct(object$seurat_clusters),
    assays = paste(names(object@assays), collapse = ","),
    reductions = paste(names(object@reductions), collapse = ","),
    graphs = paste(names(object@graphs), collapse = ","),
    object_size_bytes = as.numeric(object.size(object))
  )
}

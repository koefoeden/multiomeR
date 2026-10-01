comparison_aggregation_names <- function() {
  paste0("comparison_", c(1, 2, 5, 10, 20), "x")
}

comparison_configure_parallelism <- function() {
  workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
  if (is.na(workers) || workers < 1L) {
    workers <- 1L
  }
  if (workers > 1L) {
    options(future.globals.maxSize = 160 * 1024^3)
    future::plan(future::multicore, workers = workers)
  }
  invisible(workers)
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

comparison_create_sample_object <- function(counts_list, sample_id, fragment_file, cells) {
  required_assays <- c("Gene Expression", "Peaks")
  if (!all(required_assays %in% names(counts_list))) {
    stop("Input does not contain Gene Expression and Peaks matrices.", call. = FALSE)
  }

  object <- SeuratObject::CreateSeuratObject(
    counts = counts_list[["Gene Expression"]][, cells, drop = FALSE],
    assay = "RNA",
    project = sample_id,
    min.cells = 0L,
    min.features = 0L
  )
  object[["ATAC"]] <- Signac::CreateChromatinAssay5(
    counts = counts_list[["Peaks"]][, cells, drop = FALSE],
    fragments = fragment_file,
    validate.fragments = TRUE,
    verbose = FALSE
  )
  object$GEM_well_ID <- sample_id
  SeuratObject::RenameCells(object, add.cell.id = sample_id)
}

comparison_read_inputs <- function(config) {
  sample_objects <- purrr::map2(
    config$sample_tibble$sample_id,
    config$sample_tibble$matrix_h5,
    \(sample_id, matrix_h5) {
      counts_list <- Seurat::Read10X_h5(matrix_h5, use.names = TRUE)
      cells <- colnames(counts_list[["Gene Expression"]])
      fragment_file <- config$sample_tibble$fragment_file[
        match(sample_id, config$sample_tibble$sample_id)
      ]
      comparison_create_sample_object(counts_list, sample_id, fragment_file, cells)
    }
  )

  object <- purrr::reduce(
    sample_objects,
    \(x, y) merge(x, y, merge.data = FALSE, merge.dr = FALSE)
  )
  object[["RNA"]] <- SeuratObject::JoinLayers(object[["RNA"]])
  object[["ATAC"]] <- SeuratObject::JoinLayers(object[["ATAC"]])
  object
}

comparison_preprocess_RNA <- function(object, config) {
  comparison_configure_parallelism()
  SeuratObject::DefaultAssay(object) <- "RNA"
  object <- Seurat::SCTransform(
    object,
    conserve.memory = TRUE,
    seed.use = config$random_seed,
    verbose = TRUE
  )
  Seurat::RunPCA(
    object,
    npcs = max(config$RNA_dims),
    seed.use = config$random_seed,
    verbose = TRUE
  )
}

comparison_cluster_RNA <- function(object, config) {
  comparison_configure_parallelism()
  object <- Seurat::FindNeighbors(
    object,
    reduction = "pca",
    dims = config$RNA_dims,
    k.param = config$WNN_k,
    verbose = TRUE
  )
  object <- Seurat::FindClusters(
    object,
    resolution = config$GEX_cluster_resolution,
    algorithm = 4L,
    random.seed = config$random_seed,
    verbose = TRUE
  )
  object$RNA_clusters <- object$seurat_clusters
  Seurat::RunUMAP(
    object,
    reduction = "pca",
    dims = config$RNA_dims,
    n.neighbors = config$UMAP_n_neighbors,
    min.dist = config$UMAP_min_dist,
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_",
    seed.use = config$random_seed,
    verbose = TRUE
  )
}

comparison_call_consensus_peaks <- function(
  fragment_files,
  cluster_names,
  config,
  blacklist_GRanges
) {
  comparison_configure_parallelism()
  fragment_files <- unlist(fragment_files, use.names = FALSE)
  cluster_names <- as.character(cluster_names)
  if (length(fragment_files) != length(cluster_names)) {
    stop("Peak-calling fragment files and cluster names must have equal lengths.", call. = FALSE)
  }

  within_cluster_peaks <- future.apply::future_lapply(
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
  )

  combine_collapse_GRanges_ArchR(
    within_cluster_peaks,
    by = "fold_change"
  )
}

comparison_rebuild_ATAC_assay <- function(object, peaks) {
  comparison_configure_parallelism()
  old_ATAC_assay <- object[["ATAC"]]
  counts <- Signac::FeatureMatrix(
    object = old_ATAC_assay,
    features = peaks,
    cells = colnames(object),
    fragtk = FALSE,
    process_n = 2000L,
    verbose = TRUE
  )
  object[["ATAC"]] <- Signac::CreateChromatinAssay5(
    counts = counts,
    fragments = Signac::Fragments(old_ATAC_assay),
    validate.fragments = FALSE,
    verbose = FALSE
  )
  object
}

comparison_preprocess_ATAC <- function(object, config) {
  comparison_configure_parallelism()
  SeuratObject::DefaultAssay(object) <- "ATAC"
  object <- Signac::RunTFIDF(object)
  object <- Signac::FindTopFeatures(object, min.cutoff = "q0")
  object <- Signac::RunSVD(
    object,
    n = max(config$ATAC_dims)
  )
  if (is.null(object[["lsi"]])) {
    stop("Signac::RunSVD() did not create the expected lsi reduction.", call. = FALSE)
  }
  object
}

comparison_run_WNN <- function(object, config) {
  comparison_configure_parallelism()
  object <- Seurat::FindMultiModalNeighbors(
    object,
    reduction.list = list("pca", "lsi"),
    dims.list = list(config$RNA_dims, config$ATAC_dims),
    k.nn = config$WNN_k,
    verbose = TRUE
  )
  object <- Seurat::RunUMAP(
    object,
    nn.name = "weighted.nn",
    min.dist = config$UMAP_min_dist,
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    seed.use = config$random_seed,
    verbose = TRUE
  )
  Seurat::FindClusters(
    object,
    graph.name = "wsnn",
    algorithm = 3L,
    resolution = config$WNN_cluster_resolution,
    random.seed = config$random_seed,
    verbose = TRUE
  )
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

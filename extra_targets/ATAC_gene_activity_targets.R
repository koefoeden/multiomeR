rlang::list2(
  targets::tar_target(
    name = ATAC_cells_per_GEM_well_tibble.ATAC,
    description = "Pair each GEM well's ATAC-retained nuclei with its fragment directory for per-well gene scores",
    command = split_cells_by_GEM_well(
      metadata_w_cell_types_tibble.ATAC,
      GEM_well_IDs = aggregation_GEM_well_IDs,
      fragment_dirs = aggregation_fragments_w_prefix_bpcells_dir_syms
    ),
    iteration = "vector"
  ),
  tarchetypes::tar_file(
    name = gene_score_archr_BPCells_matrix_dirs_per_GEM_well.ATAC,
    description = "Compute ArchR-style ATAC gene activity scores for one GEM well's nuclei with BPCells and write them to disk",
    command = {
      out_dir <- get_structured_file_path(override_suffix = ATAC_cells_per_GEM_well_tibble.ATAC$GEM_well_ID)
      tile_matrix_dir <- tempfile(pattern = "gene_score_tiles_", tmpdir = dirname(out_dir))
      on.exit(if (fs::dir_exists(tile_matrix_dir)) fs::dir_delete(tile_matrix_dir), add = TRUE)

      fragments <- open_GEM_well_fragments(ATAC_cells_per_GEM_well_tibble.ATAC, aggregated_cellranger_ref_list$genomes[[1]]) |>
        BPCells::select_cells(ATAC_cells_per_GEM_well_tibble.ATAC$barcodes[[1]])
      genes <- marker_validated_Ensembl_annotations_GRanges_list$genes
      chromosome_sizes <- get_chrom_sizes_for_BPCells_tile_calling(aggregated_cellranger_ref_list$genomes[[1]])
      genes <- filter_GRanges_to_chrom_sizes(genes, chromosome_sizes)
      blacklist <- filter_GRanges_to_chrom_sizes(blacklist_GRanges.ATAC, chromosome_sizes, allow_empty = TRUE)

      if (!"gene_name" %in% colnames(S4Vectors::mcols(genes))) {
        genes$gene_name <- names(genes)
      }
      # Put genes with no retained ArchR tiles first so BPCells can keep a
      # full zero row for them after blacklist filtering.
      genes <- order_genes_for_BPCells_gene_score_archr(
        genes = genes,
        chromosome_sizes = chromosome_sizes,
        blacklist = blacklist
      )
      gene_scores <- BPCells::gene_score_archr(
        fragments = fragments,
        genes = genes,
        chromosome_sizes = chromosome_sizes,
        blacklist = blacklist,
        gene_name_column = "gene_name",
        tile_matrix_path = tile_matrix_dir
      )

      BPCells::write_matrix_dir(gene_scores, out_dir, overwrite = TRUE)
      out_dir
    },
    pattern = map(ATAC_cells_per_GEM_well_tibble.ATAC),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = gene_score_archr_BPCells_matrix.ATAC,
    description = "Open the per-GEM-well BPCells ArchR-style ATAC gene activity scores as one matrix",
    command = do.call(cbind, lapply(gene_score_archr_BPCells_matrix_dirs_per_GEM_well.ATAC, open_BPCells_dir)),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = marker_gene_activity_dot_plot.7_ATAC_QC,
    description = "Dot plot of marker gene activity scores per ATAC cell type. [checkpoint:7_ATAC-QC]",
    command = {
      plot <- plot_marker_gene_activity_dot_BPCells(
        feature_matrix = gene_score_archr_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        marker_genes_list = UCell_GEX_marker_genes_list,
        group_col = "LSI_harmony_SNN_cluster_cell_type"
      )
      save_plots_structured(plot, width = max(10, 4 + 0.35 * nlevels(plot$data$marker_feature)))
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)

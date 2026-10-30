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
  targets::tar_target(
    name = gene_score_archr_weights.ATAC,
    description = "Build the ArchR-style gene-by-tile weights that every GEM well's gene activity scores share [resource_observation: peak_RAM_GiB=2.8; genes=32k; tiles=6.18M; aggregation=unpublished; method=VmHWM; date=2026-10-21]",
    command = get_gene_score_archr_weights(
      genes = marker_validated_Ensembl_annotations_GRanges_list$genes,
      blacklist = blacklist_GRanges.ATAC,
      genome = aggregated_cellranger_ref_list$genomes[[1]]
    ),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = gene_score_archr_BPCells_matrix_dirs_per_GEM_well.ATAC,
    description = "Compute ArchR-style ATAC gene activity scores for one GEM well's nuclei with BPCells and write them to disk [resource_observation: peak_RAM_GiB=0.7; n_nuclei=15k per well; aggregation=unpublished; method=VmHWM; date=2026-10-21]",
    command = {
      out_dir <- get_structured_file_path(override_suffix = ATAC_cells_per_GEM_well_tibble.ATAC$GEM_well_ID)
      unlink(out_dir, recursive = TRUE) # BPCells writes into an existing folder twice
      open_GEM_well_fragments(ATAC_cells_per_GEM_well_tibble.ATAC, aggregated_cellranger_ref_list$genomes[[1]]) |>
        BPCells::select_cells(ATAC_cells_per_GEM_well_tibble.ATAC$barcodes[[1]]) |>
        write_gene_score_archr(gene_score_archr_weights.ATAC, out_dir)
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

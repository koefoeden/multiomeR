rlang::list2(
  tarchetypes::tar_file(
    name = peaks_per_cluster_narrowPeaks.peaks.ATAC,
    description = "Call full-genome ATAC peaks per cluster using capped discovery cells [part_of_graph:ATAC] [part_of_graph:seurat_export]",
    command = call_peaks_w_BPCells_tile(
      ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
      ATAC_BCs_per_peak_cluster = peak_calling_cluster_discovery_tibble.ATAC$BCs_for_peak_discovery[[1]],
      ATAC_peak_calling_cluster_names = peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]],
      genome = aggregated_cellranger_ref_list$genomes[[1]],
      output_suffix = peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]],
      allow_no_peaks = TRUE
    ),
    pattern = map(peak_calling_cluster_discovery_tibble.ATAC),
    resources = get_tar_resources(RAM_GB_req = 60)
  )
)

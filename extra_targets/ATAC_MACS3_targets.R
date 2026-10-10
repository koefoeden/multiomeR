rlang::list2(
  tarchetypes::tar_file(
    name = fragments_per_peak_calling_cluster_discovery.fragments.ATAC,
    description = "Export the full-genome ATAC fragments of each half of a cluster's peak-discovery nuclei",
    command = write_ATAC_fragments_for_peak_calling_cluster(
      ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
      BCs_per_peak_cluster = peak_calling_cluster_discovery_tibble.ATAC$BCs_for_peak_discovery[[1]],
      peak_calling_cluster_name = peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]],
      output_suffix = paste0(peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]], "__discovery")
    ),
    pattern = map(peak_calling_cluster_discovery_tibble.ATAC)
  ),
  tarchetypes::tar_file(
    name = peaks_per_cluster_narrowPeaks.peaks.ATAC,
    description = "Call full-genome ATAC peaks in each half of a cluster's capped discovery cells [part_of_graph:ATAC] [part_of_graph:seurat_export]",
    command = call_peaks_w_MACS3(
      ATAC_fragments_per_cluster = fragments_per_peak_calling_cluster_discovery.fragments.ATAC,
      ATAC_peak_calling_cluster_names = peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]],
      genome = aggregated_cellranger_ref_list$genomes[[1]],
      output_suffix = peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]],
      allow_no_peaks = TRUE
    ),
    pattern = map(fragments_per_peak_calling_cluster_discovery.fragments.ATAC, peak_calling_cluster_discovery_tibble.ATAC),
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)

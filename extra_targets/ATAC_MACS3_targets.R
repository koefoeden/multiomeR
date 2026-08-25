rlang::list2(
  tarchetypes::tar_file(
    name = fragments_per_peak_calling_cluster_discovery.fragments.ATAC,
    description = "Export full-genome ATAC fragments used for per-cluster peak discovery",
    command = write_ATAC_fragments_for_peak_calling_cluster(
      ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
      BCs_per_peak_cluster = peak_calling_cluster_discovery_tibble.ATAC$BCs_for_peak_discovery[[1]],
      peak_calling_cluster_name = peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]],
      output_suffix = paste0(peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]], "__discovery")
    ),
    pattern = map(peak_calling_cluster_discovery_tibble.ATAC)
  )
)

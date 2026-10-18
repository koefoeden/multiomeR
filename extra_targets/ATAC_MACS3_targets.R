rlang::list2(
  tarchetypes::tar_file(
    name = peak_calling_cut_sites_per_GEM_well.fragments.ATAC,
    description = "Write one GEM well's ATAC cut sites of each half of the clusters' peak-discovery nuclei in one pass",
    command = write_peak_calling_cut_sites_for_GEM_well(
      GEM_well_cells = GEX_cells_per_GEM_well_tibble.ATAC,
      discovery_tibble = peak_calling_cluster_discovery_tibble.ATAC,
      genome = aggregated_cellranger_ref_list$genomes[[1]]
    ),
    pattern = map(GEX_cells_per_GEM_well_tibble.ATAC),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peaks_per_cluster_narrowPeaks.peaks.ATAC,
    description = "Call full-genome ATAC peaks in each half of a cluster's capped discovery cells [part_of_graph:ATAC] [part_of_graph:seurat_export] [resource_observation: peak_RAM_GiB=6.4; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10] [resource_observation: peak_RAM_GiB=13.6; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-01]",
    command = call_peaks_w_MACS3(
      ATAC_cut_site_files = get_cut_site_files(
        peak_calling_cut_sites_per_GEM_well.fragments.ATAC,
        peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]]
      ),
      ATAC_peak_calling_cluster_names = peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]],
      genome = aggregated_cellranger_ref_list$genomes[[1]],
      output_suffix = peak_calling_cluster_discovery_tibble.ATAC$replicate_name[[1]],
      allow_no_peaks = TRUE
    ),
    pattern = map(peak_calling_cluster_discovery_tibble.ATAC),
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)

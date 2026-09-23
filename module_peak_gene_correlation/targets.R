peak_gene_correlation_tibble <- build_module_tibble("peak_gene_correlation", aggregation_tibble, aggregation_tibble_all_from_yaml) |>
  dplyr::mutate(peak_gene_correlation_target_suffix =
    paste("peak_gene_correlation", aggregation, sep = ".")) |>
  add_aggregation_target_syms(c(
    "marker_validated_Ensembl_annotations_GRanges_list",
    "aggregated_counts_BPCells_matrix.GEX",
    "consensus_peak_GRanges.ATAC",
    "consensus_peak_BPCells_matrix.ATAC",
    "metadata_w_cell_types_tibble.WNN",
    "harmony_embeddings_matrix.ATAC",
    "combined_BPCells_fragment_obj.ATAC"
  ))

if (nrow(peak_gene_correlation_tibble) == 0L) list() else tarchetypes::tar_map(
  values = peak_gene_correlation_tibble,
  names = peak_gene_correlation_target_suffix,
  descriptions = NULL,
  delimiter = ".",
  source("module_peak_gene_correlation/correlation_targets.R")$value
)

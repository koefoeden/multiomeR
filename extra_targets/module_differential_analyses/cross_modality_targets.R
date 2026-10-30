rlang::list2(
  targets::tar_target(
    name = significant_elements_modality_distribution_tibble,
    description = "Compare each contrast's share of significant elements across pseudobulk differential modalities",
    command = dplyr::bind_rows(
      gene_expression = significant_elements_tibble.gene_expression,
      chromatin_accessibility = significant_elements_tibble.chromatin_accessibility,
      motif_family_accessibility = significant_elements_tibble.motif_family_accessibility,
      TF_motif_accessibility = significant_elements_tibble.TF_motif_accessibility,
      .id = "modality"
    ) |>
      dplyr::summarise(n_significant = sum(n_significant), .by = c(modality, model, contrast)) |>
      dplyr::mutate(
        total_significant = sum(n_significant),
        prop_significant = dplyr::if_else(total_significant > 0, n_significant / total_significant, 0),
        .by = c(modality, model)
      )
  ),
  tarchetypes::tar_file(
    name = significant_elements_modality_distribution_plots,
    description = "Plot each contrast's share of significant elements across pseudobulk differential modalities per model. [checkpoint:differential_analyses]",
    command = significant_elements_modality_distribution_tibble |>
      plot_pseudobulk_differential_significant_elements_modality_distribution() |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = CollecTRI_JASPAR_family_map,
    description = "Map CollecTRI TF and retained TF-complex regulons to JASPAR2026 motif-similarity families [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_family_map(
      CollecTRI_network_tibble = CollecTRI_human_network_tibble,
      JASPAR_motif_family_members_tibble = JASPAR_motif_family_members_tibble
    )
  ),
  targets::tar_target(
    name = CollecTRI_TF_motif_matrix,
    description = "Match the own JASPAR2026 motifs of CollecTRI-mapped regulators to the chromVAR peaks [part_of_graph:differential_analyses] [resource_observation: peak_RAM_GiB=4.7; n_motifs=745; n_peaks=683k; aggregation=unpublished; method=VmHWM; date=2026-10-29]",
    command = {
      peak_ranges <- SummarizedExperiment::rowRanges(chromVAR_obj.ATAC)
      motif_ids <- unlist(strsplit(CollecTRI_JASPAR_family_map$motif_ids, ";", fixed = TRUE))
      get_motif_matrix_from_ATAC_peak_names(
        ATAC_peak_names = get_peak_names_from_GRanges(peak_ranges),
        ATAC_peak_GRanges = peak_ranges,
        motif_matrix_list = JASPAR_motif_matrix_list[
          JASPAR_motif_family_members_tibble$motif_feature[JASPAR_motif_family_members_tibble$motif_id %in% motif_ids]
        ],
        genome = aggregated_cellranger_ref_list$genomes[[1]]
      )
    },
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = TF_motif_root_correlation_tibble,
    description = "Correlate CollecTRI regulators' own-motif accessibility with their family root's across cell types, choosing the own motif below 0.5 [part_of_graph:differential_analyses] [resource_observation: peak_RAM_GiB=17.9; n_pseudobulks=1893; aggregation=unpublished; method=VmHWM; date=2026-10-29]",
    command = get_TF_motif_root_correlation_tibble(
      pseudobulk_ATAC_counts_matrix = pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      TF_motif_matrix = CollecTRI_TF_motif_matrix,
      motif_family_matrix = chromVAR_TF_motif_family_matrix.ATAC,
      motif_family_members_tibble = JASPAR_motif_family_members_tibble
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = pseudobulk_TF_motif_accessibility_matrix,
    description = "Compute pseudobulk accessibility scores of the own motifs that replace their family root in the CollecTRI crosswalk [part_of_graph:differential_analyses] [resource_observation: peak_RAM_GiB=18.3; n_pseudobulks=1893; aggregation=unpublished; method=VmHWM; date=2026-10-29]",
    command = get_pseudobulk_motif_family_accessibility_matrix(
      pseudobulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      chromVAR_motif_family_matrix = CollecTRI_TF_motif_matrix[,
        TF_motif_root_correlation_tibble$motif_feature[TF_motif_root_correlation_tibble$use_own_motif],
        drop = FALSE
      ]
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = regulon_accessibility_source_comparison_tibble.CollecTRI_JASPAR,
    description = "Compare regulon target enrichment, JASPAR motif-family or own-motif accessibility and TF expression results [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_comparison_tibble(
      regulon_results_tibble = results.CollecTRI.regulon_target_enrichment.gene_expression,
      motif_family_accessibility_results_tibble = results_tibble.motif_family_accessibility,
      TF_motif_accessibility_results_tibble = results_tibble.TF_motif_accessibility,
      gene_expression_results_tibble = results_tibble.gene_expression,
      CollecTRI_JASPAR_family_map = CollecTRI_JASPAR_family_map,
      TF_motif_root_correlation_tibble = TF_motif_root_correlation_tibble
    )
  ),
  targets::tar_target(
    name = regulon_accessibility_family_comparison_tibble.CollecTRI_JASPAR,
    description = "Collapse regulon-motif comparisons to motif families or own motifs per model and contrast [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_family_comparison_tibble(
      regulon_accessibility_source_comparison_tibble.CollecTRI_JASPAR
    )
  ),
  targets::tar_target(
    name = regulon_accessibility_concordance_tibble.CollecTRI_JASPAR,
    description = "Summarize regulon-motif correlation, direction concordance, joint significance and a null range per contrast [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_concordance_tibble(
      regulon_accessibility_family_comparison_tibble.CollecTRI_JASPAR
    )
  ),
  tarchetypes::tar_file(
    name = regulon_accessibility_concordance_plots.CollecTRI_JASPAR,
    description = "Plot the concordance of regulon target enrichment and motif-family accessibility per contrast. [checkpoint:differential_analyses]",
    command = regulon_accessibility_concordance_tibble.CollecTRI_JASPAR |>
      plot_CollecTRI_JASPAR_concordance() |>
      save_plots_structured()
  ),
)

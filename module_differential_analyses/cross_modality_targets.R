rlang::list2(
  targets::tar_target(
    name = significant_elements_modality_distribution_tibble,
    description = "Compare each contrast's share of significant elements across pseudobulk differential modalities",
    command = dplyr::bind_rows(
      gene_expression = significant_elements_tibble.gene_expression,
      chromatin_accessibility = significant_elements_tibble.chromatin_accessibility,
      motif_family_accessibility = significant_elements_tibble.motif_family_accessibility,
      transcription_factor_activity = significant_elements_tibble.transcription_factor_activity,
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
    name = activity_accessibility_source_comparison_tibble.CollecTRI_JASPAR,
    description = "Compare CollecTRI activity, JASPAR motif-family accessibility, and TF expression results [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_comparison_tibble(
      CollecTRI_results_tibble = results_tibble.transcription_factor_activity,
      motif_family_accessibility_results_tibble = results_tibble.motif_family_accessibility,
      gene_expression_results_tibble = results_tibble.gene_expression,
      CollecTRI_JASPAR_family_map = CollecTRI_JASPAR_family_map
    )
  ),
  targets::tar_target(
    name = activity_accessibility_family_comparison_tibble.CollecTRI_JASPAR,
    description = "Collapse CollecTRI-JASPAR comparisons to motif families per model and contrast [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_family_comparison_tibble(
      activity_accessibility_source_comparison_tibble.CollecTRI_JASPAR
    )
  ),
  targets::tar_target(
    name = activity_accessibility_concordance_tibble.CollecTRI_JASPAR,
    description = "Summarize CollecTRI-JASPAR correlation, direction concordance, and joint significance per contrast [part_of_graph:differential_analyses]",
    command = get_CollecTRI_JASPAR_concordance_tibble(
      activity_accessibility_family_comparison_tibble.CollecTRI_JASPAR
    )
  ),
  tarchetypes::tar_file(
    name = activity_accessibility_concordance_plots.CollecTRI_JASPAR,
    description = "Plot expression- versus accessibility-derived TF activity concordance per contrast. [checkpoint:differential_analyses]",
    command = activity_accessibility_concordance_tibble.CollecTRI_JASPAR |>
      plot_CollecTRI_JASPAR_concordance() |>
      save_plots_structured()
  ),
)

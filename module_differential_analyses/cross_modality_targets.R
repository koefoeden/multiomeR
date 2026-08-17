rlang::list2(
  targets::tar_target(
    name = significant_elements_modality_distribution_tibble,
    description = "Compare each contrast's share of significant elements across pseudobulk DX modalities",
    command = dplyr::bind_rows(
      DGE = significant_elements_tibble.DGE,
      DCA = significant_elements_tibble.DCA,
      DTFA = significant_elements_tibble.DTFA,
      DCTA = significant_elements_tibble.DCTA,
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
    description = "Plot each contrast's share of significant elements across pseudobulk DX modalities per model. [checkpoint:differential_analyses]",
    command = significant_elements_modality_distribution_tibble |>
      plot_psbulk_DX_significant_elements_modality_distribution() |>
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
    name = TF_activity_source_comparison_tibble.CollecTRI_DTFA,
    description = "Compare CollecTRI activity, JASPAR motif-family accessibility, and TF expression results [part_of_graph:differential_analyses]",
    command = get_CollecTRI_DTFA_comparison_tibble(
      CollecTRI_results_tibble = results_tibble.DCTA,
      DTFA_results_tibble = results_tibble.DTFA,
      DGE_results_tibble = results_tibble.DGE,
      CollecTRI_JASPAR_family_map = CollecTRI_JASPAR_family_map
    )
  ),
  targets::tar_target(
    name = TF_activity_family_comparison_tibble.CollecTRI_DTFA,
    description = "Collapse CollecTRI-DTFA comparisons to motif families per model and contrast [part_of_graph:differential_analyses]",
    command = get_CollecTRI_DTFA_family_comparison_tibble(
      TF_activity_source_comparison_tibble.CollecTRI_DTFA
    )
  ),
  targets::tar_target(
    name = TF_activity_concordance_tibble.CollecTRI_DTFA,
    description = "Summarize CollecTRI-DTFA correlation, direction concordance, and joint significance per contrast [part_of_graph:differential_analyses]",
    command = get_CollecTRI_DTFA_concordance_tibble(
      TF_activity_family_comparison_tibble.CollecTRI_DTFA
    )
  ),
  tarchetypes::tar_file(
    name = TF_activity_concordance_plots.CollecTRI_DTFA,
    description = "Plot expression- versus accessibility-derived TF activity concordance per contrast. [checkpoint:differential_analyses]",
    command = TF_activity_concordance_tibble.CollecTRI_DTFA |>
      plot_CollecTRI_DTFA_concordance() |>
      save_plots_structured()
  ),
)

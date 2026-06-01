rlang::list2(
  targets::tar_target(
    name = significant_elements_modality_distribution_tibble,
    description = "Compare each contrast's share of significant elements across pseudobulk DX modalities",
    command = dplyr::bind_rows(
      DGE = significant_elements_tibble.DGE,
      DCA = significant_elements_tibble.DCA,
      DTFA = significant_elements_tibble.DTFA,
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
)

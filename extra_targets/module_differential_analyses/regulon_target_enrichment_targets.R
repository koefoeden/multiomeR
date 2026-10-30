rlang::list2(
  targets::tar_target(
    name = results.CollecTRI.regulon_target_enrichment.gene_expression,
    description = "Test CollecTRI regulons as signed gene sets with cameraPR on each gene-expression model's moderated t, with leading targets and TF expression [part_of_graph:differential_analyses]",
    command = get_regulon_target_enrichment_results(feature_matrix_fit.gene_expression, dynamic_tibble.gene_expression,
      results_tibble.gene_expression, filtered_mat_per_model.gene_expression, CollecTRI_human_network_tibble, organism_chr),
    pattern = map(feature_matrix_fit.gene_expression, dynamic_tibble.gene_expression, results_tibble.gene_expression, filtered_mat_per_model.gene_expression)
  ),
  tarchetypes::tar_file(
    name = enrichment_plots.CollecTRI.regulon_target_enrichment.gene_expression,
    description = "Save regulon target enrichment plots per model and contrast to file. [checkpoint:differential_analyses]",
    command = save_plots_structured(
      plots = plot_gene_set_enrichment_results(results.CollecTRI.regulon_target_enrichment.gene_expression, title = "Regulon target enrichment"),
      override_suffix = dynamic_tibble.gene_expression$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(results.CollecTRI.regulon_target_enrichment.gene_expression, dynamic_tibble.gene_expression)
  )
)

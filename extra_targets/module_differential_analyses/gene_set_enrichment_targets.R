rlang::list2(
  targets::tar_target(
    name = gene_sets,
    description = "Load one selected MSigDB gene-set collection for competitive testing",
    command = get_msigdb_gene_sets(
      organism_chr = organism_chr,
      collection_chr = map_MSigDB_collection,
      subcollection_chr = map_MSigDB_subcollection
    ),
    deployment = "main"
  ),
  targets::tar_target(
    name = results,
    description = "Run competitive camera enrichment, with estimated inter-gene correlations, on pseudobulk contrasts for each gene-set subcollection [part_of_graph:differential_analyses]",
    command = get_gene_set_enrichment_results(feature_matrix_fit.gene_expression, dynamic_tibble.gene_expression, gene_sets),
    pattern = map(feature_matrix_fit.gene_expression, dynamic_tibble.gene_expression)
  ),
  tarchetypes::tar_file(
    name = enrichment_plots,
    description = "Save competitive camera plots per model and contrast to file. [checkpoint:differential_analyses]",
    command = save_plots_structured(
      plots = plot_gene_set_enrichment_results(results),
      override_suffix = dynamic_tibble.gene_expression$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(results, dynamic_tibble.gene_expression)
  )
)

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
    description = "Run competitive cameraPR enrichment on pseudobulk contrasts for each gene-set subcollection [part_of_graph:differential_analyses]",
    command = get_gene_set_enrichment_results(
      pseudobulk_feature_matrix_fit = feature_matrix_fit.gene_expression,
      gene_sets = gene_sets,
      pseudobulk_feature_dynamic_tibble = dynamic_tibble.gene_expression
    ),
    pattern = map(dynamic_tibble.gene_expression, feature_matrix_fit.gene_expression)
  ),
  tarchetypes::tar_file(
    name = enrichment_plots,
    description = "Save competitive cameraPR plots per model and contrast to file. [checkpoint:differential_analyses]",
    command = save_plots_structured(
      plots = plot_gene_set_enrichment_results(results),
      override_suffix = dynamic_tibble.gene_expression$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(results, dynamic_tibble.gene_expression)
  )
)

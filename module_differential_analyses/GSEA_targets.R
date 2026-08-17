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
    command = get_GSEA_results(
      psbulk_feature_matrix_fit = feature_matrix_fit.DGE,
      gene_sets = gene_sets,
      psbulk_feature_dynamic_tibble = dynamic_tibble.DGE
    ),
    pattern = map(dynamic_tibble.DGE, feature_matrix_fit.DGE)
  ),
  targets::tar_target(
    name = plots,
    description = "Plot competitive cameraPR results per contrast",
    command = plot_GSEA_results(results),
    pattern = map(results)
  ),
  tarchetypes::tar_file(
    name = plots_file,
    description = "Save competitive cameraPR plots per model and contrast to file. [checkpoint:differential_analyses]",
    command = save_plots_structured(
      plots = plots,
      override_suffix = dynamic_tibble.DGE$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(plots, dynamic_tibble.DGE)
  )
)

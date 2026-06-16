rlang::list2(
  targets::tar_target(
    name = gene_idxs_per_gene_list,
    description = "Map gene indices and annotations for each GSEA gene-set subcollection",
    command = get_gene_idxs_and_annot_per_gene_list(
      DGE_list = feature_matrix_fit.DGE,
      gene_set_subcollection = combined_gene_sets[[map_psbulk_DX_GSEA_subcollection]]
    ),
    pattern = map(feature_matrix_fit.DGE)
  ),
  targets::tar_target(
    name = results,
    description = "Run GSEA on pseudobulk DX models for each gene-set subcollection [part_of_graph:differential_analyses]",
    command = get_GSEA_results(
      psbulk_feature_matrix_fit = feature_matrix_fit.DGE,
      gene_indices_per_gene_list = gene_idxs_per_gene_list,
      psbulk_feature_dynamic_tibble = dynamic_tibble.DGE
    ),
    pattern = map(gene_idxs_per_gene_list, dynamic_tibble.DGE, feature_matrix_fit.DGE)
  ),
  targets::tar_target(
    name = plots,
    description = "Plot GSEA results per contrast",
    command = plot_GSEA_results(results),
    pattern = map(results)
  ),
  tarchetypes::tar_file(
    name = plots_file,
    description = "Save GSEA results plots per model and contrast to file. [checkpoint:differential_analyses]",
    command = save_plots_structured(
      plots = plots,
      override_suffix = dynamic_tibble.DGE$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(plots, dynamic_tibble.DGE)
  )
)

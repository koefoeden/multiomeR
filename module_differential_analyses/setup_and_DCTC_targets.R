rlang::list2(
  tarchetypes::tar_file(
    name = donor_id_metadata_tsv.extended,
    description = "Locate the extended donor ID metadata TSV for differential analyses",
    command = differential_analyses_extended_donor_id_metadata_tsv %||%
      aggregation_donor_id_metadata_tsv,
    deployment = "main"
  ),
  targets::tar_target(
    name = donor_id_metadata_tibble.extended,
    description = "Read the extended donor ID metadata TSV into a tibble [part_of_graph:differential_analyses]",
    command = read_keyed_metadata_tibble(donor_id_metadata_tsv.extended, "donor_id")
  ),
  targets::tar_target(
    name = models,
    description = "Build named list of configured pseudobulk DX model specifications [part_of_graph:differential_analyses]",
    command = normalize_psbulk_feature_models(differential_analyses_psbulk_DX_models)
  ),
  tarchetypes::tar_file(
    name = pseudobulk_depth_distribution_plot,
    description = "Plot pseudobulk count depths and detected features per cluster-donor sample. [checkpoint:differential_analyses]",
    command = dplyr::bind_rows(
      pseudobulk_depth_tibble.GEX,
      pseudobulk_depth_tibble.ATAC
    ) |>
      plot_pseudobulk_depth_distribution(
        min_ATAC_sample_counts = differential_analyses_DTFA_min_ATAC_counts
      ) |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = plot_phenotype_vars.DCTC,
    description = "Assert and load the phenotype variable names for differential cell-type composition plots",
    command = assert_cfg_is_set(differential_analyses_DCTC_plot_phenotype_vars, "differential_analyses_DCTC_plot_phenotype_vars")
  ),
  targets::tar_target(
    name = donor_id_metadata_tibble.analysis,
    description = "Project extended donor metadata to aggregation donors and configured differential-analysis variables [part_of_graph:differential_analyses]",
    command = {
      aggregation_donor_ids <- metadata_w_cell_types_tibble.WNN |>
        dplyr::distinct(donor_id) |>
        dplyr::pull(donor_id) |>
        sort()
      aggregation_donor_metadata <- subset_keyed_metadata_tibble(
        donor_id_metadata_tibble.extended,
        "donor_id",
        aggregation_donor_ids,
        donor_id_metadata_tsv.extended
      )
      project_keyed_metadata_tibble(
        aggregation_donor_metadata,
        "donor_id",
        get_differential_analysis_metadata_columns(
          models = models,
          DCTC_formula_chr = differential_analyses_DCTC_formula_chr,
          DCTC_plot_phenotype_vars = plot_phenotype_vars.DCTC,
          DCTC_color_by_categorical_metadata_column = differential_analyses_DCTC_color_by_categorical_metadata_column
        ),
        strict = TRUE,
        canonical = TRUE
      )
    }
  ),
  tarchetypes::tar_file(
    name = by_phenotype_per_cluster_plots.DCTC,
    description = "Save DCTC scatter plot objects per cluster for each phenotype variable. [checkpoint:differential_analyses]",
    command = {
      metadata_tibble <- metadata_w_cell_types_tibble.WNN |>
        dplyr::left_join(donor_id_metadata_tibble.analysis)

      plot_DCTC_by_phenotype_per_cluster(
        metadata_tibble = metadata_tibble,
        DCTC_plot_phenotype_vars = plot_phenotype_vars.DCTC,
        cluster_col = "PCA_harmony_SNN_cluster_cell_type",
        DCTC_color_by_categorical_metadata_column = differential_analyses_DCTC_color_by_categorical_metadata_column
      ) |>
        save_plots_structured()
    },
    pattern = map(plot_phenotype_vars.DCTC)
  ),
  targets::tar_target(
    name = model_results.DCTC,
    description = "Fit the configured differential cell-type composition model [part_of_graph:differential_analyses]",
    command = metadata_w_cell_types_tibble.WNN |>
      get_DCTC_model_results(
        extended_donor_id_metadata_tibble = donor_id_metadata_tibble.analysis,
        cluster_col = "PCA_harmony_SNN_cluster_cell_type",
        DCTC_formula_chr = differential_analyses_DCTC_formula_chr
      )
  ),
  tarchetypes::tar_file(
    name = model_change_per_unit_plot.DCTC,
    description = "Plot modelled change in cell-type composition per unit of each phenotype and save to file. [checkpoint:differential_analyses]",
    command = model_results.DCTC |>
      plot_DCTC_model_change_per_unit() |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = model_coefs_forest.DCTC,
    description = "Plot forest plot of DCTC model coefficients per cell type and save to file. [checkpoint:differential_analyses]",
    command = model_results.DCTC |>
      plot_DCTC_model_coefs_forest() |>
      save_plots_structured()
  )
)

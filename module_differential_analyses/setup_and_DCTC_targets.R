rlang::list2(
  tarchetypes::tar_file(
    name = donor_id_metadata_tsv.extended,
    description = "Locate the extended donor metadata for differential analyses",
    command = differential_analyses_extended_donor_id_metadata_tsv %||% aggregation_donor_id_metadata_tsv,
    deployment = "main"
  ),
  targets::tar_target(
    name = donor_id_metadata_tibble.extended,
    description = "Read keyed differential-analysis donor metadata [part_of_graph:differential_analyses]",
    command = read_keyed_metadata_tibble(donor_id_metadata_tsv.extended, "donor_id")
  ),
  targets::tar_target(
    name = models,
    description = "Normalize named feature models [part_of_graph:differential_analyses]",
    command = normalize_psbulk_feature_models(differential_analyses_psbulk_DX_models)
  ),
  targets::tar_target(
    name = models.DCTC,
    description = "Normalize named abundance models [part_of_graph:differential_analyses]",
    command = normalize_differential_models(differential_analyses_DCTC_models)
  ),
  targets::tar_target(
    name = donor_id_metadata_tibble.analysis,
    description = "Project metadata to aggregation donors and all configured model variables [part_of_graph:differential_analyses]",
    command = project_keyed_metadata_tibble(
      subset_keyed_metadata_tibble(donor_id_metadata_tibble.extended, "donor_id",
        sort(unique(metadata_w_cell_types_tibble.WNN$donor_id)), donor_id_metadata_tsv.extended),
      "donor_id", get_differential_analysis_metadata_columns(models, models.DCTC),
      strict = TRUE, canonical = TRUE
    )
  ),
  tarchetypes::tar_file(
    name = pseudobulk_depth_distribution_plot,
    description = "Plot depth and detected features per cluster-donor sample. [checkpoint:differential_analyses]",
    command = dplyr::bind_rows(pseudobulk_depth_tibble.GEX, pseudobulk_depth_tibble.ATAC) |>
      plot_pseudobulk_depth_distribution(min_ATAC_sample_counts = differential_analyses_DTFA_min_ATAC_counts) |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = dynamic_tibble.DCTC,
    description = "Branch over named abundance models",
    command = tibble::enframe(models.DCTC, name = "model_name", value = "model") |>
      dplyr::arrange(model_name),
    iteration = "vector"
  ),
  targets::tar_target(
    name = model_data.DCTC,
    description = "Select model donors and wells; construct complete cell-type counts with fixed denominators [part_of_graph:differential_analyses]",
    command = prepare_DCTC_model_data(metadata_w_cell_types_tibble.WNN,
      donor_id_metadata_tibble.analysis, dynamic_tibble.DCTC$model[[1]], "PCA_harmony_SNN_cluster_cell_type"),
    pattern = map(dynamic_tibble.DCTC),
    iteration = "list"
  ),
  targets::tar_target(
    name = model_results.DCTC,
    description = "Fit separate beta-binomial abundance models and named contrasts [part_of_graph:differential_analyses]",
    command = fit_DCTC_model(model_data.DCTC, dynamic_tibble.DCTC$model[[1]], dynamic_tibble.DCTC$model_name),
    pattern = map(model_data.DCTC, dynamic_tibble.DCTC)
  ),
  tarchetypes::tar_file(
    name = cohort_tsv.DCTC,
    description = "Export donor eligibility, exclusions, selected wells and denominators by abundance model",
    command = save_differential_model_table(model_data.DCTC$cohort, dynamic_tibble.DCTC$model_name),
    pattern = map(model_data.DCTC, dynamic_tibble.DCTC)
  ),
  tarchetypes::tar_file(
    name = counts_tsv.DCTC,
    description = "Export the exact donor-cell-type counts used in abundance fits and plots",
    command = save_differential_model_table(model_data.DCTC$counts, dynamic_tibble.DCTC$model_name),
    pattern = map(model_data.DCTC, dynamic_tibble.DCTC)
  ),
  tarchetypes::tar_file(
    name = results_tsv.DCTC,
    description = "Export named abundance contrasts, BH FDR and fit diagnostics",
    command = save_differential_model_table(model_results.DCTC, dynamic_tibble.DCTC$model_name),
    pattern = map(model_results.DCTC, dynamic_tibble.DCTC)
  ),
  tarchetypes::tar_file(
    name = model_plots.DCTC,
    description = "Plot observed donor proportions and abundance contrasts per model. [checkpoint:differential_analyses]",
    command = plot_DCTC_model(model_data.DCTC, dynamic_tibble.DCTC$model[[1]],
      dynamic_tibble.DCTC$model_name, model_results.DCTC) |>
      save_plots_structured(override_suffix = dynamic_tibble.DCTC$model_name, width = 12, height = 10),
    pattern = map(model_data.DCTC, model_results.DCTC, dynamic_tibble.DCTC)
  )
)

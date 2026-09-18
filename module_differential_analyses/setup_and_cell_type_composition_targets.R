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
    command = normalize_pseudobulk_feature_models(differential_analyses_pseudobulk_models)
  ),
  targets::tar_target(
    name = models.cell_type_composition,
    description = "Normalize named abundance models [part_of_graph:differential_analyses]",
    command = normalize_differential_models(differential_analyses_cell_type_composition_models)
  ),
  targets::tar_target(
    name = donor_id_metadata_tibble.analysis,
    description = "Project metadata to aggregation donors and all configured model variables [part_of_graph:differential_analyses]",
    command = project_keyed_metadata_tibble(
      subset_keyed_metadata_tibble(donor_id_metadata_tibble.extended, "donor_id",
        sort(unique(metadata_w_cell_types_tibble.WNN$donor_id)), donor_id_metadata_tsv.extended),
      "donor_id", get_differential_analysis_metadata_columns(models, models.cell_type_composition),
      strict = TRUE, canonical = TRUE
    )
  ),
  tarchetypes::tar_file(
    name = pseudobulk_depth_distribution_plot,
    description = "Plot depth and detected features per cluster-donor sample. [checkpoint:differential_analyses]",
    command = dplyr::bind_rows(pseudobulk_depth_tibble.GEX, pseudobulk_depth_tibble.ATAC) |>
      plot_pseudobulk_depth_distribution(min_ATAC_sample_counts = differential_analyses_motif_family_accessibility_min_ATAC_counts) |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = dynamic_tibble.cell_type_composition,
    description = "Branch over named abundance models",
    command = tibble::enframe(models.cell_type_composition, name = "model_name", value = "model") |>
      dplyr::arrange(model_name),
    iteration = "vector"
  ),
  targets::tar_target(
    name = model_data.cell_type_composition,
    description = "Select model donors and wells; construct complete cell-type counts with fixed denominators [part_of_graph:differential_analyses]",
    command = prepare_cell_type_composition_model_data(metadata_w_cell_types_tibble.WNN,
      donor_id_metadata_tibble.analysis, dynamic_tibble.cell_type_composition$model[[1]], "WNN_harmony_SNN_cluster_cell_type"),
    pattern = map(dynamic_tibble.cell_type_composition),
    iteration = "list"
  ),
  targets::tar_target(
    name = model_results.cell_type_composition,
    description = "Fit separate beta-binomial abundance models and named contrasts [part_of_graph:differential_analyses]",
    command = fit_cell_type_composition_model(model_data.cell_type_composition, dynamic_tibble.cell_type_composition$model[[1]], dynamic_tibble.cell_type_composition$model_name),
    pattern = map(model_data.cell_type_composition, dynamic_tibble.cell_type_composition)
  ),
  tarchetypes::tar_file(
    name = cohort_tsv.cell_type_composition,
    description = "Export donor eligibility, exclusions, selected wells and denominators by abundance model",
    command = save_differential_model_table(model_data.cell_type_composition$cohort, dynamic_tibble.cell_type_composition$model_name),
    pattern = map(model_data.cell_type_composition, dynamic_tibble.cell_type_composition)
  ),
  tarchetypes::tar_file(
    name = counts_tsv.cell_type_composition,
    description = "Export the exact donor-cell-type counts used in abundance fits and plots",
    command = save_differential_model_table(model_data.cell_type_composition$counts, dynamic_tibble.cell_type_composition$model_name),
    pattern = map(model_data.cell_type_composition, dynamic_tibble.cell_type_composition)
  ),
  tarchetypes::tar_file(
    name = results_tsv.cell_type_composition,
    description = "Export named abundance contrasts, BH FDR and fit diagnostics",
    command = save_differential_model_table(model_results.cell_type_composition, dynamic_tibble.cell_type_composition$model_name),
    pattern = map(model_results.cell_type_composition, dynamic_tibble.cell_type_composition)
  ),
  tarchetypes::tar_file(
    name = model_plots.cell_type_composition,
    description = "Plot observed donor proportions and abundance contrasts per model. [checkpoint:differential_analyses]",
    command = plot_cell_type_composition_model(model_data.cell_type_composition, dynamic_tibble.cell_type_composition$model[[1]],
      dynamic_tibble.cell_type_composition$model_name, model_results.cell_type_composition) |>
      save_plots_structured(override_suffix = dynamic_tibble.cell_type_composition$model_name, width = 12, height = 10),
    pattern = map(model_data.cell_type_composition, model_results.cell_type_composition, dynamic_tibble.cell_type_composition)
  )
)

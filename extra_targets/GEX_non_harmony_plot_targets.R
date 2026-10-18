# Plots on the uncorrected GEX UMAP. _targets.R maps them only over aggregations with
# Harmony covariates; without them the non-Harmony UMAP repeats the Harmony one.
rlang::list2(
  targets::tar_target(
    name = non_harmony_continuous_UMAP_spec.GEX,
    description = "Continuous non-Harmony GEX UMAP variables to plot one at a time",
    command = {
      feature_vars <- intersect(
        unique(c(GEX_marker_genes_vec, interesting_genes)),
        rownames(aggregated_counts_BPCells_matrix.GEX)
      )
      tibble::tibble(
        variable = c(aggregation_non_peak_based_continuous_vars, feature_vars),
        value_source = c(
          rep("metadata", length(aggregation_non_peak_based_continuous_vars)),
          rep("feature", length(feature_vars))
        )
      )
    },
    iteration = "vector"
  ),
  tarchetypes::tar_file(
    name = non_harmony.categorical.UMAPs.3_GEX_QC,
    description = "UMAPs colored by categorical metadata variables on the uncorrected PCA embedding. [checkpoint:3_GEX-QC]",
    command = metadata_w_cell_types_analysis_tibble.GEX |>
      plot_UMAP_from_metadata(variable = categorical_UMAP_var.GEX, umap_cols = c("GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(categorical_UMAP_var.GEX, "[/\\\\]", "_")
      ),
    pattern = map(categorical_UMAP_var.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = non_harmony.continuous.UMAPs.3_GEX_QC,
    description = "UMAPs colored by continuous features on the uncorrected PCA embedding. [checkpoint:3_GEX-QC]",
    command = plot_UMAP_from_metadata(
      metadata_tibble = continuous_UMAP_inputs.GEX$metadata_tibble,
      variable = non_harmony_continuous_UMAP_spec.GEX$variable,
      value_source = non_harmony_continuous_UMAP_spec.GEX$value_source,
      feature_matrix = continuous_UMAP_inputs.GEX$feature_matrix,
      umap_cols = c("GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")
    ) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(non_harmony_continuous_UMAP_spec.GEX$variable, "[/\\\\]", "_")
      ),
    pattern = map(non_harmony_continuous_UMAP_spec.GEX)
  )
)

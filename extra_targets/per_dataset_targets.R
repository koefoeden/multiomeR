rlang::list2(
  targets::tar_target(
    name = per_dataset_unfiltered_cells_n,
    description = "Sum full QC metadata barcode counts across all GEM wells in this dataset",
    command = sum(unlist(dataset_unfiltered_cells_n_vecs_syms)),
  ),
  targets::tar_target(
    name = per_dataset_excluded_BCs_list,
    description = "Combine per GEM well QC exclusion lists into a single dataset-level list",
    command = collapse_duplicate_names(purrr::flatten(
      dataset_excluded_barcodes_by_type_list_syms
    ))
  ),
  tarchetypes::tar_file(
    name = per_dataset_excluded_upset,
    description = "Plot an UpSet plot of dataset-level QC exclusion overlaps across all GEM wells and save to file",
    command = plot_upset_from_excluded_BCs_list(
      QC_excluded_BCs_list = per_dataset_excluded_BCs_list,
      n_total = per_dataset_unfiltered_cells_n
    ) |>
      save_plots_structured(width = 20, height = 7),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = per_dataset_excluded_cellranger_only_BCs_list,
    description = "Combine per GEM well CellRanger-only QC exclusion lists into a dataset-level list",
    command = collapse_duplicate_names(purrr::flatten(
      dataset_excluded_cellranger_only_barcodes_by_type_list_syms
    ))
  ),
  tarchetypes::tar_file(
    name = per_dataset_excluded_cellranger_only_upset,
    description = "Plot an UpSet plot of dataset-level CellRanger-only QC exclusion overlaps and save to file",
    command = plot_upset_from_excluded_BCs_list(
      QC_excluded_BCs_list = per_dataset_excluded_cellranger_only_BCs_list,
      n_total = nrow(per_dataset_cellranger_kept_metadata_tibble)
    ) |>
      save_plots_structured(width = 20, height = 7),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = per_dataset_cellranger_kept_metadata_tibble,
    description = "Combine CellRanger-kept metadata across GEM wells, selecting QC variables for per-dataset QC plots",
    command = dplyr::bind_rows(dataset_cellranger_kept_metadata_tibble_syms) |>
      dplyr::select(dplyr::any_of(c("GEM_well_ID", "dataset", PROCESSING_PER_GEM_well_QC_VARS))),
  ),
  tarchetypes::tar_file(
    name = per_dataset_QC_violins,
    description = "Plot violin plots of QC metrics per GEM well for this dataset and save to file",
    command = {
      threshold_tibble <- get_QC_exclude_threshold_tibble(
        QC_exclude_vector = dataset_QC_exclude_list_per_GEM_well,
        feature_names = PROCESSING_PER_GEM_well_QC_VARS
      )

      plot_data <- per_dataset_cellranger_kept_metadata_tibble |>
        tidyr::pivot_longer(
          cols = dplyr::any_of(PROCESSING_PER_GEM_well_QC_VARS),
          names_to = "feature",
          values_to = "value"
        ) |>
        dplyr::filter(
          .by = c(GEM_well_ID, feature),
          value >= stats::quantile(value, probs = 0.02, na.rm = TRUE),
          value <= stats::quantile(value, probs = 0.98, na.rm = TRUE)
        )

      plots <- plot_data$feature |>
        unique() |>
        purrr::set_names() |>
        purrr::map(\(feature) {
          feature_thresholds <- threshold_tibble |>
            dplyr::filter(.data$feature == .env$feature)

          plot <- plot_data |>
            dplyr::filter(.data$feature == .env$feature) |>
            ggplot2::ggplot(ggplot2::aes(x = GEM_well_ID, y = value, fill = dataset))

          if (nrow(feature_thresholds) > 0) {
            plot <- plot +
              ggplot2::geom_rect(
                data = feature_thresholds,
                ggplot2::aes(ymin = ymin, ymax = ymax),
                xmin = -Inf,
                xmax = Inf,
                fill = "grey60",
                alpha = 0.25,
                inherit.aes = FALSE
              ) +
              ggplot2::geom_hline(
                data = feature_thresholds,
                ggplot2::aes(yintercept = threshold),
                color = "grey35",
                linetype = "dashed",
                inherit.aes = FALSE
              )
          }

          plot +
            ggplot2::geom_violin(scale = "width") +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
              axis.title.x = ggplot2::element_blank(),
              legend.position = "none"
            )
        })

      save_plots_structured(plots)
    }
  )
)

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
    command = plot_per_dataset_QC_violins(
      metadata_tibble = per_dataset_cellranger_kept_metadata_tibble,
      QC_exclude_vector = dataset_QC_exclude_list,
      feature_names = PROCESSING_PER_GEM_well_QC_VARS
    ) |>
      save_plots_structured()
  )
)

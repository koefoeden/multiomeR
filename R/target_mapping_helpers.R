#' Validate processing and aggregation config
#'
#' Check that aggregation YAML rows reference valid, non-empty GEM well sets.
#'
#' @param GEM_well_tibble GEM well config tibble containing the valid `GEM_well_ID`
#'   values.
#' @param aggregation_tibble_from_yaml Aggregation config tibble with
#'   `aggregation` labels and list-column `aggregation_GEM_well_IDs`.
#' @param GEM_well_config_file Path to the GEM well TSV config.
#' @param aggregation_config_file Path to the aggregation YAML config.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

validate_processing_and_aggregation_config <- function(
  GEM_well_tibble,
  aggregation_tibble_from_yaml,
  GEM_well_config_file = "cfg_GEM_wells.tsv",
  aggregation_config_file = "cfg_aggregations.yaml"
) {
  configured_aggregation_GEM_well_IDs <- unique(unlist(
    aggregation_tibble_from_yaml$aggregation_GEM_well_IDs
  ))

  unknown_aggregation_GEM_wells_vec <- base::setdiff(
    configured_aggregation_GEM_well_IDs,
    GEM_well_tibble$GEM_well_ID
  )
  if (length(unknown_aggregation_GEM_wells_vec) > 0) {
    stop(
      aggregation_config_file,
      " references GEM_well_ID value(s) not defined in ",
      GEM_well_config_file,
      ": ",
      paste(unknown_aggregation_GEM_wells_vec, collapse = ", "),
      call. = FALSE
    )
  }

  empty_aggregation_vec <- aggregation_tibble_from_yaml$aggregation[
    lengths(aggregation_tibble_from_yaml$aggregation_GEM_well_IDs) == 0
  ]
  if (length(empty_aggregation_vec) > 0) {
    stop(
      aggregation_config_file,
      " contains aggregation(s) without GEM_well_ID values: ",
      paste(empty_aggregation_vec, collapse = ", "),
      call. = FALSE
    )
  }

  inactive_GEM_well_tibble <- GEM_well_tibble |>
    dplyr::filter(!purrr::map_lgl(is_active, isTRUE)) |>
    dplyr::select(GEM_well_ID, dataset)
  inactive_reference_details <- aggregation_tibble_from_yaml |>
    dplyr::select(aggregation, aggregation_GEM_well_IDs) |>
    tidyr::unnest_longer(aggregation_GEM_well_IDs, values_to = "GEM_well_ID") |>
    dplyr::inner_join(inactive_GEM_well_tibble, by = "GEM_well_ID") |>
    dplyr::summarise(
      GEM_well_IDs = paste(GEM_well_ID, collapse = ", "),
      .by = c(aggregation, dataset)
    ) |>
    dplyr::transmute(
      details = paste0(
        "aggregation '", aggregation,
        "' -> inactive GEM well(s) in dataset '", dataset,
        "' (GEM_well_ID value(s): ", GEM_well_IDs, ")"
      )
    ) |>
    dplyr::pull(details)
  if (length(inactive_reference_details) > 0) {
    stop(
      aggregation_config_file,
      " contains active aggregation references to inactive GEM wells: ",
      paste(inactive_reference_details, collapse = "; "),
      ". Activate each GEM well in cfg_GEM_wells.tsv or remove it from the active aggregation.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

read_keyed_metadata_tibble <- function(metadata_tsv, key_col) {
  metadata_tibble <- readr::read_tsv(metadata_tsv, show_col_types = FALSE)

  if (!key_col %in% colnames(metadata_tibble)) {
    stop("Metadata file must contain key column '", key_col, "': ", metadata_tsv, call. = FALSE)
  }

  duplicate_keys <- metadata_tibble |>
    dplyr::count(.data[[key_col]]) |>
    dplyr::filter(.data$n > 1) |>
    dplyr::pull(.data[[key_col]])
  if (length(duplicate_keys) > 0) {
    stop(
      metadata_tsv,
      " contains duplicated ",
      key_col,
      " value(s): ",
      paste(duplicate_keys, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_tibble |>
    dplyr::mutate(dplyr::across(dplyr::all_of(key_col), as.character))
}

assert_donor_GEM_well_metadata_column_ownership <- function(donor_id_metadata_tibble, GEM_well_metadata_tibble) {
  overlapping_cols <- intersect(
    setdiff(colnames(donor_id_metadata_tibble), "donor_id"),
    setdiff(colnames(GEM_well_metadata_tibble), "GEM_well_ID")
  )

  if (length(overlapping_cols) > 0) {
    stop(
      "Donor and GEM well metadata contain overlapping non-key column(s): ",
      paste(overlapping_cols, collapse = ", "),
      ". Each metadata column must belong to exactly one table.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Build GEM well mapping tibble
#'
#' Read GEM well-level processing configuration for the GEM well `tar_map()`.
#'
#' The ` ;; ` separator preserves multiple QC expressions inside the tabular
#' configuration file. User-facing GEM well column names are normalized to the
#' existing internal target-command symbols so the configuration migration does
#' not invalidate cached processing targets.
#'
#' @param GEM_well_config_file Path to the GEM well TSV config.
#' @return A tibble with one row per GEM well and normalized config columns.
#' @keywords internal

parse_GEM_well_QC_exclude_list <- function(value, GEM_well_ID) {
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    return(NULL)
  }

  parsed <- strsplit(value, ";;", fixed = TRUE)[[1]] |>
    trimws() |>
    purrr::discard(\(item) !nzchar(item))
  purrr::walk(parsed, \(expression) {
    tryCatch(
      rlang::parse_expr(expression),
      error = function(error) {
        stop(
          "Invalid filter expression for GEM well '", GEM_well_ID, "': ",
          expression, " (", conditionMessage(error), ")",
          call. = FALSE
        )
      }
    )
  })
  parsed
}

build_GEM_well_tibble <- function(GEM_well_config_file = "cfg_GEM_wells.tsv") {
  required_columns <- c(
    "GEM_well_ID",
    "dataset",
    "GEM_well_donor_id",
    "GEM_well_n_donors",
    "GEM_well_cellranger_arc_count_dir",
    "GEM_well_cellbender_h5_file",
    "GEM_well_donors_VCF_file",
    "GEM_well_cellranger_arc_reference_json",
    "GEM_well_QC_exclude_list",
    "GEM_well_run_amulet",
    "GEM_well_is_active"
  )
  GEM_well_tibble <- readr::read_tsv(GEM_well_config_file, show_col_types = FALSE) |>
    dplyr::mutate(dplyr::across(c(GEM_well_ID, GEM_well_donor_id), as.character))

  missing_columns <- setdiff(required_columns, colnames(GEM_well_tibble))
  if (length(missing_columns) > 0L) {
    stop(
      GEM_well_config_file, " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  duplicated_GEM_well_IDs <- GEM_well_tibble |>
    dplyr::count(.data$GEM_well_ID) |>
    dplyr::filter(.data$n > 1L) |>
    dplyr::pull(.data$GEM_well_ID)
  if (length(duplicated_GEM_well_IDs) > 0L) {
    stop(
      GEM_well_config_file, " contains duplicated GEM_well_ID value(s): ",
      paste(duplicated_GEM_well_IDs, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(GEM_well_tibble$GEM_well_cellranger_arc_reference_json) ||
      any(!nzchar(GEM_well_tibble$GEM_well_cellranger_arc_reference_json))) {
    stop("Every GEM well must define GEM_well_cellranger_arc_reference_json.", call. = FALSE)
  }
  if (anyNA(GEM_well_tibble$GEM_well_run_amulet) ||
      anyNA(GEM_well_tibble$GEM_well_is_active)) {
    stop("Every GEM well must define GEM_well_run_amulet and GEM_well_is_active.", call. = FALSE)
  }

  GEM_well_tibble |>
    dplyr::mutate(
      dataset_cellranger_arc_reference_json = purrr::map(
        .data$GEM_well_cellranger_arc_reference_json,
        identity
      ),
      dataset_QC_exclude_list_per_GEM_well = purrr::map2(
        .data$GEM_well_QC_exclude_list,
        .data$GEM_well_ID,
        parse_GEM_well_QC_exclude_list
      ),
      dataset_run_amulet = purrr::map(.data$GEM_well_run_amulet, identity),
      is_active = purrr::map(.data$GEM_well_is_active, identity)
    ) |>
    dplyr::select(
      dplyr::all_of(c(
        "GEM_well_ID",
        "dataset",
        "GEM_well_donor_id",
        "GEM_well_n_donors",
        "GEM_well_cellranger_arc_count_dir",
        "GEM_well_cellbender_h5_file",
        "GEM_well_donors_VCF_file",
        "dataset_cellranger_arc_reference_json",
        "dataset_QC_exclude_list_per_GEM_well",
        "dataset_run_amulet",
        "is_active"
      ))
    )
}

#' Build the active GEM well mapping tibble
#'
#' Keep only GEM wells whose GEM well-level config enables graph construction.
#'
#' @param GEM_well_tibble Complete GEM well mapping tibble created by
#'   `build_GEM_well_tibble()`.
#' @return The active GEM well rows, with columns unchanged.
#' @keywords internal

build_active_GEM_well_tibble <- function(GEM_well_tibble) {
  GEM_well_tibble |>
    dplyr::filter(purrr::map_lgl(is_active, isTRUE))
}

#' Build aggregation mapping tibble
#'
#' Filter active aggregations, validate their GEM well references, and add the
#' list-columns used to splice upstream target symbols into aggregation targets.
#'
#' @param aggregation_tibble_all_from_yaml Aggregation config tibble created
#'   from `cfg_aggregations.yaml`.
#' @param GEM_well_tibble GEM well mapping tibble created by
#'   `build_GEM_well_tibble()`.
#' @param aggregation_config_file Path to the aggregation YAML config.
#' @param GEM_well_config_file Path to the GEM well TSV config.
#' @return A tibble with one row per active aggregation, QC feature columns, and
#'   target symbol list-columns.
#' @keywords internal

build_aggregation_tibble <- function(
  aggregation_tibble_all_from_yaml,
  GEM_well_tibble,
  aggregation_config_file = "cfg_aggregations.yaml",
  GEM_well_config_file = "cfg_GEM_wells.tsv"
) {
  aggregation_tibble_from_yaml <- aggregation_tibble_all_from_yaml |>
    dplyr::filter(purrr::map_lgl(is_active, isTRUE))

  aggregation_tibble <- aggregation_tibble_from_yaml |>
    add_target_sym_cols(
      aggregation_GEX_counts_BPCells_matrix_syms = target_sym_col("GEX_counts_BPCells_matrix", "aggregation_GEM_well_IDs"),
      aggregation_fragments_w_prefix_bpcells_syms = target_sym_col("fragments_w_prefix_bpcells", "aggregation_GEM_well_IDs"),
      aggregation_cellranger_summary_file_syms = target_sym_col("cellranger_summary_file", "aggregation_GEM_well_IDs"),
      aggregation_cellranger_kept_metadata_tibble_syms = target_sym_col("cellranger_kept_metadata_tibble", "aggregation_GEM_well_IDs"),
      aggregation_unfiltered_cells_n_vecs_syms = target_sym_col("unfiltered_cells_n_vecs", "aggregation_GEM_well_IDs"),
      aggregation_excluded_barcodes_by_type_list_syms = target_sym_col("excluded_barcodes_by_type_list", "aggregation_GEM_well_IDs"),
      aggregation_excluded_cellranger_only_barcodes_by_type_list_syms = target_sym_col("excluded_cellranger_only_barcodes_by_type_list", "aggregation_GEM_well_IDs"),
      aggregation_cellranger_ref_list_syms = target_sym_col("cellranger_ref_list", "aggregation_GEM_well_IDs"),
      aggregation_gene_features_df_syms = target_sym_col("gene_features_df", "aggregation_GEM_well_IDs", transform = \(ids) utils::head(ids, 1))
    )

  project_categorical_vars <- purrr::map(
    aggregation_tibble$aggregation_categorical_vars,
    \(aggregation_categorical_vars) {
      c(
        aggregation_categorical_vars,
        "GEM_well_ID",
        "donor_id"
      )
    }
  )
  non_peak_vars <- purrr::map(
    aggregation_tibble$aggregation_continuous_vars,
    \(aggregation_continuous_vars) {
      c(
        aggregation_continuous_vars,
        PROCESSING_CONTINUOUS_QC_FEATURES$non_peak_QC
      )
    }
  )
  peak_vars <- purrr::map(
    non_peak_vars,
    \(feature_names) c(PROCESSING_CONTINUOUS_QC_FEATURES$peak_QC, feature_names)
  )

  aggregation_tibble <- dplyr::bind_cols(
    aggregation_tibble,
    PROCESSING_QC_FEATURE_SETS$categorical |>
      purrr::map(
        \(feature_names) {
          purrr::map(
            project_categorical_vars,
            \(project_vars) c(feature_names, project_vars)
          )
        }
      ) |>
      tibble::as_tibble(),
    tibble::tibble(
      aggregation_non_peak_based_continuous_QC_vars = rep(
        list(PROCESSING_CONTINUOUS_QC_FEATURES$non_peak_QC),
        nrow(aggregation_tibble)
      ),
      aggregation_non_peak_based_continuous_vars = non_peak_vars,
      aggregation_peak_based_continuous_QC_vars = rep(
        list(PROCESSING_CONTINUOUS_QC_FEATURES$peak_QC),
        nrow(aggregation_tibble)
      ),
      aggregation_w_peaks_continuous_vars = peak_vars,
      aggregation_w_WNN_continuous_vars = purrr::map(
        peak_vars,
        \(feature_names) c(feature_names, PROCESSING_CONTINUOUS_QC_FEATURES$WNN)
      ),
      aggregation_continuous_features_vec.GEX = purrr::map(
        non_peak_vars,
        \(feature_names) c(PROCESSING_CONTINUOUS_QC_FEATURES$GEX_cell_cycle, feature_names)
      )
    )
  )

  validate_processing_and_aggregation_config(
    GEM_well_tibble = GEM_well_tibble,
    aggregation_tibble_from_yaml = aggregation_tibble_from_yaml,
    GEM_well_config_file = GEM_well_config_file,
    aggregation_config_file = aggregation_config_file
  )

  aggregation_tibble |>
    dplyr::mutate(
      aggregation_dataset_vec = purrr::map(
        aggregation_GEM_well_IDs,
        \(ids) GEM_well_tibble$dataset[match(ids, GEM_well_tibble$GEM_well_ID)] |>
          unique()
      )
    ) |>
    add_target_sym_cols(
      per_dataset_excluded_upset_syms = target_sym_col("per_dataset_excluded_upset", "aggregation_dataset_vec"),
      per_dataset_excluded_cellranger_only_upset_syms = target_sym_col("per_dataset_excluded_cellranger_only_upset", "aggregation_dataset_vec"),
      per_dataset_QC_violins_syms = target_sym_col("per_dataset_QC_violins", "aggregation_dataset_vec")
    )
}

#' Build dataset config tibble
#'
#' Collapse GEM well-level settings by dataset for interactive configuration
#' loading and dataset-level QC summaries.
#'
#' @param GEM_well_tibble GEM well mapping tibble created by
#'   `build_GEM_well_tibble()`.
#' @return A tibble with one row per dataset and collapsed GEM well settings.
#' @keywords internal

build_dataset_config_tibble <- function(GEM_well_tibble) {
  GEM_well_tibble |>
    dplyr::summarise(
      dataset_cellranger_arc_reference_json = list(
        .data$dataset_cellranger_arc_reference_json[[1]]
      ),
      dataset_QC_exclude_list_per_GEM_well = list(unique(unlist(
        .data$dataset_QC_exclude_list_per_GEM_well,
        use.names = FALSE
      ))),
      dataset_run_amulet = list(any(unlist(.data$dataset_run_amulet))),
      is_active = list(any(unlist(.data$is_active))),
      .by = dataset
    )
}

#' Build dataset mapping tibble
#'
#' Collapse GEM well rows by dataset for the dataset QC-summary `tar_map()`.
#'
#' @inheritParams build_dataset_config_tibble
#' @return A tibble with one row per dataset and GEM well target symbol
#'   list-columns.
#' @keywords internal

build_dataset_tibble <- function(GEM_well_tibble) {
  dataset_config_tibble <- build_dataset_config_tibble(GEM_well_tibble)

  GEM_well_tibble |>
    dplyr::summarise(
      dataset_GEM_well_IDs = list(GEM_well_ID),
      .by = dataset
    ) |>
    add_target_sym_cols(
      dataset_unfiltered_cells_n_vecs_syms = target_sym_col("unfiltered_cells_n_vecs", "dataset_GEM_well_IDs"),
      dataset_cellranger_kept_metadata_tibble_syms = target_sym_col("cellranger_kept_metadata_tibble", "dataset_GEM_well_IDs"),
      dataset_excluded_cellranger_only_barcodes_by_type_list_syms = target_sym_col("excluded_cellranger_only_barcodes_by_type_list", "dataset_GEM_well_IDs"),
      dataset_excluded_barcodes_by_type_list_syms = target_sym_col("excluded_barcodes_by_type_list", "dataset_GEM_well_IDs")
    ) |>
    dplyr::left_join(dataset_config_tibble, by = "dataset")
}

#' Get Roadmap EDACC names
#'
#' Extract the unique non-empty Roadmap EDACC names requested across active
#' aggregations.
#'
#' @param aggregation_tibble Aggregation mapping tibble created by
#'   `build_aggregation_tibble()`.
#' @return A named character vector of Roadmap EDACC names.
#' @keywords internal

get_roadmap_EDACC_names <- function(aggregation_tibble) {
  aggregation_tibble$aggregation_roadmap_EDACC_names |>
    unlist(use.names = FALSE) |>
    as.character() |>
    purrr::discard(is.na) |>
    purrr::discard(\(x) x == "") |>
    unique() |>
    purrr::set_names()
}

#' Normalize module names
#'
#' Convert a possibly nested or missing module config value into a clean
#' character vector.
#'
#' @param modules Module names from an aggregation config row.
#' @return A character vector of non-missing, non-empty module names.
#' @keywords internal

normalize_modules <- function(modules) {
  if (is.null(modules)) {
    return(character())
  }

  modules <- unlist(modules, use.names = FALSE)
  modules <- as.character(modules)
  modules[!is.na(modules) & modules != ""]
}

#' Test module membership for aggregations
#'
#' Identify which aggregation config rows opted into a module.
#'
#' @param modules A list-column or vector of configured module names.
#' @param module_name Module name to test for.
#' @return A logical vector with one value per `modules` element.
#' @keywords internal

aggregation_has_module <- function(modules, module_name) {
  purrr::map_lgl(as.list(modules), \(x) module_name %in% normalize_modules(x))
}

#' Validate aggregation module names
#'
#' Check that module names configured in `cfg_aggregations.yaml` are known to
#' the root target graph.
#'
#' @param aggregation_tibble Aggregation mapping tibble with a `modules` column.
#' @param known_modules Character vector of supported module names.
#' @param aggregation_config_file Path to the aggregation YAML config.
#' @return Invisibly returns `NULL`; errors on unknown configured modules.
#' @keywords internal

validate_aggregation_module_names <- function(
  aggregation_tibble,
  known_modules,
  aggregation_config_file = "cfg_aggregations.yaml"
) {
  configured_modules <- aggregation_tibble$modules |>
    as.list() |>
    purrr::map(normalize_modules) |>
    unlist(use.names = FALSE) |>
    unique()

  unknown_modules <- setdiff(configured_modules, known_modules)

  if (length(unknown_modules) > 0) {
    stop(
      aggregation_config_file,
      " references unknown module(s): ",
      paste(unknown_modules, collapse = ", "),
      ". Valid modules: ",
      paste(known_modules, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Read module config tibble
#'
#' Read a module YAML config and keep empty module configs compatible with the
#' downstream validation path.
#'
#' @param config_file Path to a module YAML config file.
#' @param module_name Module identifier used in error messages.
#' @param module_aggregation_tibble Aggregations that opted into the module,
#'   filtered from the main aggregation config.
#' @param aggregation_tibble Main aggregation config tibble defining all valid
#'   aggregation names.
#' @param aggregation_config_file Path to the aggregation YAML config.
#' @param manifest_file Path to the pipeline parameter manifest TSV.
#' @return A tibble with at least an `aggregation` column.
#' @keywords internal

read_module_config_tibble <- function(
  config_file,
  module_name,
  module_aggregation_tibble,
  aggregation_tibble,
  aggregation_config_file = "cfg_aggregations.yaml",
  manifest_file = "cfg_pipeline_parameters.tsv"
) {
  module_config_tibble <- read_manifest_config_tibble(
    config_file = config_file,
    manifest_file = manifest_file,
    scope = module_name,
    key_col = "aggregation"
  )

  if (!"aggregation" %in% names(module_config_tibble)) {
    module_config_tibble <- tibble::tibble(aggregation = character())
  }

  missing_module_rows <- setdiff(module_aggregation_tibble$aggregation, module_config_tibble$aggregation)
  if (length(missing_module_rows) > 0) {
    stop(
      aggregation_config_file,
      " opts aggregation(s) into module '",
      module_name,
      "' that are missing from ",
      config_file,
      ": ",
      paste(missing_module_rows, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_module_rows <- setdiff(module_config_tibble$aggregation, aggregation_tibble$aggregation)
  if (length(unknown_module_rows) > 0) {
    stop(
      config_file,
      " contains aggregation(s) not defined in ",
      aggregation_config_file,
      ": ",
      paste(unknown_module_rows, collapse = ", "),
      call. = FALSE
    )
  }

  module_config_tibble
}

#' Define a target symbol column
#'
#' Create a compact specification for adding one list-column of target symbols
#' from a source column of target suffixes.
#'
#' @param target Base target name without the mapped suffix.
#' @param from Name of the source column containing target suffix values.
#' @param sep Separator used between `target` and each suffix.
#' @param transform Function applied to each source value before target names
#'   are constructed.
#' @return A list consumed by `add_target_sym_cols()`.
#' @keywords internal

target_sym_col <- function(target, from, sep = ".", transform = identity) {
  list(target = target, from = from, sep = sep, transform = transform)
}

#' Add target symbol columns
#'
#' Add one or more list-columns of `rlang::syms()` values from compact target
#' symbol specifications.
#'
#' @param .data Mapping tibble to mutate.
#' @param ... Named target symbol specifications created by `target_sym_col()`.
#' @return `.data` with the requested target symbol list-columns added.
#' @keywords internal

add_target_sym_cols <- function(.data, ...) {
  target_sym_specs <- rlang::list2(...)
  target_sym_cols <- purrr::imap(target_sym_specs, \(spec, col_name) {
    purrr::map(
      .data[[spec$from]],
      \(x) {
        target_suffixes <- spec$transform(x)
        rlang::syms(stringr::str_c(spec$target, target_suffixes, sep = spec$sep))
      }
    )
  })

  .data |>
    dplyr::mutate(!!!target_sym_cols)
}

#' Add aggregation target symbols
#'
#' Add module-local target symbol columns where each symbol is suffixed by the
#' aggregation name from the module config row.
#'
#' @param module_tibble Module mapping tibble with an `aggregation` column.
#' @param target_names Character vector of base target names to suffix by
#'   aggregation.
#' @return `module_tibble` with one symbol column per target name.
#' @keywords internal

add_aggregation_target_syms <- function(module_tibble, target_names) {
  target_sym_cols <- target_names |>
    purrr::set_names() |>
    purrr::map(\(target_name) {
      purrr::map(
        module_tibble$aggregation,
        \(aggregation) rlang::sym(stringr::str_c(target_name, aggregation, sep = "."))
      )
    })

  module_tibble |>
    dplyr::mutate(!!!target_sym_cols)
}

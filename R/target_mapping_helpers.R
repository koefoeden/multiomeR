#' Validate processing and aggregation config
#'
#' Check that aggregation YAML rows reference valid, non-empty reaction sets.
#'
#' @param reaction_tibble Reaction config tibble containing the valid `reaction_ID`
#'   values.
#' @param aggregation_tibble_from_yaml Aggregation config tibble with
#'   `aggregation` labels and list-column `aggregation_reaction_IDs`.
#' @param reaction_config_file Path to the reaction TSV config.
#' @param aggregation_config_file Path to the aggregation YAML config.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

validate_processing_and_aggregation_config <- function(
  reaction_tibble,
  aggregation_tibble_from_yaml,
  reaction_config_file = "cfg_reactions.tsv",
  aggregation_config_file = "cfg_aggregations.yaml"
) {
  configured_aggregation_reaction_IDs <- unique(unlist(
    aggregation_tibble_from_yaml$aggregation_reaction_IDs
  ))

  unknown_aggregation_reactions_vec <- base::setdiff(
    configured_aggregation_reaction_IDs,
    reaction_tibble$reaction_ID
  )
  if (length(unknown_aggregation_reactions_vec) > 0) {
    stop(
      aggregation_config_file,
      " references reaction ID(s) not defined in ",
      reaction_config_file,
      ": ",
      paste(unknown_aggregation_reactions_vec, collapse = ", "),
      call. = FALSE
    )
  }

  empty_aggregation_vec <- aggregation_tibble_from_yaml$aggregation[
    lengths(aggregation_tibble_from_yaml$aggregation_reaction_IDs) == 0
  ]
  if (length(empty_aggregation_vec) > 0) {
    stop(
      aggregation_config_file,
      " contains aggregation(s) without reaction IDs: ",
      paste(empty_aggregation_vec, collapse = ", "),
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

assert_donor_reaction_metadata_column_ownership <- function(donor_id_metadata_tibble, reaction_ID_metadata_tibble) {
  overlapping_cols <- intersect(
    setdiff(colnames(donor_id_metadata_tibble), "donor_id"),
    setdiff(colnames(reaction_ID_metadata_tibble), "TENX_reaction_ID")
  )

  if (length(overlapping_cols) > 0) {
    stop(
      "Donor and reaction metadata contain overlapping non-key column(s): ",
      paste(overlapping_cols, collapse = ", "),
      ". Each metadata column must belong to exactly one table.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Build reaction mapping tibble
#'
#' Read reaction-level config and attach dataset-level YAML config values for
#' use by the reaction `tar_map()`.
#'
#' @param dataset_tibble_from_yaml Dataset config tibble created from
#'   `cfg_datasets.yaml`.
#' @param reaction_config_file Path to the reaction TSV config.
#' @param dataset_config_file Path to the dataset YAML config.
#' @return A tibble with one row per reaction and joined dataset config columns.
#' @keywords internal

build_reaction_tibble <- function(
  dataset_tibble_from_yaml,
  reaction_config_file = "cfg_reactions.tsv",
  dataset_config_file = "cfg_datasets.yaml"
) {
  reaction_tibble <- readr::read_tsv(reaction_config_file, show_col_types = FALSE) |>
    dplyr::mutate(dplyr::across(c(reaction_ID, reaction_donor_id), as.character))

  unknown_datasets <- setdiff(reaction_tibble$dataset, dataset_tibble_from_yaml$dataset)
  if (length(unknown_datasets) > 0) {
    stop(
      reaction_config_file,
      " references dataset(s) not defined in ",
      dataset_config_file,
      ": ",
      paste(unknown_datasets, collapse = ", "),
      call. = FALSE
    )
  }

  reaction_tibble |>
    dplyr::left_join(dataset_tibble_from_yaml, by = "dataset")
}

#' Build aggregation mapping tibble
#'
#' Filter active aggregations, validate their reaction references, and add the
#' list-columns used to splice upstream target symbols into aggregation targets.
#'
#' @param aggregation_tibble_all_from_yaml Aggregation config tibble created
#'   from `cfg_aggregations.yaml`.
#' @param reaction_tibble Reaction mapping tibble created by
#'   `build_reaction_tibble()`.
#' @param aggregation_config_file Path to the aggregation YAML config.
#' @param reaction_config_file Path to the reaction TSV config.
#' @return A tibble with one row per active aggregation, QC feature columns, and
#'   target symbol list-columns.
#' @keywords internal

build_aggregation_tibble <- function(
  aggregation_tibble_all_from_yaml,
  reaction_tibble,
  aggregation_config_file = "cfg_aggregations.yaml",
  reaction_config_file = "cfg_reactions.tsv"
) {
  aggregation_tibble_from_yaml <- aggregation_tibble_all_from_yaml |>
    dplyr::filter(purrr::map_lgl(is_active, isTRUE))

  aggregation_tibble <- aggregation_tibble_from_yaml |>
    add_target_sym_cols(
      aggregation_GEX_counts_BPCells_matrix_syms = target_sym_col("GEX_counts_BPCells_matrix", "aggregation_reaction_IDs"),
      aggregation_fragments_w_prefix_bpcells_syms = target_sym_col("fragments_w_prefix_bpcells", "aggregation_reaction_IDs"),
      aggregation_cellranger_summary_file_syms = target_sym_col("cellranger_summary_file", "aggregation_reaction_IDs"),
      aggregation_cellranger_kept_metadata_tibble_syms = target_sym_col("cellranger_kept_metadata_tibble", "aggregation_reaction_IDs"),
      aggregation_unfiltered_cells_n_vecs_syms = target_sym_col("unfiltered_cells_n_vecs", "aggregation_reaction_IDs"),
      aggregation_excluded_barcodes_by_type_list_syms = target_sym_col("excluded_barcodes_by_type_list", "aggregation_reaction_IDs"),
      aggregation_excluded_cellranger_only_barcodes_by_type_list_syms = target_sym_col("excluded_cellranger_only_barcodes_by_type_list", "aggregation_reaction_IDs"),
      aggregation_cellranger_ref_list_syms = target_sym_col("cellranger_ref_list", "aggregation_reaction_IDs", transform = \(ids) utils::head(ids, 1)),
      aggregation_gene_features_df_syms = target_sym_col("gene_features_df", "aggregation_reaction_IDs", transform = \(ids) utils::head(ids, 1))
    )

  project_categorical_vars <- purrr::map(
    aggregation_tibble$aggregation_categorical_vars,
    \(aggregation_categorical_vars) {
      c(
        aggregation_categorical_vars,
        "TENX_reaction_ID",
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
    reaction_tibble = reaction_tibble,
    aggregation_tibble_from_yaml = aggregation_tibble_from_yaml,
    reaction_config_file = reaction_config_file,
    aggregation_config_file = aggregation_config_file
  )

  aggregation_tibble |>
    dplyr::mutate(
      aggregation_dataset_vec = purrr::map(
        aggregation_reaction_IDs,
        \(ids) reaction_tibble$dataset[match(ids, reaction_tibble$reaction_ID)] |>
          unique()
      )
    ) |>
    add_target_sym_cols(
      per_dataset_excluded_upset_syms = target_sym_col("per_dataset_excluded_upset", "aggregation_dataset_vec"),
      per_dataset_excluded_cellranger_only_upset_syms = target_sym_col("per_dataset_excluded_cellranger_only_upset", "aggregation_dataset_vec"),
      per_dataset_QC_violins_syms = target_sym_col("per_dataset_QC_violins", "aggregation_dataset_vec")
    )
}

#' Build dataset mapping tibble
#'
#' Collapse reaction rows by dataset and attach dataset-level YAML config values
#' for use by the dataset `tar_map()`.
#'
#' @param reaction_tibble Reaction mapping tibble created by
#'   `build_reaction_tibble()`.
#' @param dataset_tibble_from_yaml Dataset config tibble created from
#'   `cfg_datasets.yaml`.
#' @return A tibble with one row per dataset and per-reaction target symbol
#'   list-columns.
#' @keywords internal

build_dataset_tibble <- function(reaction_tibble, dataset_tibble_from_yaml) {
  reaction_tibble |>
    dplyr::summarise(
      dataset_reaction_IDs = list(reaction_ID),
      .by = dataset
    ) |>
    add_target_sym_cols(
      dataset_unfiltered_cells_n_vecs_syms = target_sym_col("unfiltered_cells_n_vecs", "dataset_reaction_IDs"),
      dataset_cellranger_kept_metadata_tibble_syms = target_sym_col("cellranger_kept_metadata_tibble", "dataset_reaction_IDs"),
      dataset_excluded_cellranger_only_barcodes_by_type_list_syms = target_sym_col("excluded_cellranger_only_barcodes_by_type_list", "dataset_reaction_IDs"),
      dataset_excluded_barcodes_by_type_list_syms = target_sym_col("excluded_barcodes_by_type_list", "dataset_reaction_IDs")
    ) |>
    dplyr::left_join(dataset_tibble_from_yaml, by = "dataset")
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
#' @return A tibble with at least an `aggregation` column.
#' @keywords internal

read_module_config_tibble <- function(
  config_file,
  module_name,
  module_aggregation_tibble,
  aggregation_tibble,
  aggregation_config_file = "cfg_aggregations.yaml"
) {
  module_config_tibble <- read_config_tibble(
    config_file = config_file,
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

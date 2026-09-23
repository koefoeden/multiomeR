#' Validate processing and aggregation config
#'
#' Check that aggregation YAML rows reference defined, active GEM wells.
#'
#' @param GEM_well_tibble GEM well config tibble containing the valid `GEM_well_ID`
#'   values.
#' @param aggregation_tibble_from_yaml Aggregation config tibble with
#'   `aggregation` labels and list-column `aggregation_GEM_well_IDs`.
#' @return Invisibly returns after the validation or setup side effect succeeds.
#' @keywords internal

validate_processing_and_aggregation_config <- function(GEM_well_tibble, aggregation_tibble_from_yaml) {
  aggregation_config_file <- configuration_path("cfg_aggregations.yaml")
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
      configuration_path("cfg_GEM_wells.tsv"),
      ": ",
      paste(unknown_aggregation_GEM_wells_vec, collapse = ", "),
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

subset_keyed_metadata_tibble <- function(metadata_tibble, key_col, keys, source_label) {
  keys <- as.character(keys)
  missing_keys <- setdiff(keys, metadata_tibble[[key_col]])
  if (length(missing_keys) > 0L) {
    stop(
      source_label,
      " is missing configured ",
      key_col,
      " value(s): ",
      paste(missing_keys, collapse = ", "),
      call. = FALSE
    )
  }

  metadata_tibble[
    match(keys, metadata_tibble[[key_col]]),
    ,
    drop = FALSE
  ]
}

project_keyed_metadata_tibble <- function(
  metadata_tibble,
  key_col,
  requested_columns = NULL,
  strict = FALSE,
  canonical = FALSE
) {
  requested_columns <- requested_columns %||% character()
  requested_columns <- unique(as.character(unlist(requested_columns, use.names = FALSE)))
  requested_columns <- requested_columns[
    !is.na(requested_columns) & nzchar(requested_columns) & requested_columns != "NULL"
  ]

  if (strict) {
    missing_columns <- setdiff(
      requested_columns,
      c(key_col, colnames(metadata_tibble))
    )
    if (length(missing_columns) > 0L) {
      stop(
        "Metadata is missing requested column(s): ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (canonical) {
    requested_columns <- sort(setdiff(requested_columns, key_col))
  }

  projected_metadata <- metadata_tibble |>
    dplyr::select(
      dplyr::all_of(key_col),
      dplyr::any_of(requested_columns)
    )

  if (canonical) {
    projected_metadata <- projected_metadata |>
      dplyr::arrange(.data[[key_col]])
  }

  projected_metadata
}

metadata_columns_in_filter_expressions <- function(filter_expressions) {
  filter_expressions <- filter_expressions %||% character()
  filter_expressions |>
    purrr::map(rlang::parse_expr) |>
    purrr::map(base::all.vars) |>
    unlist(use.names = FALSE) |>
    unique()
}

get_GEM_well_annotation_metadata_tibble <- function(GEM_well_metadata_tibble) {
  processing_columns <- c(
    "GEM_well_donor_id",
    "GEM_well_n_donors",
    "GEM_well_cellranger_arc_count_dir",
    "GEM_well_add_cellbender",
    "GEM_well_cellbender_h5_file",
    "GEM_well_donors_VCF_file",
    "GEM_well_QC_exclude_list",
    "GEM_well_is_active"
  )

  GEM_well_metadata_tibble |>
    dplyr::select(-dplyr::any_of(processing_columns))
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

build_GEM_well_tibble <- function(GEM_well_config_file = configuration_path("cfg_GEM_wells.tsv")) {
  required_columns <- c(
    "GEM_well_ID",
    "GEM_well_dataset",
    "GEM_well_donor_id",
    "GEM_well_n_donors",
    "GEM_well_cellranger_arc_count_dir",
    "GEM_well_add_cellbender",
    "GEM_well_cellbender_h5_file",
    "GEM_well_donors_VCF_file",
    "GEM_well_QC_exclude_list",
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
  if (anyNA(GEM_well_tibble$GEM_well_add_cellbender) ||
      anyNA(GEM_well_tibble$GEM_well_is_active)) {
    stop("Every GEM well must define GEM_well_add_cellbender and GEM_well_is_active.", call. = FALSE)
  }
  inconsistent_cellbender <- xor(
    GEM_well_tibble$GEM_well_add_cellbender,
    !is.na(GEM_well_tibble$GEM_well_cellbender_h5_file)
  )
  if (any(inconsistent_cellbender)) {
    stop(
      "GEM_well_add_cellbender and GEM_well_cellbender_h5_file disagree for GEM well(s): ",
      paste(GEM_well_tibble$GEM_well_ID[inconsistent_cellbender], collapse = ", "),
      call. = FALSE
    )
  }

  if (length(validation_aggregations())) {
    config <- read_aggregation_config_tibble()
    selected <- unlist(config$aggregation_GEM_well_IDs[vapply(config$is_active, isTRUE, logical(1))])
    missing <- setdiff(selected, GEM_well_tibble$GEM_well_ID)
    if (length(missing)) stop("Missing validation GEM wells: ", paste(missing, collapse = ", "))
    GEM_well_tibble$GEM_well_is_active <- GEM_well_tibble$GEM_well_ID %in% selected
  }

  GEM_well_tibble |>
    dplyr::mutate(
      dataset = .data$GEM_well_dataset,
      GEM_well_QC_exclude_list = purrr::map2(
        .data$GEM_well_QC_exclude_list,
        .data$GEM_well_ID,
        parse_GEM_well_QC_exclude_list
      ),
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
        "GEM_well_QC_exclude_list",
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
#' @return A tibble with one row per active aggregation, QC feature columns, and
#'   target symbol list-columns.
#' @keywords internal

build_aggregation_tibble <- function(aggregation_tibble_all_from_yaml, GEM_well_tibble) {
  aggregation_tibble_from_yaml <- aggregation_tibble_all_from_yaml |>
    dplyr::filter(purrr::map_lgl(is_active, isTRUE))

  aggregation_tibble <- aggregation_tibble_from_yaml |>
    add_GEM_well_target_syms(c(
      "GEX_counts_BPCells_matrix",
      "vireo_donor_ids_tibble",
      "cellranger_barcodes_tsv",
      "fragments_w_prefix_bpcells",
      "cellranger_summary_file",
      "cellranger_kept_metadata_tibble",
      "unfiltered_cells_n_vecs",
      "excluded_barcodes_by_type_list",
      "excluded_cellranger_only_barcodes_by_type_list",
      "cellranger_ref_list",
      "gene_features_df"
    ))

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
      aggregation_w_peaks_continuous_vars = peak_vars,
      aggregation_w_WNN_continuous_vars = purrr::map(
        peak_vars,
        \(feature_names) c(feature_names, PROCESSING_CONTINUOUS_QC_FEATURES$WNN)
      )
    )
  )

  validate_processing_and_aggregation_config(
    GEM_well_tibble = GEM_well_tibble,
    aggregation_tibble_from_yaml = aggregation_tibble_from_yaml
  )

  aggregation_tibble |>
    dplyr::mutate(
      # Scalar component counts keep PCA and LSI commands unchanged when only
      # the downstream dimension selection changes.
      aggregation_GEX_PCA_n_components = purrr::map_dbl(aggregation_GEX_data_PCs, max),
      aggregation_ATAC_LSI_n_components = purrr::map_dbl(aggregation_ATAC_data_PCs, max),
      aggregation_GEM_well_QC_exclude_list = purrr::map(
        aggregation_GEM_well_IDs,
        \(ids) {
          stats::setNames(
            GEM_well_tibble$GEM_well_QC_exclude_list[
              match(ids, GEM_well_tibble$GEM_well_ID)
            ],
            ids
          )
        }
      )
    )
}

#' Build dataset config tibble
#'
#' Collapse GEM well-level settings by dataset for interactive configuration
#' loading.
#'
#' @param GEM_well_tibble GEM well mapping tibble created by
#'   `build_GEM_well_tibble()`.
#' @return A tibble with one row per dataset and collapsed GEM well settings.
#' @keywords internal

build_dataset_config_tibble <- function(GEM_well_tibble) {
  GEM_well_tibble |>
    dplyr::summarise(
      dataset_QC_exclude_list = list(unique(unlist(
        .data$GEM_well_QC_exclude_list,
        use.names = FALSE
      ))),
      is_active = list(any(unlist(.data$is_active))),
      .by = dataset
    )
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

#' Build module mapping values
#'
#' Keep the active aggregations that opt into a module and join the settings
#' from the module's configuration file, which must cover exactly those
#' aggregations that are defined in the aggregation config.
#'
#' @param module_name Module identifier, as listed under `modules`.
#' @param aggregation_tibble Active aggregation mapping tibble.
#' @param aggregation_tibble_all_from_yaml All configured aggregations.
#' @return The module's aggregation rows with the module settings joined.
#' @keywords internal

build_module_tibble <- function(module_name, aggregation_tibble, aggregation_tibble_all_from_yaml) {
  module_aggregation_tibble <- aggregation_tibble |>
    dplyr::filter(purrr::map_lgl(modules, \(modules) module_name %in% modules))
  if (nrow(module_aggregation_tibble) == 0L) {
    return(module_aggregation_tibble)
  }

  config_file <- configuration_path(paste0("cfg_module_", module_name, ".yaml"), must_exist = FALSE)
  if (!file.exists(config_file)) {
    stop("Missing configuration for enabled module '", module_name, "': ", config_file, call. = FALSE)
  }
  module_config_tibble <- read_manifest_config_tibble(
    config_file = config_file,
    manifest_file = "cfg_pipeline_parameters.tsv",
    scope = module_name,
    key_col = "aggregation"
  )
  missing_module_rows <- setdiff(module_aggregation_tibble$aggregation, module_config_tibble[["aggregation"]])
  if (length(missing_module_rows) > 0) {
    stop(
      "cfg_aggregations.yaml opts aggregation(s) into module '", module_name, "' that are missing from ",
      config_file, ": ", paste(missing_module_rows, collapse = ", "),
      call. = FALSE
    )
  }
  unknown_module_rows <- setdiff(module_config_tibble$aggregation, aggregation_tibble_all_from_yaml$aggregation)
  if (length(unknown_module_rows) > 0) {
    stop(
      config_file, " contains aggregation(s) not defined in cfg_aggregations.yaml: ",
      paste(unknown_module_rows, collapse = ", "),
      call. = FALSE
    )
  }

  dplyr::left_join(module_aggregation_tibble, module_config_tibble, by = "aggregation")
}

#' Add GEM well target symbols
#'
#' Add one list-column per base target name, named
#' `aggregation_<target>_syms`, holding that target's symbols for each of the
#' aggregation's GEM wells.
#'
#' @param aggregation_tibble Aggregation mapping tibble with the list-column
#'   `aggregation_GEM_well_IDs`.
#' @param target_names Character vector of base target names to suffix by
#'   GEM well ID.
#' @return `aggregation_tibble` with one symbol list-column per target name.
#' @keywords internal

add_GEM_well_target_syms <- function(aggregation_tibble, target_names) {
  target_sym_cols <- purrr::map(
    purrr::set_names(target_names, paste0("aggregation_", target_names, "_syms")),
    \(target_name) {
      purrr::map(
        aggregation_tibble$aggregation_GEM_well_IDs,
        \(GEM_well_IDs) rlang::syms(stringr::str_c(target_name, GEM_well_IDs, sep = "."))
      )
    }
  )
  dplyr::mutate(aggregation_tibble, !!!target_sym_cols)
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

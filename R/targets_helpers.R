tar_name_wo_suffixes <- function(target_name = targets::tar_name()) {
  target_name %>%
    stringr::str_remove("_[[:alnum:]]{16}$") %>%
    stringr::str_split_1("\\.") %>%
    utils::head(n = 1)
}

#' Get structured output path
#'
#' Build a deterministic target-derived output path under the targets store for plots or file artifacts.
#'
#' @param kind Top-level output category, usually `plots`, `files`, or another structured-output directory.
#' @param filetype Output file extension without a leading dot; currently constrained by the helper to supported graphics formats.
#' @param override_suffix Optional replacement for the target-derived filename suffix; use `NULL` to keep the default target name.
#' @param full_target_name Targets name, including branch suffix when present, used to derive output paths.
#' @param suffix_in_subdir Logical; when TRUE, place target suffixes in subdirectories instead of appending them to filenames.
#' @param list_output Logical; when TRUE, return/create a directory for multiple files rather than a single file path.
#' @return A length-one filesystem path, or a directory path when `list_output = TRUE`.
#' @keywords internal

get_structured_output_path <- function(
  kind = c("files", "plots"),
  filetype = NULL,
  override_suffix = NULL,
  full_target_name = targets::tar_name(),
  suffix_in_subdir = FALSE,
  list_output = FALSE
) {
  kind <- match.arg(kind)
  dyn_suffix <- full_target_name %>%
    stringr::str_extract("_[[:alnum:]]{16}$") %>%
    stringr::str_remove("^_")
  target_name <- full_target_name %>%
    stringr::str_remove("_[[:alnum:]]{16}$")
  target_parts <- target_name %>%
    stringr::str_split_1("\\.")
  stem <- target_parts[[1]]
  parent_dirs <- rev(target_parts[-1])
  parent_dir <- do.call(file.path, as.list(c(targets::tar_config_get("store"), kind, parent_dirs)))
  suffix <- dplyr::coalesce(override_suffix, dyn_suffix)
  has_suffix <- !is.na(suffix) && nzchar(suffix)

  if (suffix_in_subdir) {
    if (!has_suffix) {
      stop("A non-empty suffix is required when `suffix_in_subdir = TRUE`.", call. = FALSE)
    }
    out_dir <- file.path(parent_dir, stem, suffix)
    if (list_output || is.null(filetype)) {
      return(out_dir)
    }
    return(file.path(dirname(out_dir), paste0(basename(out_dir), ".", filetype)))
  }

  leaf <- stem
  if (has_suffix) {
    leaf <- paste0(leaf, "_", suffix)
  }
  out_path <- file.path(parent_dir, leaf)
  if (list_output || is.null(filetype)) {
    return(out_path)
  }
  paste0(out_path, ".", filetype)
}
#' Return the preferred graphics device for a file type.
#'
#' @keywords internal
get_plot_device <- function(filetype) {
  switch(
    filetype,
    "svg" = svglite::svglite,
    "png" = ragg::agg_png
  )
}

#' Open a graphics device for non-ggplot plot saving.
#'
#' @keywords internal
open_plot_graphics_device <- function(path, filetype, plot_save_args) {
  device_args <- plot_save_args
  width <- device_args$width
  height <- device_args$height
  dpi <- dplyr::coalesce(device_args$dpi, device_args$res, 300)
  units <- dplyr::coalesce(device_args$units, "in")
  device_args$width <- NULL
  device_args$height <- NULL
  device_args$dpi <- NULL
  device_args$res <- NULL
  device_args$units <- NULL

  fs::dir_create(dirname(path))

  device_fn <- get_plot_device(filetype)

  switch(
    filetype,
    "svg" = do.call(
      device_fn,
      c(
        list(
          filename = path,
          width = width,
          height = height
        ),
        device_args
      )
    ),
    "png" = do.call(
      device_fn,
      c(
        list(
          filename = path,
          width = width,
          height = height,
          units = units,
          res = dpi
        ),
        device_args
      )
    )
  )
}

#' Draw a plot-like object on the active graphics device.
#'
#' @keywords internal
draw_plot_object <- function(plot) {
  if (inherits(plot, c("grob", "gTree", "gtable"))) {
    grid::grid.newpage()
    grid::grid.draw(plot)
    return(invisible(NULL))
  }
  if (inherits(plot, "recordedplot")) {
    grDevices::replayPlot(plot)
    return(invisible(NULL))
  }
  print(plot)
  invisible(NULL)
}

get_discrete_axis_break_count <- function(panel_params, axis) {
  axis_params <- panel_params[[axis]]
  if (is.null(axis_params) || !isTRUE(axis_params$scale_is_discrete)) {
    return(0L)
  }

  breaks <- axis_params$breaks
  if (is.function(breaks)) {
    breaks <- axis_params$get_breaks()
  }
  if (is.null(breaks) || inherits(breaks, "waiver")) {
    return(0L)
  }

  length(unique(as.character(breaks[!is.na(breaks)])))
}

get_ggplot_auto_save_dimensions <- function(plot) {
  plot_build <- ggplot2::ggplot_build(plot)
  plot_layout <- plot_build$layout$layout
  n_facet_cols <- length(unique(plot_layout$COL))
  n_facet_rows <- length(unique(plot_layout$ROW))

  x_breaks <- plot_build$layout$panel_params |>
    purrr::map_int(get_discrete_axis_break_count, axis = "x") |>
    max(0L)
  y_breaks <- plot_build$layout$panel_params |>
    purrr::map_int(get_discrete_axis_break_count, axis = "y") |>
    max(0L)

  list(
    width = max(10, 3 * n_facet_cols, 4 + 0.35 * x_breaks * n_facet_cols),
    height = max(10, 3 * n_facet_rows, 4 + 0.25 * y_breaks * n_facet_rows)
  )
}

add_ggplot_title_if_missing <- function(plot, title) {
  if (!inherits(plot, "ggplot") || is.null(title) || !nzchar(title)) {
    return(plot)
  }
  if (!is.null(plot$labels$title)) {
    return(plot)
  }

  plot + ggplot2::labs(title = title)
}

#' Cull dense discrete legends
#'
#' Remove discrete ggplot legends whose number of breaks exceeds a configured limit.
#'
#' @param plot ggplot, patchwork, recordedplot, or compatible plot object to inspect or save.
#' @param n_distinct_max Maximum number of discrete scale breaks to keep; scales
#'   with more breaks are hidden to avoid unreadable legends in high-cardinality
#'   metadata plots.
#' @return The input ggplot with selected discrete guides removed when they exceed the break limit.
#' @keywords internal

cull_dense_discrete_legends <- function(plot, n_distinct_max) {
  if (!inherits(plot, "ggplot") || is.infinite(n_distinct_max)) {
    return(plot)
  }

  legend_aesthetics <- c("colour", "color", "fill", "shape", "linetype", "size", "alpha")
  scales <- ggplot2::ggplot_build(plot)$plot$scales$scales |>
    as.list()
  aesthetics_to_cull <- scales |>
    purrr::keep(\(scale) {
      if (
        length(intersect(scale$aesthetics, legend_aesthetics)) == 0 ||
          !isTRUE(scale$is_discrete()) ||
          identical(scale$guide, "none")
      ) {
        return(FALSE)
      }

      breaks <- scale$get_breaks()
      length(unique(as.character(breaks[!is.na(breaks)]))) > n_distinct_max
    }) |>
    purrr::map(\(scale) intersect(scale$aesthetics, legend_aesthetics)) |>
    unlist(use.names = FALSE) |>
    unique()

  if (length(aesthetics_to_cull) == 0) {
    return(plot)
  }

  guide_args <- rep(list("none"), length(aesthetics_to_cull))
  names(guide_args) <- aesthetics_to_cull
  plot + do.call(ggplot2::guides, guide_args)
}

#' Save plots structured
#'
#' Save one plot or a named list of plots to structured target-derived paths, adding titles and suppressing dense legends where needed.
#'
#' @param plots Single plot object or named list of plot objects to save; list names become filename suffixes and optional plot titles.
#' @param filetype Output file extension without a leading dot; currently constrained by the helper to supported graphics formats.
#' @param override_suffix Optional replacement for the target-derived filename suffix; use `NULL` to keep the default target name.
#' @param dyn_suffix_in_subdir Logical; when TRUE, dynamic-branch suffixes are used as subdirectories instead of filename suffixes.
#' @param target_name Full targets name used to derive structured output paths; defaults to the currently running target.
#' @param discrete_legend_n_distinct_max Maximum number of discrete legend entries to keep before replacing dense legends with `guide = 'none'`.
#' @param ... Additional arguments forwarded to the graphics device or
#'   `ggplot2::ggsave()`, depending on the plot object.
#' @return A file path for a single plot, or an output directory path for a list of plots.
#' @keywords internal

save_plots_structured <- skip_invalidate(function(
  plots,
  filetype = "png",
  override_suffix = NULL,
  dyn_suffix_in_subdir = FALSE,
  target_name = targets::tar_name(),
  discrete_legend_n_distinct_max = 20,
  ...
) {
  save_args <- list(...)
  filetype <- match.arg(filetype, c("svg", "png"))
  if (
    !is.numeric(discrete_legend_n_distinct_max) ||
      length(discrete_legend_n_distinct_max) != 1 ||
      is.na(discrete_legend_n_distinct_max) ||
      discrete_legend_n_distinct_max < 0
  ) {
    stop("`discrete_legend_n_distinct_max` must be a non-negative number or Inf.")
  }
  if (
    !is.null(override_suffix) &&
      (!is.character(override_suffix) || length(override_suffix) != 1)
  ) {
    stop("`override_suffix` must be NULL or a length-1 character vector.")
  }
  if (
    !is.logical(dyn_suffix_in_subdir) ||
      length(dyn_suffix_in_subdir) != 1 ||
      is.na(dyn_suffix_in_subdir)
  ) {
    stop("`dyn_suffix_in_subdir` must be TRUE or FALSE.")
  }
  is_plain_list <- is.list(plots) && identical(class(plots), "list")
  is_empty_plot_list <- inherits(plots, "empty_plot_list")
  if (is_empty_plot_list) {
    return(character())
  }
  is_single_plot <- !is.null(plots) && !is_plain_list
  is_plot_list <- is_plain_list && length(plots) > 0 && all(!purrr::map_lgl(plots, is.null))
  if (!is_single_plot && !is_plot_list) {
    stop("`plots` must be a plot object or a non-empty list of plot objects.")
  }
  n_plots <- if (is_plot_list) length(plots) else 1L
  for (dimension_arg in c("width", "height")) {
    if (
      !is.null(save_args[[dimension_arg]]) &&
        length(save_args[[dimension_arg]]) != 1 &&
        (!is_plot_list || length(save_args[[dimension_arg]]) != n_plots)
    ) {
      if (is_plot_list) {
        stop(stringr::str_glue("`{dimension_arg}` must be length 1 or length {n_plots} when saving a plot list."))
      }
      stop(stringr::str_glue("`{dimension_arg}` must be length 1 when saving a single plot."))
    }
  }
  dyn_suffix <- stringr::str_match(target_name, "_([[:alnum:]]{16})$")[, 2]
  if (dyn_suffix_in_subdir && is.na(dyn_suffix)) {
    stop("`dyn_suffix_in_subdir = TRUE` requires a dynamically suffixed target name.")
  }
  save_one_plot <- function(plot, path, plot_index = 1L) {
    plot_save_args <- save_args
    for (dimension_arg in c("width", "height")) {
      if (!is.null(plot_save_args[[dimension_arg]]) && length(plot_save_args[[dimension_arg]]) > 1) {
        plot_save_args[[dimension_arg]] <- plot_save_args[[dimension_arg]][[plot_index]]
      }
    }
    if (inherits(plot, "ggplot")) {
      plot <- cull_dense_discrete_legends(plot, discrete_legend_n_distinct_max)
      auto_dimensions <- get_ggplot_auto_save_dimensions(plot)
      plot_save_args$width <- dplyr::coalesce(plot_save_args$width, auto_dimensions$width)
      plot_save_args$height <- dplyr::coalesce(plot_save_args$height, auto_dimensions$height)
      if (identical(filetype, "png")) {
        plot_save_args$res <- dplyr::coalesce(plot_save_args$res, plot_save_args$dpi, 300)
      }
      do.call(
        ggplot2::ggsave,
        c(
          list(
            filename = path,
            plot = plot,
            device = get_plot_device(filetype)
          ),
          limitsize = FALSE,
          create.dir = TRUE,
          plot_save_args
        )
      )
      return(path)
    } else {
      plot_save_args$width <- dplyr::coalesce(plot_save_args$width, 10)
      plot_save_args$height <- dplyr::coalesce(plot_save_args$height, 10)
    }
    open_plot_graphics_device(path, filetype, plot_save_args)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_plot_object(plot)
    path
  }
  if (is_single_plot) {
    out_path <- get_structured_output_path(
      kind = "plots",
      filetype = filetype,
      override_suffix = override_suffix,
      full_target_name = target_name,
      suffix_in_subdir = dyn_suffix_in_subdir
    )
    return(save_one_plot(plots, out_path))
  }
  out_dir <- get_structured_output_path(
    kind = "plots",
    override_suffix = override_suffix,
    full_target_name = target_name,
    suffix_in_subdir = dyn_suffix_in_subdir,
    list_output = TRUE
  )
  plot_names <- names(plots)
  if (is.null(plot_names)) {
    plot_names <- rep("", length(plots))
  }
  index_width <- max(2, nchar(length(plots)))
  file_indices <- sprintf(paste0("%0", index_width, "d"), seq_along(plots))
  sanitized_names <- stringr::str_replace_all(plot_names, "[/\\\\]", "_")
  file_stems <- ifelse(
    nzchar(sanitized_names),
    paste0(file_indices, "_", sanitized_names),
    file_indices
  )
  purrr::pmap_chr(
    list(plots, file_stems, plot_names, seq_along(plots)),
    \(plot, file_stem, plot_name, plot_index) {
      plot <- add_ggplot_title_if_missing(plot, plot_name)
      save_one_plot(plot, file.path(out_dir, paste0(file_stem, ".", filetype)), plot_index)
    }
  )
})


get_structured_file_path <- function(filetype = NULL, override_suffix = NULL) {
  file_path <- get_structured_output_path(
    kind = "files",
    filetype = filetype,
    override_suffix = override_suffix
  )

  if (is.null(filetype)) {
    fs::dir_create(file_path)
  } else {
    fs::dir_create(dirname(file_path))
  }
  return(file_path)
}

#' Get dataset specific skip patterns
#'
#' Build regex skip patterns that prevent dataset-specific targets outside the requested dataset set from running.
#'
#' @param dataset_pattern Regular expression selecting datasets whose targets should remain runnable during skip-pattern construction.
#' @return A single regex string combining configured `cfg_tar_make_skip_regex_patterns`
#'   and cross-dataset exclusions, or an empty character value when no configured
#'   skip patterns apply.
#' @keywords internal

get_dataset_specific_skip_patterns <- function(dataset_pattern = ".") {
  # intelligently find targets that should be skipped based on the config for the active datasets

  config_file <- "cfg_datasets.yaml"
  full_cfg_tibble <- get_cfg_per_dataset_tibble(config_file = config_file, verbose = FALSE)
  relevant_cfg_tibble <- full_cfg_tibble %>% dplyr::filter(stringr::str_detect(dataset, dataset_pattern))

  per_dataset_skip_pattern_vec <- relevant_cfg_tibble %>%
    dplyr::mutate(
      combined_tar_make_skip_regex_patterns = purrr::map_chr(
        cfg_tar_make_skip_regex_patterns,
        ~ if (is.null(.x)) "" else stringr::str_c("(", .x, ")") %>% stringr::str_flatten("|")
      ),
      dataset_pattern = stringr::str_c(".*", dataset, "$"),
      combined_skip_pattern = stringr::str_c(combined_tar_make_skip_regex_patterns, dataset_pattern)
    ) %>%
    dplyr::filter(combined_tar_make_skip_regex_patterns != "") %>%
    dplyr::pull(combined_skip_pattern)

  if (length(per_dataset_skip_pattern_vec) > 0) {
    cat(stringr::str_glue("Will skip targets matching the following patterns based on cfg_tar_make_skip_regex_patterns in {config_file}:\n\n"))
    purrr::walk(per_dataset_skip_pattern_vec, ~ cat("- ", .x, "\n"))
    cat("\n")
    combined_cfg_skip_patterns <- per_dataset_skip_pattern_vec %>% stringr::str_flatten("|")
  } else {
    combined_cfg_skip_patterns <- character(0)
  }

  # Create skip-pattterns to handle dataset-specific targets, i.e. always skip dataset-specific targets in all other datasets, e.g.:

  relevant_datasets <- relevant_cfg_tibble %>% dplyr::pull(dataset)
  all_datasets <- full_cfg_tibble %>% dplyr::pull(dataset)

  dataset_specific_skip_patterns <- all_datasets %>%
    purrr::map(\(dataset_name) {
      other_datasets <- purrr::discard(relevant_datasets, ~ . == dataset_name)
      if (length(other_datasets) == 0) {
        return(NULL)
      }
      other_datasets_pattern <- other_datasets %>% stringr::str_flatten("|")
      skip_pattern <- stringr::str_glue("^{dataset_name}_.*({other_datasets_pattern})$")
      return(skip_pattern)
    }) %>%
    unlist()

  cat(stringr::str_glue("Will skip data-specific targets in the non-relevant datasets using the following regex patterns:\n\n"))
  purrr::walk(dataset_specific_skip_patterns, ~ cat("- ", .x, "\n"))
  cat("\n")

  # combine all skip patterns
  final_skip_pattern <- combined_cfg_skip_patterns %>% c(dataset_specific_skip_patterns) %>% stringr::str_flatten("|")

  return(final_skip_pattern)
}

# Helpers to retain provenance when dynamic mapping over targets
# Should error if the thing is a tibble, as this will be unpacked.
# or alternatively handle so it is added to the tibble as a column.
#' Set name attr
#'
#' Attach a stable `name` attribute to vectors or list elements so dynamic branches keep provenance labels.
#'
#' @param x Object to annotate with a `name` attribute; lists are processed element-wise.
#' @param name Attribute value to attach; defaults to the names of `x` for lists.
#' @return The input object with a `name` attribute attached; lists are returned with attributes set element-wise.
#' @keywords internal

set_name_attr <- function(x, name = NULL) {
  if (is.null(name)) {
    name <- names(x)
  }

  # vectorize
  if (is.list(x)) {
    purrr::map2(
      .x = x,
      .y = name,
      ~ {
        attr(.x, "name") <- .y
        .x
      }
    )
  } else {
    attr(x, "name") <- name
    x
  }
}

# This function should fail if it returns NULL, i.e. no name is set or the custom attribute is not set either.
get_name_attr <- function(x) {
  if (identical(class(x), "list")) {
    purrr::map_chr(x, ~ attr(.x, "name"))
  } else {
    attr(x, "name")
  }
}

skip_w_dummy_file_if <- function(condition_expr, targets_list, skip_expr = fs::file_create("NULL.txt")) {
  # capture (quote) both pieces without evaluating now
  condition_lang <- substitute(condition_expr)
  skip_lang <- substitute(skip_expr)

  # build the hook: if (<condition>) <skip_expr> else .x
  hook_lang <- bquote(if (.(condition_lang)) .(skip_lang) else .x)

  # hand the *quoted* hook to the raw variant
  tarchetypes::tar_hook_outer_raw(targets_list, hook = hook_lang)
}


#' Tar read inner filename
#'
#' Resolve file paths inside a file target, optionally expanding over selected dynamic branches.
#'
#' @param name File target name whose stored inner file paths should be read.
#' @param branches Dynamic branch selection passed through to targets metadata/path resolution.
#' @param meta Targets metadata table used to locate branch paths; defaults are usually read from the active store.
#' @param path_store Targets store path containing metadata and file-target objects.
#' @return Character vector of file paths stored inside the requested file target or branches.
#' @keywords internal

tar_read_inner_filename <- function(name, branches, meta, path_store) {
  run_time <- targets:::tar_runtime
  old_store <- run_time$store
  run_time$store <- path_store
  on.exit(run_time$store <- old_store)
  index <- meta$name == name
  if (!any(index)) {
    targets::tar_throw_validate("target ", name, " not found")
  }
  row <- meta[max(which(index)), , drop = FALSE] # nolint
  record <- targets:::record_from_row(row = row, path_store = path_store)
  targets:::if_any(
    record$type %in% c("stem", "branch"),
    targets:::read_builder(record),
    record$children
  )
}


# Called by tar_load()
w_def <- function(packages) {
  targets::tar_option_get("packages") %>%
    c(packages) %>%
    unique()
}

#' Save a ggplot as both PDF (cairo_pdf) and PNG at explicit journal dimensions.
#'
#' Returns a character vector of both output paths so that a target with
#' format = "file" tracks both files for change detection.
#'
#' @param plot      A ggplot object.
#' @param filename  Base filename including extension (e.g. "fig1.pdf").
#'                  The extension is replaced automatically for each format.
#' @param width_mm  Figure width in mm. Use PUB_COL1_MM or PUB_COL2_MM.
#' @param height_mm Figure height in mm. Use values up to PUB_MAX_H_MM.
#' @param figures_dir Output directory. Defaults to <store>/plots/pub_figures/.
#' @param ...       Additional arguments forwarded to ggplot2::ggsave().

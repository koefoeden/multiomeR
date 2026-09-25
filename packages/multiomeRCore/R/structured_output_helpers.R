#' Append relevant configuration values to a plot caption
#'
#' @param plot A ggplot, patchwork composite, or list of plots. Empty results are preserved.
#' @param config_file Configuration filename displayed in the caption.
#' @param ... Named resolved parameter values; names identify configuration keys.
#' @param .max_value_chars Maximum characters per displayed value, including the
#'   trailing `...` when truncated. Applied after collapsing vectors and whitespace.
#' @return The plot with its existing caption followed by the parameter summary.
#' @keywords internal
add_plot_parameters <- function(plot, config_file, ..., .max_value_chars = 40L) {
  parameters <- list(...)
  if (inherits(plot, "empty_plot_list")) return(plot)
  if (is.list(plot) && !inherits(plot, c("ggplot", "patchwork"))) {
    return(lapply(plot, \(item) do.call(add_plot_parameters,
      c(list(plot = item, config_file = config_file, .max_value_chars = .max_value_chars), parameters))))
  }
  values <- vapply(parameters, \(value) {
    if (length(value) == 0L) "(none)" else if (is.list(value))
      paste(deparse(value), collapse = " ") else paste(value, collapse = ", ")
  }, character(1)) |>
    stringr::str_squish() |>
    stringr::str_trunc(width = .max_value_chars, ellipsis = "...")
  parameter_caption <- stringr::str_wrap(paste0(
    "Relevant parameters (", config_file, "): ",
    paste(names(parameters), values, sep = ": ", collapse = "; ")
  ), width = 110)
  is_composite <- inherits(plot, "patchwork")
  existing_caption <- if (is_composite) plot$patches$annotation$caption else plot$labels$caption
  caption <- paste(c(existing_caption, parameter_caption), collapse = "\n\n")
  plot <- plot + if (is_composite) patchwork::plot_annotation(caption = caption) else ggplot2::labs(caption = caption)
  align_plot_captions(plot)
}

align_plot_captions <- function(plot) {
  # ggtext textboxes set their own alignment and cannot be merged with element_text().
  caption_theme <- function(theme) ggplot2::theme(
    plot.caption = if (inherits(theme$plot.caption, c("element_blank", "element_textbox")))
      theme$plot.caption else ggplot2::element_text(hjust = 0, vjust = 1),
    plot.caption.position = "plot"
  )
  if (inherits(plot, "patchwork")) {
    for (i in seq_len(length(plot))) plot[[i]] <- align_plot_captions(plot[[i]])
    return(plot + patchwork::plot_annotation(theme = caption_theme(plot$patches$annotation$theme)))
  }
  if (inherits(plot, "ggplot")) return(plot + caption_theme(plot$theme))
  plot
}

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

get_ggplot_auto_save_dimensions <- function(plot_build) {
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
#' @param plot_build Result of `ggplot2::ggplot_build(plot)`.
#' @return The input ggplot with selected discrete guides removed when they exceed the break limit.
#' @keywords internal

cull_dense_discrete_legends <- function(plot, n_distinct_max, plot_build) {
  if (is.infinite(n_distinct_max)) {
    return(plot)
  }

  legend_aesthetics <- c("colour", "color", "fill", "shape", "linetype", "size", "alpha")
  scales <- plot_build$plot$scales$scales |>
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
#' @param ... Additional arguments forwarded to `ggplot2::ggsave()` and the
#'   graphics device. ggplot objects without explicit dimensions are sized from
#'   their facets and discrete breaks; other grid objects default to 10 x 10 in.
#' @details Rendering completes in staging before published files change.
#' A per-target inventory removes only previously owned files, including when
#' the new result is empty. On first use, existing ownership is recovered from
#' file paths registered in the targets store when available.
#' @return Image paths for the saved plots, or `character()` for an empty result.
#' @keywords internal

save_plots_structured <- function(
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
  is_plain_list <- is.list(plots) && identical(class(plots), "list")
  is_empty_plot_list <- inherits(plots, "empty_plot_list") || (is_plain_list && length(plots) == 0L)
  is_single_plot <- !is.null(plots) && !is_plain_list && !is_empty_plot_list
  is_plot_list <- is_plain_list && length(plots) > 0 && all(!purrr::map_lgl(plots, is.null))
  if (!is_single_plot && !is_plot_list && !is_empty_plot_list) {
    stop("`plots` must be a plot object or a list of plot objects (possibly empty).")
  }
  n_plots <- if (is_single_plot) 1L else length(plots)
  for (dimension_arg in c("width", "height")) {
    if (
      !is.null(save_args[[dimension_arg]]) &&
        length(save_args[[dimension_arg]]) != 1 &&
        length(save_args[[dimension_arg]]) != n_plots
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
  save_one_plot <- function(plot, image_path, plot_index = 1L) {
    plot <- align_plot_captions(plot)
    plot_save_args <- save_args
    for (dimension_arg in c("width", "height")) {
      if (!is.null(plot_save_args[[dimension_arg]]) && length(plot_save_args[[dimension_arg]]) > 1) {
        plot_save_args[[dimension_arg]] <- plot_save_args[[dimension_arg]][[plot_index]]
      }
    }
    auto_dimensions <- list(width = 10, height = 10)
    if (inherits(plot, "ggplot")) {
      plot_build <- ggplot2::ggplot_build(plot)
      plot <- cull_dense_discrete_legends(plot, discrete_legend_n_distinct_max, plot_build)
      auto_dimensions <- get_ggplot_auto_save_dimensions(plot_build)
    }
    plot_save_args$width <- dplyr::coalesce(plot_save_args$width, auto_dimensions$width)
    plot_save_args$height <- dplyr::coalesce(plot_save_args$height, auto_dimensions$height)
    if (identical(filetype, "png")) {
      plot_save_args$res <- dplyr::coalesce(plot_save_args$res, plot_save_args$dpi, 300)
    }
    do.call(
      ggplot2::ggsave,
      c(
        list(
          filename = image_path,
          plot = plot,
          device = get_plot_device(filetype)
        ),
        limitsize = FALSE,
        create.dir = TRUE,
        plot_save_args
      )
    )
    image_path
  }
  # Build the complete destination inventory before touching published files.
  if (is_single_plot) {
    image_paths <- get_structured_output_path(
      kind = "plots", filetype = filetype, override_suffix = override_suffix,
      full_target_name = target_name, suffix_in_subdir = dyn_suffix_in_subdir
    )
    plots <- list(plots)
  } else {
    out_dir <- get_structured_output_path(
      kind = "plots", override_suffix = override_suffix,
      full_target_name = target_name, suffix_in_subdir = dyn_suffix_in_subdir,
      list_output = TRUE
    )
    plot_names <- names(plots)
    if (is.null(plot_names)) plot_names <- rep("", length(plots))
    file_indices <- sprintf(paste0("%0", max(2, nchar(length(plots))), "d"), seq_along(plots))
    sanitized_names <- stringr::str_replace_all(plot_names, "[/\\\\]", "_")
    file_stems <- ifelse(nzchar(sanitized_names), paste0(file_indices, "_", sanitized_names), file_indices)
    image_paths <- if (length(plots)) file.path(out_dir, paste0(file_stems, ".", filetype)) else character()
    plots <- purrr::map2(plots, plot_names, add_ggplot_title_if_missing)
  }
  store <- targets::tar_config_get("store")
  inventory_dir <- file.path(store, "plot_inventory")
  fs::dir_create(inventory_dir)
  inventory_file <- file.path(inventory_dir, paste0(target_name, ".rds"))
  if (file.exists(inventory_file)) {
    previous_paths <- file.path(store, readRDS(inventory_file))
  } else {
    # Older saves had no inventory. Never infer ownership from directory contents.
    registered <- tryCatch(
      suppressMessages(suppressWarnings(targets::tar_read_raw(target_name, store = store))),
      error = function(e) character()
    )
    previous_images <- if (is.character(registered)) registered else character()
    relative_paths <- fs::path_rel(previous_images, start = file.path(store, "plots"))
    previous_images <- previous_images[!grepl("^\\.\\.(/|$)", relative_paths) & grepl("\\.(png|svg)$", previous_images)]
    previous_paths <- previous_images
  }

  staging_dir <- tempfile(".staging-", tmpdir = inventory_dir)
  fs::dir_create(staging_dir)
  on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)
  staged_images <- file.path(staging_dir, paste0(seq_along(plots), ".", filetype))
  for (i in seq_along(plots)) {
    save_one_plot(plots[[i]], staged_images[[i]], i)
  }
  # Nothing above this point replaces or deletes the last successful output.
  for (i in seq_along(image_paths)) {
    fs::dir_create(dirname(image_paths[[i]]))
    if (!file.rename(staged_images[[i]], image_paths[[i]])) {
      stop("Could not publish plot output: ", image_paths[[i]], call. = FALSE)
    }
  }
  obsolete <- setdiff(fs::path_abs(previous_paths), fs::path_abs(image_paths))
  if (length(obsolete)) fs::file_delete(obsolete[file.exists(obsolete)])
  staged_inventory <- file.path(staging_dir, "inventory.rds")
  saveRDS(as.character(fs::path_rel(image_paths, start = store)), staged_inventory)
  if (!file.rename(staged_inventory, inventory_file)) {
    stop("Could not publish plot inventory: ", inventory_file, call. = FALSE)
  }
  image_paths
}


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

#' Derive serialized plot object path
#'
#' Convert a target-tracked plot image path into its mirrored serialized plot
#' object path.
#'
#' @param image_path Path to an image saved under a structured `plots`
#'   directory.
#' @return Character vector with the corresponding `.rds` path under
#'   `plot_objects`.
#' @keywords internal
shared_plot_object_path <- function(image_path) {
  image_path |>
    sub("/plots/", "/plot_objects/", x = _, fixed = TRUE) |>
    sub("[.](png|svg)$", ".rds", x = _)
}

#' Add outer panel tag
#'
#' Add a stable panel letter to a plot without relying on nested patchwork
#' automatic tagging.
#'
#' @param plot Plot-like object coercible with `ggplotify::as.ggplot()`.
#' @param tag Panel tag.
#' @return A ggplot object with the requested panel tag.
#' @keywords internal
shared_tag_panel <- function(plot, tag) {
  ggplotify::as.ggplot(plot) +
    ggplot2::labs(tag = tag) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
}

#' Pick plot output path
#'
#' Select one target output path from a dynamic/list file target by matching on
#' the file basename.
#'
#' @param image_paths Character vector or nested list of image paths.
#' @param pick_regex Optional regular expression required to match the basename.
#' @param exclude_regex Optional regular expression excluded from the basename.
#' @return One image path.
#' @keywords internal
shared_pick_plot_path <- function(image_paths, pick_regex = NULL, exclude_regex = NULL) {
  image_paths <- sort(unlist(image_paths, use.names = FALSE))

  if (!is.null(pick_regex) && nzchar(pick_regex)) {
    image_paths <- image_paths[grepl(pick_regex, basename(image_paths))]
  }
  if (!is.null(exclude_regex) && nzchar(exclude_regex)) {
    image_paths <- image_paths[!grepl(exclude_regex, basename(image_paths))]
  }
  if (!length(image_paths)) {
    stop("No image paths matched.", call. = FALSE)
  }
  image_paths[[1]]
}

#' Read target plot object
#'
#' Load a serialized plot object from the structured plot-object sidecar path
#' matching a target-tracked image path.
#'
#' @param target_call Expression that returns one or more target-tracked image
#'   paths.
#' @param pick_regex Optional basename regex passed to
#'   `shared_pick_plot_path()`.
#' @param exclude_regex Optional basename exclusion regex passed to
#'   `shared_pick_plot_path()`.
#' @return A plot object.
#' @keywords internal
shared_read_target_plot_object <- function(target_call, pick_regex = NULL, exclude_regex = NULL) {
  image_paths <- eval.parent(substitute(target_call))
  image_path <- shared_pick_plot_path(
    image_paths,
    pick_regex = pick_regex,
    exclude_regex = exclude_regex
  )

  object_path <- shared_plot_object_path(image_path)
  if (!file.exists(object_path)) {
    stop("Missing serialized plot object: ", object_path, call. = FALSE)
  }

  readRDS(object_path)
}

#' Prepare manuscript figure panels
#'
#' Build tagged figure panels and captions from a compact panel specification
#' table.
#'
#' @param panel_specs Tibble with `plot_call`, `pick_regex`, `title`, `tag`,
#'   `caption`, and optionally `plot_modifier`.
#' @param figure_theme Theme added to each panel.
#' @param show_tags Logical; if `TRUE`, add panel tags and combine panel
#'   captions.
#' @param envir Evaluation environment for plot calls and modifiers.
#' @return List with `plots` and `caption`.
#' @keywords internal
shared_prepare_plot_panels <- function(
  panel_specs,
  figure_theme = ggplot2::theme(),
  show_tags = TRUE,
  envir = parent.frame()
) {
  fig_caption <- if (show_tags) {
    paste(
      paste0("<b>", panel_specs$tag, ")</b> ", panel_specs$caption),
      collapse = " "
    )
  } else {
    panel_specs$caption[[1]]
  }
  if (!"plot_modifier" %in% names(panel_specs)) {
    panel_specs$plot_modifier <- rep(list(NA), nrow(panel_specs))
  }

  read_panel_plot <- function(plot_call, pick_regex, title, tag, plot_modifier) {
    pick_regex <- if (is.na(pick_regex)) NULL else pick_regex
    title <- if (is.na(title)) NULL else title
    plot_or_paths <- eval(plot_call, envir = envir)

    if (inherits(plot_or_paths, "ggplot")) {
      plot <- plot_or_paths
    } else {
      image_path <- shared_pick_plot_path(
        plot_or_paths,
        pick_regex = pick_regex
      )
      object_path <- shared_plot_object_path(image_path)
      if (!file.exists(object_path)) {
        stop("Missing serialized plot object: ", object_path, call. = FALSE)
      }
      plot <- readRDS(object_path)
    }

    if (inherits(plot, "patchwork")) {
      plot <- plot & figure_theme
    } else {
      plot <- plot + figure_theme
    }

    if (!is.null(title)) {
      if (inherits(plot, "patchwork")) {
        plot <- plot + patchwork::plot_annotation(title = title)
      } else {
        plot <- plot + ggplot2::labs(title = title)
      }
    }

    missing_modifier <- is.null(plot_modifier) ||
      (length(plot_modifier) == 1L && is.atomic(plot_modifier) && is.na(plot_modifier))
    if (!missing_modifier) {
      plot <- eval(plot_modifier, envir = list2env(list(plot = plot), parent = envir))
    }

    if (show_tags) {
      shared_tag_panel(plot, tag)
    } else {
      ggplotify::as.ggplot(plot)
    }
  }

  panel_plots <- purrr::pmap(
    panel_specs[c("plot_call", "pick_regex", "title", "tag", "plot_modifier")],
    read_panel_plot
  ) |>
    purrr::set_names(panel_specs$tag)

  list(plots = panel_plots, caption = fig_caption)
}

#' Render manuscript figures
#'
#' Render one or more patchwork figures from panel and figure specification
#' tables.
#'
#' @param panel_specs Panel specification tibble.
#' @param figure_specs Figure specification tibble with `figure_number`,
#'   `height`, `layout`, and `subtitle`.
#' @param output_dir Output directory for PNG files.
#' @param figure_theme Theme added to each panel.
#' @param figure_label Figure caption label.
#' @param width Output figure width in inches.
#' @param dpi Output PNG resolution.
#' @param envir Evaluation environment for layouts and panel calls.
#' @return Character vector of written PNG paths.
#' @keywords internal
shared_render_figures <- function(
  panel_specs,
  figure_specs,
  output_dir,
  figure_theme = ggplot2::theme(),
  figure_label = "Figure",
  width = 7,
  dpi = 300,
  envir = parent.frame()
) {
  render_one_figure <- function(figure_number) {
    figure_spec <- figure_specs[
      figure_specs$figure_number == figure_number,
      ,
      drop = FALSE
    ]
    if (nrow(figure_spec) != 1L) {
      stop(
        "Expected one figure spec for figure_number: ",
        figure_number,
        call. = FALSE
      )
    }

    figure_panel_specs <- panel_specs[
      panel_specs$figure_number == figure_number,
      ,
      drop = FALSE
    ]
    figure_panel_specs$figure_number <- NULL

    prepared_panels <- shared_prepare_plot_panels(
      figure_panel_specs,
      figure_theme = figure_theme,
      show_tags = nrow(figure_panel_specs) > 1L,
      envir = envir
    )

    layout_env <- list2env(as.list(prepared_panels$plots), parent = envir)
    subtitle <- figure_spec$subtitle[[1]]
    caption_text <- c(subtitle, prepared_panels$caption)
    caption_text <- caption_text[!is.na(caption_text) & nzchar(caption_text)]
    fig_caption <- if (length(caption_text)) {
      paste0(
        "<b>",
        figure_label,
        " ",
        figure_number,
        ":</b> ",
        paste(caption_text, collapse = " ")
      )
    } else {
      NULL
    }

    combined_plot <- eval(figure_spec$layout[[1]], envir = layout_env) +
      patchwork::plot_annotation(
        caption = fig_caption,
        theme = ggplot2::theme(
          plot.caption = ggtext::element_textbox_simple(
            hjust = 0,
            halign = 0,
            lineheight = 1.05,
            margin = ggplot2::margin(t = 6),
            width = grid::unit(1, "npc")
          ),
          plot.caption.position = "plot"
        )
      )

    output_png <- file.path(output_dir, paste0(figure_number, ".png"))
    ggplot2::ggsave(
      filename = output_png,
      plot = combined_plot,
      width = width,
      height = figure_spec$height[[1]],
      dpi = dpi,
      limitsize = FALSE,
      create.dir = TRUE
    )

    output_png
  }

  purrr::map_chr(figure_specs$figure_number, render_one_figure)
}

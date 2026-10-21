#' Wrap a plot as one figure panel
#'
#' Patchworks are wrapped rather than converted to grobs, which would measure
#' text on whichever device is open instead of the one that saves the figure.
#'
#' @param plot Plot-like object coercible with `ggplotify::as.ggplot()`.
#' @return A ggplot object.
#' @keywords internal
as_figure_panel <- function(plot) {
  if (inherits(plot, "patchwork")) {
    patchwork::wrap_elements(plot = plot)
  } else if (inherits(plot, "ggplot")) {
    plot
  } else {
    ggplotify::as.ggplot(plot)
  }
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
tag_panel <- function(plot, tag) {
  as_figure_panel(plot) +
    ggplot2::labs(tag = tag) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 1.03),
      plot.margin = ggplot2::margin(t = 4)
    )
}

#' Prepare manuscript figure panels
#'
#' Build tagged figure panels and captions from a compact panel specification
#' table.
#'
#' @param panel_specs Tibble with `tag` and a `plot_input` list-column of plot
#'   objects. Optional columns are `title`, `caption`, and `plot_modifier`;
#'   modifiers must be functions accepting and returning one plot.
#' @param figure_theme Theme added to each panel. By default, plot subtitles and
#'   plot captions are omitted because explanatory detail belongs in the
#'   assembled figure caption.
#' @param show_tags Logical; if `TRUE`, add panel tags and combine panel
#'   captions.
#' @return List with `plots` and `caption`.
#' @keywords internal
prepare_plot_panels <- function(
  panel_specs,
  figure_theme = ggplot2::theme(
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  ),
  show_tags = TRUE
) {
  required_columns <- c("tag", "plot_input")
  missing_columns <- setdiff(required_columns, colnames(panel_specs))
  if (length(missing_columns)) {
    stop(
      "Missing required panel specification columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(panel_specs)) {
    stop("`panel_specs` must contain at least one panel.", call. = FALSE)
  }
  if (anyNA(panel_specs$tag) || any(!nzchar(panel_specs$tag))) {
    stop("Panel tags must be non-empty strings.", call. = FALSE)
  }
  if (anyDuplicated(panel_specs$tag)) {
    stop("Panel tags must be unique.", call. = FALSE)
  }

  for (column_name in c("title", "caption")) {
    if (!column_name %in% colnames(panel_specs)) {
      panel_specs[[column_name]] <- rep(NA_character_, nrow(panel_specs))
    }
  }
  if (!"plot_modifier" %in% colnames(panel_specs)) {
    panel_specs$plot_modifier <- rep(list(NULL), nrow(panel_specs))
  }

  caption_is_present <- !is.na(panel_specs$caption) & nzchar(panel_specs$caption)
  fig_caption <- if (!any(caption_is_present)) {
    NULL
  } else if (show_tags) {
    paste(
      paste0(
        "<b>",
        panel_specs$tag[caption_is_present],
        ")</b> ",
        panel_specs$caption[caption_is_present]
      ),
      collapse = " "
    )
  } else {
    panel_specs$caption[[which(caption_is_present)[[1]]]]
  }
  prepare_panel_plot <- function(plot_input, title, tag, plot_modifier) {
    title <- if (is.na(title)) NULL else title
    if (!inherits(plot_input, c("ggplot", "grob", "gTree", "gtable", "recordedplot"))) {
      stop("Panel `plot_input` entries must be plot objects.", call. = FALSE)
    }
    plot <- plot_input

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

    if (!is.null(plot_modifier)) {
      if (!is.function(plot_modifier)) {
        stop("Panel plot modifiers must be functions or NULL.", call. = FALSE)
      }
      plot <- plot_modifier(plot)
    }

    if (show_tags) tag_panel(plot, tag) else as_figure_panel(plot)
  }

  panel_plots <- purrr::pmap(
    panel_specs[c("plot_input", "title", "tag", "plot_modifier")],
    prepare_panel_plot
  ) |>
    purrr::set_names(panel_specs$tag)

  list(plots = panel_plots, caption = fig_caption)
}

#' Compose a manuscript figure
#'
#' Prepare explicitly supplied plot inputs and arrange them as a patchwork
#' figure without writing files.
#'
#' @param panel_specs Panel specification tibble accepted by
#'   `prepare_plot_panels()`.
#' @param layout Quoted patchwork layout expression referring to panels by tag.
#' @param figure_caption Optional figure-level caption text.
#' @param figure_theme Theme added to every panel. By default, plot subtitles
#'   and plot captions are omitted because explanatory detail belongs in the
#'   assembled figure caption.
#' @param show_tags Logical; if `TRUE`, add panel tags and tagged panel captions.
#' @param caption_label Optional caption prefix, such as `"Figure 2"`.
#' @param envir Parent environment used to evaluate the patchwork layout.
#' @return A composed plot or patchwork figure.
#' @keywords internal
compose_manuscript_figure <- function(
  panel_specs,
  layout,
  figure_caption = NULL,
  figure_theme = ggplot2::theme(
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  ),
  show_tags = nrow(panel_specs) > 1L,
  caption_label = NULL,
  envir = parent.frame()
) {
  loadNamespace("patchwork")
  prepared_panels <- prepare_plot_panels(
    panel_specs,
    figure_theme = figure_theme,
    show_tags = show_tags
  )

  layout_env <- list2env(as.list(prepared_panels$plots), parent = envir)
  caption_text <- c(figure_caption, prepared_panels$caption)
  caption_text <- caption_text[!is.na(caption_text) & nzchar(caption_text)]
  figure_caption <- if (length(caption_text)) {
    paste(caption_text, collapse = " ")
  } else {
    NULL
  }
  if (!is.null(figure_caption) && !is.null(caption_label) && nzchar(caption_label)) {
    figure_caption <- paste0("<b>", caption_label, ":</b> ", figure_caption)
  }

  combined_plot <- eval(layout, envir = layout_env)
  if (is.null(figure_caption)) {
    return(combined_plot)
  }

  combined_plot +
    patchwork::plot_annotation(
      caption = figure_caption,
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
}

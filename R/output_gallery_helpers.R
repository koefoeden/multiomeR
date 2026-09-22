#' Check whether a path is absolute
#'
#' @keywords internal
is_absolute_gallery_path <- function(path) {
  !is.na(path) & grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

#' Read the output gallery manifest
#'
#' @param manifest_file YAML manifest path.
#' @return A tibble with one row per gallery item.
#' @keywords internal
read_output_gallery_manifest <- function(manifest_file = file.path("website", "output_gallery.yaml")) {
  items <- yaml::read_yaml(manifest_file)
  if (length(items) == 0) {
    return(tibble::tibble())
  }

  purrr::map_dfr(items, \(item) {
    tibble::tibble(
      id = item$id %||% NA_character_,
      section = item$section %||% NA_character_,
      subsection = item$subsection %||% NA_character_,
      title = item$title %||% NA_character_,
      description = item$description %||% NA_character_,
      target = item$target %||% NA_character_,
      asset = item$asset %||% NA_character_,
      status = item$status %||% NA_character_,
      source_file = item$source_file %||% NA_character_
    )
  })
}

#' Check output gallery asset status
#'
#' @param manifest_file YAML manifest path.
#' @param gallery_root Root directory for relative gallery asset paths.
#' @return Manifest tibble with resolved asset paths and existence flags.
#' @keywords internal
check_output_gallery_assets <- function(manifest_file = file.path("website", "output_gallery.yaml"), gallery_root = "website") {
  manifest <- read_output_gallery_manifest(manifest_file)
  gallery_ids <- stats::na.omit(manifest$id)

  if (anyDuplicated(gallery_ids)) {
    stop("Output gallery IDs must be unique.", call. = FALSE)
  }
  if (any(!grepl("^[a-z][a-z0-9-]*$", gallery_ids))) {
    stop(
      "Output gallery IDs must start with a lowercase letter and contain only lowercase letters, numbers, and hyphens.",
      call. = FALSE
    )
  }

  manifest |>
    dplyr::mutate(
      asset_path = dplyr::if_else(
        is.na(.data$asset) | is_absolute_gallery_path(.data$asset),
        .data$asset,
        file.path(gallery_root, .data$asset)
      ),
      asset_exists = !is.na(.data$asset_path) & file.exists(.data$asset_path)
    )
}

#' Render one output gallery card
#'
#' @param item One-row gallery item represented as a list.
#' @keywords internal
render_gallery_item <- function(item) {
  id_attribute <- if (is.na(item$id)) {
    ""
  } else {
    sprintf(' id="%s"', htmltools::htmlEscape(item$id))
  }
  title <- htmltools::htmlEscape(item$title)
  description <- htmltools::htmlEscape(item$description)
  target <- htmltools::htmlEscape(item$target)
  asset <- item$asset
  asset_exists <- isTRUE(item$asset_exists)

  cat(sprintf('<article%s class="output-gallery-card">\n', id_attribute))
  if (asset_exists) {
    cat(sprintf('<a href="%s"><img src="%s" alt="%s"></a>\n', asset, asset, title))
  } else {
    cat('<div class="output-gallery-placeholder">Preview pending</div>\n')
  }
  cat(sprintf('<p class="output-gallery-title">%s</p>\n', title))
  cat(sprintf('<p>%s</p>\n', description))
  cat(sprintf('<p class="output-gallery-target"><code>%s</code></p>\n', target))
  cat('</article>\n')
}

#' Render output gallery cards in a responsive grid
#'
#' @param items Gallery manifest rows to render.
#' @keywords internal
render_gallery_grid <- function(items) {
  cat('<div class="output-gallery-grid">\n')
  purrr::pwalk(items, \(...) render_gallery_item(list(...)))
  cat('</div>\n')
}

#' Render selected output gallery cards
#'
#' @param gallery_items Gallery manifest rows from `check_output_gallery_assets()`.
#' @param ids Gallery item IDs to render, in display order.
#' @keywords internal
render_gallery_cards <- function(gallery_items, ids) {
  missing_ids <- setdiff(ids, gallery_items$id)
  if (length(missing_ids) > 0) {
    stop("Unknown output gallery IDs: ", paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  render_gallery_grid(gallery_items[match(ids, gallery_items$id), , drop = FALSE])
  invisible(NULL)
}

#' Render one output gallery section
#'
#' @param gallery_items Gallery manifest rows.
#' @param section Section name to render.
#' @param subsection_descriptions Optional named character vector of text to
#'   render below matching subsection headings.
#' @keywords internal
render_gallery_section <- function(gallery_items, section, subsection_descriptions = NULL) {
  section_items <- dplyr::filter(gallery_items, .data$section == .env$section)
  if (nrow(section_items) == 0) {
    stop("No output gallery items found for section: ", section, call. = FALSE)
  }

  if (all(is.na(section_items$subsection))) {
    render_gallery_grid(section_items)
    return(invisible(NULL))
  }

  for (subsection_name in unique(stats::na.omit(section_items$subsection))) {
    subsection_items <- dplyr::filter(section_items, .data$subsection == subsection_name)
    cat(sprintf("\n## %s\n\n", htmltools::htmlEscape(subsection_name)))
    if (!is.null(subsection_descriptions) && subsection_name %in% names(subsection_descriptions)) {
      cat(sprintf("%s\n\n", htmltools::htmlEscape(subsection_descriptions[[subsection_name]])))
    }
    render_gallery_grid(subsection_items)
  }

  invisible(NULL)
}

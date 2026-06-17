#' Default for NULL values
#'
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

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
  read_output_gallery_manifest(manifest_file) |>
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
  title <- htmltools::htmlEscape(item$title)
  description <- htmltools::htmlEscape(item$description)
  target <- htmltools::htmlEscape(item$target)
  asset <- item$asset
  asset_exists <- isTRUE(item$asset_exists)

  cat('<article class="output-gallery-card">\n')
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

#' Render one output gallery section
#'
#' @param gallery_items Gallery manifest rows.
#' @param section Section name to render.
#' @keywords internal
render_gallery_section <- function(gallery_items, section) {
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
    render_gallery_grid(subsection_items)
  }

  invisible(NULL)
}

#' Resolve source paths for output gallery targets
#'
#' @param manifest Manifest tibble from `read_output_gallery_manifest()`.
#' @param meta Targets metadata with at least `name`, `path`, and `error`.
#' @return Manifest tibble with one resolved source path per item.
#' @keywords internal
resolve_output_gallery_source_paths <- function(
  manifest = read_output_gallery_manifest(),
  meta = targets::tar_meta(fields = c("name", "path", "error"))
) {
  path_records <- meta |>
    dplyr::select("name", "path", "error") |>
    tidyr::unnest_longer("path", values_to = "target_path", keep_empty = TRUE)

  manifest |>
    dplyr::mutate(.gallery_row = dplyr::row_number()) |>
    dplyr::left_join(path_records, by = c("target" = "name")) |>
    dplyr::slice_head(n = 1, by = ".gallery_row") |>
    dplyr::select(-".gallery_row") |>
    dplyr::mutate(
      source_path = dplyr::if_else(
        !is.na(.data$source_file) & dir.exists(.data$target_path),
        file.path(.data$target_path, .data$source_file),
        .data$target_path
      ),
      source_exists = !is.na(.data$source_path) & file.exists(.data$source_path),
      target_error = !is.na(.data$error)
    )
}

#' Copy completed pipeline outputs into the docs gallery
#'
#' This is a maintainer helper. Run it after the public example targets have
#' completed, then commit the copied files under `website/gallery_assets/`.
#'
#' @param manifest_file YAML manifest path.
#' @param gallery_root Root directory for relative gallery asset paths.
#' @param overwrite Whether to overwrite existing gallery assets.
#' @param dry_run If `TRUE`, return the copy plan without copying files.
#' @return A tibble describing copied and skipped assets.
#' @keywords internal
collect_output_gallery_assets <- function(
  manifest_file = file.path("website", "output_gallery.yaml"),
  gallery_root = "website",
  overwrite = FALSE,
  dry_run = FALSE
) {
  copy_plan <- resolve_output_gallery_source_paths(
    manifest = read_output_gallery_manifest(manifest_file)
  ) |>
    dplyr::mutate(
      dest_path = dplyr::if_else(
        is.na(.data$asset) | is_absolute_gallery_path(.data$asset),
        .data$asset,
        file.path(gallery_root, .data$asset)
      ),
      dest_exists = !is.na(.data$dest_path) & file.exists(.data$dest_path),
      source_is_dir = !is.na(.data$source_path) & dir.exists(.data$source_path),
      dest_has_extension = !is.na(.data$dest_path) & nzchar(tools::file_ext(.data$dest_path)),
      can_copy = .data$source_exists &
        !is.na(.data$dest_path) &
        (overwrite | !.data$dest_exists) &
        (!.data$source_is_dir | !.data$dest_has_extension | !is.na(.data$source_file))
    )

  if (dry_run) {
    return(copy_plan)
  }

  purrr::pwalk(
    dplyr::select(copy_plan, "source_path", "dest_path", "can_copy"),
    \(source_path, dest_path, can_copy) {
      if (!isTRUE(can_copy)) {
        return(invisible(NULL))
      }
      dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(source_path)) {
        if (file.exists(dest_path) && overwrite) {
          unlink(dest_path, recursive = TRUE)
        }
        file.copy(source_path, dest_path, recursive = TRUE)
      } else {
        file.copy(source_path, dest_path, overwrite = overwrite)
      }
    }
  )

  check_output_gallery_assets(manifest_file = manifest_file, gallery_root = gallery_root) |>
    dplyr::left_join(
      dplyr::select(copy_plan, "target", "source_path", "source_exists", "dest_path", "can_copy"),
      by = "target"
    )
}

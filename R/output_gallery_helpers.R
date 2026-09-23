#' Render the output gallery from the manifest written by dev/refresh_output_gallery.py
#'
#' Each group becomes a level-2 heading and each checkpoint section a level-3
#' heading whose anchor is the checkpoint name, so guides can link to it.
#'
#' @param manifest_file Gallery manifest path.
#' @param gallery_root Directory that relative preview paths resolve against.
#' @keywords internal
render_output_gallery <- function(manifest_file = "output_gallery.yaml", gallery_root = ".") {
  items <- purrr::map_dfr(yaml::read_yaml(manifest_file)$items, \(item) tibble::tibble(
    checkpoint = item$checkpoint, group = item$group, section = item$section %||% NA_character_,
    target = item$target, description = item$description, asset = item$asset
  ))
  missing <- items$asset[!file.exists(file.path(gallery_root, items$asset))]
  if (length(missing)) stop("Missing gallery previews: ", paste(missing, collapse = ", "), call. = FALSE)
  anchor <- \(checkpoint) gsub("[^a-z0-9]+", "-", tolower(checkpoint))
  for (group in unique(items$group)) {
    group_items <- items[items$group == group, ]
    sectioned <- !anyNA(group_items$section)
    cat(sprintf("\n## %s%s\n\n", group, if (sectioned) "" else sprintf(" {#%s}", anchor(group_items$checkpoint[[1]]))))
    for (checkpoint in unique(group_items$checkpoint)) {
      section_items <- group_items[group_items$checkpoint == checkpoint, ]
      if (sectioned) cat(sprintf("\n### %s {#%s}\n\n", section_items$section[[1]], anchor(checkpoint)))
      render_gallery_grid(section_items)
    }
  }
  invisible(NULL)
}

#' Render gallery previews in a responsive grid
#'
#' @param items Gallery manifest rows with `target`, `description` and `asset`.
#' @keywords internal
render_gallery_grid <- function(items) {
  escape <- htmltools::htmlEscape
  cat('<div class="output-gallery-grid">\n')
  for (i in seq_len(nrow(items))) {
    cat(sprintf(paste0(
      '<figure class="output-gallery-card" id="%s">\n',
      '<a href="%s"><img src="%s" alt="%s" loading="lazy"></a>\n',
      '<figcaption>%s <code>%s</code></figcaption>\n</figure>\n'
    ), escape(items$target[[i]]), items$asset[[i]], items$asset[[i]], escape(items$description[[i]]),
    escape(items$description[[i]]), escape(items$target[[i]])))
  }
  cat('</div>\n')
}

#' Count assignments against the complete Cell Ranger-called barcode set.
summarize_demultiplexing_assignments <- function(vireo_results, barcode_files, GEM_well_IDs) {
  stopifnot(length(vireo_results) == length(GEM_well_IDs),
    length(barcode_files) == length(GEM_well_IDs), !anyDuplicated(GEM_well_IDs))
  purrr::pmap_dfr(list(vireo_results, barcode_files, GEM_well_IDs),
    function(result, barcode_file, GEM_well_ID) {
      barcodes <- readLines(barcode_file)
      if (anyDuplicated(barcodes) || anyDuplicated(result$barcode) ||
          !setequal(barcodes, result$barcode)) {
        stop("Vireo assignments must cover every Cell Ranger-called barcode exactly once: ", GEM_well_ID)
      }
      if (anyNA(result$vireo_type) ||
          any(!result$vireo_type %in% c("singlet", "doublet", "unassigned", "not_demultiplexed"))) {
        stop("Unexpected Vireo assignment class: ", GEM_well_ID)
      }
      counts <- result |>
        dplyr::mutate(
          assignment_class = dplyr::if_else(.data$vireo_type == "not_demultiplexed", "singlet", .data$vireo_type),
          donor_id = dplyr::if_else(.data$assignment_class == "singlet", as.character(.data$donor_id), NA_character_)
        )
      if (any(counts$assignment_class == "singlet" &
          (is.na(counts$donor_id) | !nzchar(counts$donor_id)))) {
        stop("Assigned singlets require a donor ID: ", GEM_well_ID)
      }
      counts |>
        dplyr::count(.data$assignment_class, .data$donor_id, name = "n_nuclei") |>
        dplyr::mutate(GEM_well_ID = GEM_well_ID, n_called = length(barcodes),
          fraction = .data$n_nuclei / .data$n_called,
          not_demultiplexed = all(result$vireo_type == "not_demultiplexed"))
    }) |>
    dplyr::mutate(
      GEM_well_ID = factor(.data$GEM_well_ID, levels = GEM_well_IDs),
      assignment_class = factor(.data$assignment_class, levels = c("unassigned", "doublet", "singlet"))
    ) |>
    dplyr::arrange(.data$GEM_well_ID, .data$assignment_class, .data$donor_id) |>
    dplyr::mutate(
      xmin = pmax(0, cumsum(.data$fraction) - .data$fraction),
      xmax = pmin(1, cumsum(.data$fraction)),
      .by = GEM_well_ID
    )
}

#' Plot donor-specific singlet segments with three assignment-class colors.
plot_demultiplexing_assignment_bars <- function(counts) {
  counts |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,
      ymin = as.numeric(GEM_well_ID) - 0.4, ymax = as.numeric(GEM_well_ID) + 0.4,
      fill = assignment_class), color = "white", linewidth = 0.2) +
    # The target saves at 16 inches wide; narrow donor segments need smaller, rotated labels.
    ggplot2::geom_text(data = dplyr::filter(counts, .data$assignment_class == "singlet"),
      ggplot2::aes(x = (xmin + xmax) / 2, y = as.numeric(GEM_well_ID), label = donor_id,
        angle = ifelse(fraction < 0.04, 90, 0), size = pmin(2.5, 300 * fraction))) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(counts$GEM_well_ID)),
      labels = levels(counts$GEM_well_ID)) +
    ggplot2::scale_fill_manual(values = c(unassigned = "#999999", doublet = "#D55E00", singlet = "#009E73"), drop = FALSE) +
    ggplot2::labs(title = "Donor assignments per GEM well before QC",
      subtitle = "Look for excess doublet or unassigned calls and uneven donor representation before pooling wells.",
      x = "Fraction of all Cell Ranger-called nuclei", y = NULL, fill = "Assignment",
      caption = "Singlet segments are labelled by donor ID. Single-donor wells use their configured donor and are shown as singlets; genotype doublets were not assessed in those wells.") +
    ggplot2::theme(legend.position = "top", panel.grid.major.y = ggplot2::element_blank())
}

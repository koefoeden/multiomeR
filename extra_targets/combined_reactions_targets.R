rlang::list2(
  targets::tar_target(
    name = combined_donor_id_per_barcode_tibble,
    dplyr::bind_rows(rlang::syms(paste0(cellranger_kept_metadata_tibble, ".", reactions_tibble$reaction_ID)))
  ),
  tarchetypes::tar_file(
    name = combined_feature_by_reaction_QC_violins,
    description = "Plot violin plots of QC metrics per reaction across all datasets and save to file",
    command = {
      plot_data <- combined_donor_id_per_barcode_tibble |>
        dplyr::select(dplyr::any_of(c("TENX_reaction_ID", "dataset", PROCESSING_PER_REACTION_QC_VARS))) |>
        tidyr::pivot_longer(cols = dplyr::any_of(PROCESSING_PER_REACTION_QC_VARS), names_to = "feature", values_to = "value") |>
        dplyr::filter(
          .by = c(TENX_reaction_ID, feature),
          value >= stats::quantile(value, probs = 0.02, na.rm = TRUE),
          value <= stats::quantile(value, probs = 0.98, na.rm = TRUE)
        )

      plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = TENX_reaction_ID, y = value, fill = dataset)) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::facet_wrap(~feature, scales = "free", ncol = 1) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.title.x = ggplot2::element_blank(), legend.position = "none")

      save_plots_structured(plot)
    }
  ),
  tarchetypes::tar_file(
    name = combined_amulet_barplot,
    description = "Plot proportion of AMULET-discarded barcodes per reaction as a stacked bar plot and save to file",
    command = {
      proportion_tibble <- combined_donor_id_per_barcode_tibble |>
        dplyr::filter(!vireo_type %in% c("discarded_by_cellranger", "not_demultiplexed")) |>
        dplyr::count(dataset, TENX_reaction_ID, discarded_by_amulet) |>
        dplyr::mutate(prop = n / sum(n), .by = c(dataset, TENX_reaction_ID))

      if (nrow(proportion_tibble) > 0) {
        plot <- ggplot2::ggplot(proportion_tibble, ggplot2::aes(x = TENX_reaction_ID, y = prop, fill = discarded_by_amulet)) +
          ggplot2::geom_col() +
          ggplot2::facet_wrap(~dataset, scales = "free_x") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        save_plots_structured(plot)
      } else {
        # no demultiplexing was done for any reactions
        "NULL.txt"
      }
    }
  ),
  tarchetypes::tar_file(
    name = combined_vireo_types_barplot,
    description = "Plot proportion of vireo assignment types per reaction as a stacked bar plot and save to file",
    command = {
      proportion_tibble <- combined_donor_id_per_barcode_tibble |>
        dplyr::filter(!vireo_type %in% c("discarded_by_cellranger", "not_demultiplexed")) |>
        dplyr::group_by(dataset, TENX_reaction_ID, vireo_type) |>
        dplyr::tally() |>
        dplyr::mutate(prop = n / sum(n))

      if (nrow(proportion_tibble) > 0) {
        plot <- ggplot2::ggplot(proportion_tibble, ggplot2::aes(x = TENX_reaction_ID, y = prop, fill = vireo_type)) +
          ggplot2::geom_col() +
          ggplot2::facet_wrap(~dataset, scales = "free_x") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        save_plots_structured(plot)
      } else {
        # no demultiplexing was done for any reactions
        "NULL.txt"
      }
    }
  )
)

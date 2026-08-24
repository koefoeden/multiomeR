rlang::list2(
  targets::tar_target(
    name = peak_gene_correlation_gene_TSS_tibble.peak_gene_correlation.WNN,
    description = "Build the gene TSS table used for ArchR-style peak-gene correlation candidates",
    command = make_peak_gene_correlation_gene_TSS_tibble(
      reference_Ensembl_annotations_GRanges_list = marker_validated_Ensembl_annotations_GRanges_list,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN,
    description = "Pair consensus ATAC peaks with same-chromosome gene TSSs within 250 kb",
    command = make_peak_gene_correlation_candidate_pairs(
      consensus_peak_GRanges = consensus_peak_GRanges.ATAC,
      gene_TSS_tibble = peak_gene_correlation_gene_TSS_tibble.peak_gene_correlation.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_chromosome_tibble.peak_gene_correlation.WNN,
    description = "List chromosomes with candidate peak-gene pairs",
    command = make_peak_gene_correlation_chromosome_tibble(
      peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = peak_gene_correlation_cell_group_diagnostics_tibble.peak_gene_correlation.WNN,
    description = "Record broad cell groups dropped before donor-state pseudobulking",
    command = make_peak_gene_correlation_cell_group_diagnostics(
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      ATAC_counts_matrix = consensus_peak_BPCells_matrix.ATAC,
      embedding_matrix = harmony_embeddings_matrix.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_cell_groups_tibble.peak_gene_correlation.WNN,
    description = "Split eligible broad WNN cell-type groups for donor-state pseudobulking",
    command = make_peak_gene_correlation_cell_groups(
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      ATAC_counts_matrix = consensus_peak_BPCells_matrix.ATAC,
      embedding_matrix = harmony_embeddings_matrix.ATAC
    ) |>
      dplyr::group_by(cell_group) |>
      targets::tar_group(),
    iteration = "group",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_donor_state_records.peak_gene_correlation.WNN,
    description = "Build mutually exclusive donor by ATAC-state pseudobulks within one broad cell type",
    command = make_peak_gene_correlation_donor_state_record(
      cell_group_tibble = peak_gene_correlation_cell_groups_tibble.peak_gene_correlation.WNN,
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      ATAC_embedding_matrix = harmony_embeddings_matrix.ATAC
    ),
    pattern = map(peak_gene_correlation_cell_groups_tibble.peak_gene_correlation.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_donor_state_aggregates_tibble.peak_gene_correlation.WNN,
    description = "Combine retained donor by ATAC-state pseudobulk memberships",
    command = combine_peak_gene_correlation_donor_state_records(
      records = peak_gene_correlation_donor_state_records.peak_gene_correlation.WNN,
      component = "aggregates"
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_donor_state_diagnostics_tibble.peak_gene_correlation.WNN,
    description = "Combine donor-state pseudobulk eligibility diagnostics",
    command = combine_peak_gene_correlation_donor_state_records(
      records = peak_gene_correlation_donor_state_records.peak_gene_correlation.WNN,
      component = "diagnostics"
    ),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = peak_gene_correlation_group_chromosome_tibble.peak_gene_correlation.WNN,
    description = "Split peak-gene correlation work by broad cell type and chromosome",
    command = make_peak_gene_correlation_group_chromosome_tibble(
      donor_state_aggregates_tibble = peak_gene_correlation_donor_state_aggregates_tibble.peak_gene_correlation.WNN,
      chromosome_tibble = peak_gene_correlation_chromosome_tibble.peak_gene_correlation.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN
    ),
    iteration = "group",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_aggregate_matrices.peak_gene_correlation.WNN,
    description = "Aggregate chromosome-specific GEX and ATAC counts over donor-state pseudobulks",
    command = make_peak_gene_correlation_aggregate_matrices(
      group_chromosome_tibble = peak_gene_correlation_group_chromosome_tibble.peak_gene_correlation.WNN,
      donor_state_aggregates_tibble = peak_gene_correlation_donor_state_aggregates_tibble.peak_gene_correlation.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      ATAC_counts_matrix = consensus_peak_BPCells_matrix.ATAC
    ),
    pattern = map(peak_gene_correlation_group_chromosome_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN,
    description = "Normalize donor-state pseudobulk GEX and ATAC matrices as log1p CPM",
    command = normalize_peak_gene_correlation_aggregate_matrices(
      peak_gene_correlation_aggregate_matrices.peak_gene_correlation.WNN
    ),
    pattern = map(peak_gene_correlation_aggregate_matrices.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_results_tibbles.peak_gene_correlation.WNN,
    description = "Score donor-adjusted peak-gene associations for one broad cell type and chromosome",
    command = score_peak_gene_correlations_for_cell_group(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN
    ),
    pattern = map(peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = peak_gene_correlation_branch_diagnostics_tibbles.peak_gene_correlation.WNN,
    description = "Diagnose skipped or retained peak-gene correlation branches",
    command = diagnose_peak_gene_correlation_branch(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN
    ),
    pattern = map(peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_results_tibble.peak_gene_correlation.WNN,
    description = "Combine chromosome-level peak-gene correlations and add within-cell-group FDR",
    command = finalize_peak_gene_correlation_results(
      results_tibble = peak_gene_correlation_results_tibbles.peak_gene_correlation.WNN,
      aggregation = aggregation
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_links_tibble.peak_gene_correlation.WNN,
    description = "Filter positive candidate enhancer-gene links at FDR < 0.05",
    command = make_peak_gene_correlation_links(peak_gene_correlation_results_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_diagnostics_tibble.peak_gene_correlation.WNN,
    description = "Combine peak-gene correlation cell-group and branch diagnostics",
    command = dplyr::bind_rows(
      peak_gene_correlation_cell_group_diagnostics_tibble.peak_gene_correlation.WNN,
      peak_gene_correlation_donor_state_diagnostics_tibble.peak_gene_correlation.WNN,
      peak_gene_correlation_branch_diagnostics_tibbles.peak_gene_correlation.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_correlation_histogram_plot.peak_gene_correlation.WNN,
    description = "Save a facetted histogram of peak-gene correlations by cell group",
    command = {
      correlation_bin_width <- 0.025
      correlation_plot_tibble <- summarize_peak_gene_correlation_histogram(
        results_tibble = peak_gene_correlation_results_tibble.peak_gene_correlation.WNN,
        bin_width = correlation_bin_width
      )
      plot_peak_gene_correlation_histogram(
        plot_tibble = correlation_plot_tibble,
        bin_width = correlation_bin_width
      ) |>
        save_plots_structured()
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_support_counts_plot.peak_gene_correlation.WNN,
    description = "Save peak-gene correlation tested, significant, and linked pair counts by cell group",
    command = {
      support_plot_tibble <- summarize_peak_gene_correlation_support_counts(
        peak_gene_correlation_results_tibble.peak_gene_correlation.WNN
      )
      plot_peak_gene_correlation_support_counts(support_plot_tibble) |>
        save_plots_structured()
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_distance_correlation_plot.peak_gene_correlation.WNN,
    description = "Save median peak-gene correlation by absolute TSS distance and cell group",
    command = {
      distance_plot_tibble <- summarize_peak_gene_correlation_by_distance(
        peak_gene_correlation_results_tibble.peak_gene_correlation.WNN
      )
      plot_peak_gene_correlation_by_distance(distance_plot_tibble) |>
        save_plots_structured()
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_diagnostics_plot.peak_gene_correlation.WNN,
    description = "Save peak-gene correlation retained and skipped branch diagnostics by cell group",
    command = {
      diagnostics_plot_tibble <- peak_gene_correlation_diagnostics_tibble.peak_gene_correlation.WNN |>
        dplyr::mutate(skipped_reason = dplyr::coalesce(.data$skipped_reason, "retained")) |>
        dplyr::count(.data$cell_group, .data$skipped_reason, name = "n_branches")

      diagnostics_plot <- ggplot2::ggplot(
        diagnostics_plot_tibble,
        ggplot2::aes(x = .data$n_branches, y = .data$cell_group, fill = .data$skipped_reason)
      ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          x = "Branches",
          y = "Cell group",
          fill = "Status"
        )

      save_plots_structured(diagnostics_plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_links_tibble.peak_gene_correlation.WNN,
    description = "Split top peak-gene links per cell group for downstream QC plots",
    command = {
      peak_gene_correlation_links_tibble.peak_gene_correlation.WNN |>
        dplyr::filter(.data$rank_in_cell_group <= 3) |>
        dplyr::arrange(.data$cell_group, .data$rank_in_cell_group) |>
        dplyr::mutate(
          scatter_plot_name = paste(
            make.names(.data$cell_group),
            sprintf("rank%03d", .data$rank_in_cell_group),
            make.names(.data$TargetGene),
            .data$chr,
            sep = "_"
          )
        ) |>
        dplyr::left_join(
          peak_gene_correlation_candidate_pairs_tibble.peak_gene_correlation.WNN |>
            dplyr::select("peak", "TargetGeneID", "gene_matrix_feature") |>
            dplyr::distinct(),
          by = c("peak", "TargetGeneID")
        )
    },
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_aggregate_values_tibbles.peak_gene_correlation.WNN,
    description = "Extract top-link aggregate values from one normalized cell-group/chromosome branch",
    command = extract_peak_gene_correlation_top_link_aggregate_values(
      normalized_aggregate_matrices =
        peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN,
      top_links_tibble =
        peak_gene_correlation_top_links_tibble.peak_gene_correlation.WNN
    ),
    pattern = map(
      peak_gene_correlation_normalized_aggregate_matrices.peak_gene_correlation.WNN
    ),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_aggregate_values_tibble.peak_gene_correlation.WNN,
    description = "Combine compact top-link aggregate values across cell groups and chromosomes",
    command = dplyr::bind_rows(
      peak_gene_correlation_top_link_aggregate_values_tibbles.peak_gene_correlation.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_aggregate_scatter_tibble.peak_gene_correlation.WNN,
    description = "Build aggregate-level GEX and ATAC values for one top peak-gene link",
    command = {
      link_row <- peak_gene_correlation_top_links_tibble.peak_gene_correlation.WNN
      peak_gene_correlation_top_link_aggregate_values_tibble.peak_gene_correlation.WNN |>
        dplyr::filter(
          .data$scatter_plot_name == link_row$scatter_plot_name[[1]]
        )
    },
    pattern = map(peak_gene_correlation_top_links_tibble.peak_gene_correlation.WNN),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 32) # apparently 16 GB is not enough
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_top_link_aggregate_scatter_plots.peak_gene_correlation.WNN,
    description = "Save aggregate-level scatterplots for top peak-gene links per cell group",
    command = {
      plot_tibble <- peak_gene_correlation_top_link_aggregate_scatter_tibble.peak_gene_correlation.WNN
      trend_label_tibble <- split(plot_tibble, plot_tibble$cell_group) |>
        purrr::map_dfr(\(cell_group_tibble) {
          cell_group_tibble <- cell_group_tibble |>
            dplyr::filter(
              is.finite(.data$peak_accessibility_logCPM),
              is.finite(.data$gene_expression_logCPM)
            )
          label_x <- max(cell_group_tibble$peak_accessibility_logCPM, na.rm = TRUE)
          label_y <- if (
            nrow(cell_group_tibble) >= 2L &&
              dplyr::n_distinct(cell_group_tibble$peak_accessibility_logCPM) >= 2L
          ) {
            fit <- stats::lm(gene_expression_logCPM ~ peak_accessibility_logCPM, data = cell_group_tibble)
            as.numeric(stats::predict(
              fit,
              newdata = tibble::tibble(peak_accessibility_logCPM = label_x)
            ))
          } else {
            cell_group_tibble$gene_expression_logCPM[which.max(cell_group_tibble$peak_accessibility_logCPM)]
          }

          tibble::tibble(
            cell_group = cell_group_tibble$cell_group[[1]],
            label_x = label_x,
            label_y = label_y
          )
        })
      x_upper <- max(plot_tibble$peak_accessibility_logCPM, trend_label_tibble$label_x, na.rm = TRUE)
      y_upper <- max(plot_tibble$gene_expression_logCPM, trend_label_tibble$label_y, na.rm = TRUE)

      scatter_plot <- ggplot2::ggplot(
        plot_tibble,
        ggplot2::aes(
          x = .data$peak_accessibility_logCPM,
          y = .data$gene_expression_logCPM,
          color = .data$cell_group
        )
      ) +
        ggplot2::geom_jitter(
          width = 0.03,
          height = 0.03,
          size = 0.7,
          alpha = 0.45
        ) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
        ggrepel::geom_text_repel(
          data = trend_label_tibble,
          ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$cell_group, color = .data$cell_group),
          inherit.aes = FALSE,
          size = 3,
          max.overlaps = Inf,
          show.legend = FALSE
        ) +
        ggplot2::coord_cartesian(
          xlim = c(0, x_upper * 1.05 + 0.03),
          ylim = c(0, y_upper * 1.05 + 0.03)
        ) +
        ggplot2::labs(
          title = paste(plot_tibble$TargetGene[[1]], plot_tibble$peak[[1]], plot_tibble$primary_cell_group[[1]], sep = " - "),
          subtitle = paste0(
            "primary cell group = ",
            plot_tibble$primary_cell_group[[1]],
            "; r = ",
            round(plot_tibble$correlation[[1]], 3),
            "; FDR = ",
            signif(plot_tibble$FDR[[1]], 3),
            "; rank = ",
            plot_tibble$rank_in_cell_group[[1]]
          ),
          x = "Aggregate ATAC log1p CPM",
          y = "Aggregate GEX log1p CPM",
          color = "Cell group"
        )

      genome_annotation_track <- peak_gene_correlation_top_link_genome_annotation_track.peak_gene_correlation.WNN +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank()
        )
      primary_ATAC_track <- peak_gene_correlation_top_link_primary_cell_ATAC_accessibility_track.peak_gene_correlation.WNN +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank()
        )

      track_panel <- patchwork::wrap_plots(
        genome_annotation_track,
        primary_ATAC_track,
        peak_gene_correlation_top_link_peak_gene_loop_track.peak_gene_correlation.WNN,
        ncol = 1
      )
      plot <- patchwork::wrap_plots(
        track_panel,
        scatter_plot,
        ncol = 1,
        heights = c(1, 2)
      )

      plot |>
        save_plots_structured(
          override_suffix = plot_tibble$scatter_plot_name[[1]],
          dyn_suffix_in_subdir = TRUE
        )
    },
    pattern = map(
      peak_gene_correlation_top_link_aggregate_scatter_tibble.peak_gene_correlation.WNN,
      peak_gene_correlation_top_link_genome_annotation_track.peak_gene_correlation.WNN,
      peak_gene_correlation_top_link_primary_cell_ATAC_accessibility_track.peak_gene_correlation.WNN,
      peak_gene_correlation_top_link_peak_gene_loop_track.peak_gene_correlation.WNN
    )
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN,
    description = "Add plotting windows to top peak-gene links for per-locus track plots",
    command = {
      peak_gene_correlation_top_links_tibble.peak_gene_correlation.WNN |>
        dplyr::mutate(
          locus_start = pmax(
            1L,
            min(.data$start[[1]], .data$TargetGeneTSS[[1]], na.rm = TRUE) - 50000L
          ),
          locus_end = max(.data$end[[1]], .data$TargetGeneTSS[[1]], na.rm = TRUE) + 50000L
        )
    },
    pattern = map(peak_gene_correlation_top_links_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_genome_annotation_track.peak_gene_correlation.WNN,
    description = "Build the genome annotation track for one top peak-gene correlation locus",
    command = {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      genes_GRanges <- IRanges::subsetByOverlaps(
        marker_validated_Ensembl_annotations_GRanges_list$genes,
        region
      )
      gene_data_frame <- GenomicRanges::as.data.frame(genes_GRanges)
      gene_labels <- if ("gene_name" %in% colnames(gene_data_frame)) {
        gene_data_frame$gene_name
      } else if (!is.null(names(genes_GRanges))) {
        names(genes_GRanges)
      } else {
        paste0("gene_", seq_len(nrow(gene_data_frame)))
      }
      genes_tibble <- gene_data_frame |>
        tibble::as_tibble() |>
        dplyr::transmute(
          chr = as.character(.data$seqnames),
          start = .data$start,
          end = .data$end,
          gene_name = as.character(gene_labels)
        )

      gene_track <- BPCells::trackplot_genome_annotation(
        loci = genes_tibble,
        region = region,
        label_by = "gene_name",
        track_label = "Genes"
      )

      gene_track
    },
    pattern = map(peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_primary_cell_ATAC_accessibility_track.peak_gene_correlation.WNN,
    description = "Build the primary-cell-group BPCells ATAC coverage track for one top peak-gene correlation locus",
    command = {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN
      primary_cell_group <- locus_tibble$cell_group[[1]]
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      metadata_tibble <- metadata_w_cell_types_tibble.WNN
      fragments <- BPCells::select_cells(
        combined_BPCells_fragment_obj.ATAC,
        metadata_tibble$barcode_w_prefix
      )
      fragment_cell_names <- BPCells::cellNames(fragments)
      metadata <- metadata_tibble |>
        dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
        dplyr::filter(
          .data$barcode_w_prefix %in% fragment_cell_names,
          .data[["PCA_harmony_SNN_cluster_cell_type"]] == primary_cell_group
        ) |>
        dplyr::arrange(match(.data$barcode_w_prefix, fragment_cell_names))
      fragments <- BPCells::select_cells(fragments, metadata$barcode_w_prefix)

      cell_read_counts <- if ("atac_fragments" %in% colnames(metadata)) {
        metadata$atac_fragments
      } else {
        metadata$nCount_ATAC
      }

      coverage_tibble <- BPCells::trackplot_coverage(
        fragments = fragments,
        region = region,
        groups = metadata[["PCA_harmony_SNN_cluster_cell_type"]],
        cell_read_counts = cell_read_counts,
        group_order = primary_cell_group,
        bins = 500,
        return_data = TRUE
      )
      make_BPCells_ATAC_coverage_track_from_tibble(
        coverage_tibble = coverage_tibble,
        region = region
      )
    },
    pattern = map(peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_ATAC_accessibility_track.peak_gene_correlation.WNN,
    description = "Build the per-cell-group BPCells ATAC coverage track for one top peak-gene correlation locus",
    command = {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      metadata_tibble <- metadata_w_cell_types_tibble.WNN
      fragments <- BPCells::select_cells(
        combined_BPCells_fragment_obj.ATAC,
        metadata_tibble$barcode_w_prefix
      )
      fragment_cell_names <- BPCells::cellNames(fragments)
      metadata <- metadata_tibble |>
        dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
        dplyr::filter(
          .data$barcode_w_prefix %in% fragment_cell_names,
          !is.na(.data[["PCA_harmony_SNN_cluster_cell_type"]])
        ) |>
        dplyr::arrange(match(.data$barcode_w_prefix, fragment_cell_names))
      fragments <- BPCells::select_cells(fragments, metadata$barcode_w_prefix)

      groups <- metadata[["PCA_harmony_SNN_cluster_cell_type"]]
      cell_read_counts <- if ("atac_fragments" %in% colnames(metadata)) {
        metadata$atac_fragments
      } else {
        metadata$nCount_ATAC
      }
      group_order <- gtools::mixedsort(unique(as.character(groups)))

      coverage_tibble <- BPCells::trackplot_coverage(
        fragments = fragments,
        region = region,
        groups = groups,
        cell_read_counts = cell_read_counts,
        group_order = group_order,
        bins = 500,
        return_data = TRUE
      )
      coverage_track <- make_BPCells_ATAC_coverage_track_from_tibble(
        coverage_tibble = coverage_tibble,
        region = region
      )
      coverage_track
    },
    pattern = map(peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_peak_gene_loop_track.peak_gene_correlation.WNN,
    description = "Build the BPCells peak-gene loop track for filtered links at one top-link locus",
    command = {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      loops_tibble <- peak_gene_correlation_links_tibble.peak_gene_correlation.WNN |>
        dplyr::filter(
          .data$chr == locus_tibble$chr[[1]],
          .data$TargetGeneID == locus_tibble$TargetGeneID[[1]],
          .data$start <= locus_tibble$locus_end[[1]],
          .data$end >= locus_tibble$locus_start[[1]],
          .data$TargetGeneTSS >= locus_tibble$locus_start[[1]],
          .data$TargetGeneTSS <= locus_tibble$locus_end[[1]]
        ) |>
        dplyr::mutate(
          peak_position = as.integer(round((.data$start + .data$end) / 2)),
          cell_group = factor(.data$cell_group),
          abs_correlation = abs(.data$correlation)
        ) |>
        dplyr::transmute(
          chr = .data$chr,
          start = pmin(.data$peak_position, .data$TargetGeneTSS),
          end = pmax(.data$peak_position, .data$TargetGeneTSS),
          cell_group = .data$cell_group,
          abs_correlation = .data$abs_correlation
        )

      loop_track_tibble <- loops_tibble |>
        dplyr::select("chr", "start", "end", "cell_group")
      loop_width_tibble <- loops_tibble |>
        dplyr::mutate(loop_id = dplyr::row_number()) |>
        dplyr::select("loop_id", "abs_correlation")

      loop_data_tibble <- BPCells::trackplot_loop(
        loops = loop_track_tibble,
        region = region,
        color_by = "cell_group",
        track_label = "Links",
        return_data = TRUE
      )
      loop_data_tibble <- loop_data_tibble |>
        dplyr::left_join(loop_width_tibble, by = "loop_id")

      loop_track <- make_BPCells_peak_gene_loop_track_from_tibble(
        loop_data_tibble = loop_data_tibble,
        region = region
      )
      loop_track
    },
    pattern = map(peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_top_link_ATAC_tracks_plots.peak_gene_correlation.WNN,
    description = "Save combined genome annotation, ATAC coverage, and peak-gene loop tracks at top-link loci",
    command = {
      plot <- BPCells::trackplot_combine(
        tracks = list(
          peak_gene_correlation_top_link_genome_annotation_track.peak_gene_correlation.WNN,
          peak_gene_correlation_top_link_ATAC_accessibility_track.peak_gene_correlation.WNN,
          peak_gene_correlation_top_link_peak_gene_loop_track.peak_gene_correlation.WNN
        ),
        title = paste(
          peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN$TargetGene[[1]],
          peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN$peak[[1]],
          peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN$cell_group[[1]],
          sep = " - "
        )
      )

      plot |>
        save_plots_structured(
          override_suffix = peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN$scatter_plot_name[[1]],
          dyn_suffix_in_subdir = TRUE
        )
    },
    pattern = map(
      peak_gene_correlation_top_link_loci_tibble.peak_gene_correlation.WNN,
      peak_gene_correlation_top_link_genome_annotation_track.peak_gene_correlation.WNN,
      peak_gene_correlation_top_link_ATAC_accessibility_track.peak_gene_correlation.WNN,
      peak_gene_correlation_top_link_peak_gene_loop_track.peak_gene_correlation.WNN
    )
  )
)

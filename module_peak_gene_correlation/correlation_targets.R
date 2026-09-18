rlang::list2(
  targets::tar_target(
    name = peak_gene_correlation_gene_TSS_tibble.WNN,
    description = "Build the gene TSS table used for ArchR-style peak-gene correlation candidates [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_gene_TSS_tibble(
      reference_Ensembl_annotations_GRanges_list = marker_validated_Ensembl_annotations_GRanges_list,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_candidate_pairs_tibble.WNN,
    description = "Pair consensus ATAC peaks with same-chromosome gene TSSs within 250 kb [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_candidate_pairs(
      consensus_peak_GRanges = consensus_peak_GRanges.ATAC,
      gene_TSS_tibble = peak_gene_correlation_gene_TSS_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_chromosome_tibble.WNN,
    description = "List chromosomes with candidate peak-gene pairs [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_chromosome_tibble(
      peak_gene_correlation_candidate_pairs_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = peak_gene_correlation_cell_group_diagnostics_tibble.WNN,
    description = "Record broad cell groups dropped before donor-state pseudobulking [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_cell_group_diagnostics(
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      ATAC_counts_matrix = consensus_peak_BPCells_matrix.ATAC,
      embedding_matrix = harmony_embeddings_matrix.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_cell_groups_tibble.WNN,
    description = "Split eligible broad WNN cell-type groups for donor-state pseudobulking [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
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
    name = peak_gene_correlation_donor_state_records.WNN,
    description = "Build mutually exclusive donor by ATAC-state pseudobulks within one broad cell type [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_donor_state_record(
      cell_group_tibble = peak_gene_correlation_cell_groups_tibble.WNN,
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      ATAC_embedding_matrix = harmony_embeddings_matrix.ATAC
    ),
    pattern = map(peak_gene_correlation_cell_groups_tibble.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_donor_state_aggregates_tibble.WNN,
    description = "Combine retained donor by ATAC-state pseudobulk memberships [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = combine_peak_gene_correlation_donor_state_records(
      records = peak_gene_correlation_donor_state_records.WNN,
      component = "aggregates"
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_donor_state_diagnostics_tibble.WNN,
    description = "Combine donor-state pseudobulk eligibility diagnostics [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = combine_peak_gene_correlation_donor_state_records(
      records = peak_gene_correlation_donor_state_records.WNN,
      component = "diagnostics"
    ),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = peak_gene_correlation_group_chromosome_tibble.WNN,
    description = "Split peak-gene correlation work by broad cell type and chromosome [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_group_chromosome_tibble(
      donor_state_aggregates_tibble = peak_gene_correlation_donor_state_aggregates_tibble.WNN,
      chromosome_tibble = peak_gene_correlation_chromosome_tibble.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN
    ),
    iteration = "group",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_aggregate_matrices.WNN,
    description = "Aggregate chromosome-specific GEX and ATAC counts over donor-state pseudobulks [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_aggregate_matrices(
      group_chromosome_tibble = peak_gene_correlation_group_chromosome_tibble.WNN,
      donor_state_aggregates_tibble = peak_gene_correlation_donor_state_aggregates_tibble.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN,
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      ATAC_counts_matrix = consensus_peak_BPCells_matrix.ATAC
    ),
    pattern = map(peak_gene_correlation_group_chromosome_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_normalized_aggregate_matrices.WNN,
    description = "Normalize donor-state pseudobulk GEX and ATAC matrices as log1p CPM [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = normalize_peak_gene_correlation_aggregate_matrices(
      peak_gene_correlation_aggregate_matrices.WNN
    ),
    pattern = map(peak_gene_correlation_aggregate_matrices.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_results_tibbles.WNN,
    description = "Score donor-adjusted peak-gene associations for one broad cell type and chromosome [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = score_peak_gene_correlations_for_cell_group(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN
    ),
    pattern = map(peak_gene_correlation_normalized_aggregate_matrices.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = peak_gene_correlation_branch_diagnostics_tibbles.WNN,
    description = "Diagnose skipped or retained peak-gene correlation branches [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = diagnose_peak_gene_correlation_branch(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN
    ),
    pattern = map(peak_gene_correlation_normalized_aggregate_matrices.WNN),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_results_tibble.WNN,
    description = "Combine chromosome-level exploratory peak-gene correlations [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = finalize_peak_gene_correlation_results(
      results_tibbles = peak_gene_correlation_results_tibbles.WNN,
      aggregation = aggregation
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_links_tibble.WNN,
    description = "Select nonpromoter candidates at adjusted r >= 0.15 and conditional aggregate FDR < 0.05 [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_links(peak_gene_correlation_results_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_finemapping_reference.WNN,
    description = "Summarize global FDR breakpoints and gene-feature mappings for chromosome-level SuSiE inputs [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_finemapping_reference(
      finalized_results_tibble = peak_gene_correlation_results_tibble.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_finemapped_links_tibbles.WNN,
    description = "Conditionally prioritize peaks for donor-adjusted candidate enhancer genes with SuSiE [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      reference <- peak_gene_correlation_finemapping_reference.WNN
      finemap_peak_gene_correlations_for_branch(
        normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
        candidate_pairs_tibble = reference$gene_features,
        finalized_results_tibble = restore_peak_gene_correlation_FDR(
          peak_gene_correlation_results_tibbles.WNN,
          reference$FDR
        ) |>
          dplyr::arrange(.data$cell_group, .data$TargetGeneID, dplyr::desc(.data$correlation))
      )
    },
    pattern = map(
      peak_gene_correlation_normalized_aggregate_matrices.WNN,
      peak_gene_correlation_results_tibbles.WNN
    ),
    iteration = "vector",
    packages = w_def("susieR"),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_finemapped_links_tibble.WNN,
    description = "Combine compact SuSiE peak-gene prioritization records [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = dplyr::bind_rows(
      peak_gene_correlation_finemapped_links_tibbles.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_diagnostics_tibble.WNN,
    description = "Combine peak-gene correlation cell-group and branch diagnostics [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = dplyr::bind_rows(
      peak_gene_correlation_cell_group_diagnostics_tibble.WNN,
      peak_gene_correlation_donor_state_diagnostics_tibble.WNN,
      peak_gene_correlation_branch_diagnostics_tibbles.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_correlation_histogram_plot,
    description = "Save a facetted histogram of peak-gene correlations by cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      correlation_bin_width <- 0.025
      correlation_plot_tibble <- summarize_peak_gene_correlation_histogram(
        results_tibble = peak_gene_correlation_results_tibble.WNN,
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
    name = peak_gene_correlation_support_counts_plot,
    description = "Save peak-gene correlation tested, significant, and linked pair counts by cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      support_plot_tibble <- summarize_peak_gene_correlation_support_counts(
        peak_gene_correlation_results_tibble.WNN
      )
      plot_peak_gene_correlation_support_counts(support_plot_tibble) |>
        save_plots_structured()
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_distance_correlation_plot,
    description = "Save median peak-gene correlation by absolute TSS distance and cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      distance_plot_tibble <- summarize_peak_gene_correlation_by_distance(
        peak_gene_correlation_results_tibble.WNN
      )
      plot_peak_gene_correlation_by_distance(distance_plot_tibble) |>
        save_plots_structured()
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_diagnostics_plot,
    description = "Save peak-gene correlation retained and skipped branch diagnostics by cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      diagnostics_plot_tibble <- peak_gene_correlation_diagnostics_tibble.WNN |>
        dplyr::mutate(skipped_reason = dplyr::coalesce(.data$skipped_reason, "retained")) |>
        dplyr::count(.data$cell_group, .data$skipped_reason, name = "n_branches")

      diagnostics_plot <- ggplot2::ggplot(
        diagnostics_plot_tibble,
        ggplot2::aes(x = .data$n_branches, y = .data$cell_group, fill = .data$skipped_reason)
      ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = "Peak-gene analysis coverage and skipped groups",
          subtitle = stringr::str_wrap("Check skipped analyses before interpreting sparse results; insufficient data are not evidence of absent associations.", width = 100),
          caption = stringr::str_wrap("Bars count analysis branches by reported status within each cell group. A retained branch completed eligibility checks; it does not imply a significant link. Counts refer to computational branches, not cells, donors or peak-gene pairs.", width = 110),
          x = "Branches",
          y = "Cell group",
          fill = "Status"
        )

      save_plots_structured(diagnostics_plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_links_tibble.WNN,
    description = "Split top peak-gene links per cell group for downstream QC plots [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_top_links(
      links_tibble = peak_gene_correlation_links_tibble.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN,
      n_per_cell_group = peak_gene_correlation_top_links_per_cell_group
    ),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_aggregate_values_tibbles.WNN,
    description = "Extract top-link aggregate values from one normalized cell-group/chromosome branch [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = extract_peak_gene_correlation_top_link_aggregate_values(
      normalized_aggregate_matrices =
        peak_gene_correlation_normalized_aggregate_matrices.WNN,
      top_links_tibble =
        peak_gene_correlation_top_links_tibble.WNN
    ),
    pattern = map(
      peak_gene_correlation_normalized_aggregate_matrices.WNN
    ),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_aggregate_values_tibble.WNN,
    description = "Combine compact top-link aggregate values across cell groups and chromosomes [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = dplyr::bind_rows(
      peak_gene_correlation_top_link_aggregate_values_tibbles.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_aggregate_scatter_tibble.WNN,
    description = "Build aggregate-level GEX and ATAC values for one top peak-gene link [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      link_row <- peak_gene_correlation_top_links_tibble.WNN
      if (!isTRUE(link_row$is_analyzable_link[[1]])) {
        link_row |>
          dplyr::select("scatter_plot_name", "is_analyzable_link")
      } else {
        peak_gene_correlation_top_link_aggregate_values_tibble.WNN |>
          dplyr::filter(
            .data$scatter_plot_name == link_row$scatter_plot_name[[1]]
          )
      }
    },
    pattern = map(peak_gene_correlation_top_links_tibble.WNN),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 32) # apparently 16 GB is not enough
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_top_link_aggregate_scatter_plots,
    description = "Save aggregate-level scatterplots for top peak-gene links per cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      plot_tibble <- peak_gene_correlation_top_link_aggregate_scatter_tibble.WNN
      if (!isTRUE(plot_tibble$is_analyzable_link[[1]])) {
        return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links") |>
          save_plots_structured(
            override_suffix = plot_tibble$scatter_plot_name[[1]],
            dyn_suffix_in_subdir = TRUE
          ))
      }
      scatter_plot <- plot_peak_gene_correlation_aggregate_scatter(plot_tibble)

      genome_annotation_track <- peak_gene_correlation_top_link_genome_annotation_track.WNN +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank()
        )
      primary_ATAC_track <- peak_gene_correlation_top_link_primary_cell_ATAC_accessibility_track.WNN +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank()
        )

      track_panel <- patchwork::wrap_plots(
        genome_annotation_track,
        primary_ATAC_track,
        peak_gene_correlation_top_link_peak_gene_loop_track.WNN,
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
    })(),
    pattern = map(
      peak_gene_correlation_top_link_aggregate_scatter_tibble.WNN,
      peak_gene_correlation_top_link_genome_annotation_track.WNN,
      peak_gene_correlation_top_link_primary_cell_ATAC_accessibility_track.WNN,
      peak_gene_correlation_top_link_peak_gene_loop_track.WNN
    )
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_loci_tibble.WNN,
    description = "Add plotting windows to top peak-gene links for per-locus track plots [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      link_tibble <- peak_gene_correlation_top_links_tibble.WNN
      if (!isTRUE(link_tibble$is_analyzable_link[[1]])) {
        link_tibble |>
          dplyr::mutate(locus_start = 1L, locus_end = 1L)
      } else {
        link_tibble |>
          dplyr::mutate(
            locus_start = pmax(
              1L,
              min(.data$start[[1]], .data$TargetGeneTSS[[1]], na.rm = TRUE) - 50000L
            ),
            locus_end = max(.data$end[[1]], .data$TargetGeneTSS[[1]], na.rm = TRUE) + 50000L
          )
      }
    },
    pattern = map(peak_gene_correlation_top_links_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_genome_annotation_track.WNN,
    description = "Build the genome annotation track for one top peak-gene correlation locus [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.WNN
      if (!isTRUE(locus_tibble$is_analyzable_link[[1]])) {
        return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links"))
      }
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
    })(),
    pattern = map(peak_gene_correlation_top_link_loci_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_ATAC_coverage_tibble.WNN,
    description = "Compute per-cell-group BPCells ATAC coverage once for one top peak-gene correlation locus [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.WNN
      if (!isTRUE(locus_tibble$is_analyzable_link[[1]])) {
        return(tibble::tibble())
      }
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

      BPCells::trackplot_coverage(
        fragments = fragments,
        region = region,
        groups = groups,
        cell_read_counts = cell_read_counts,
        group_order = group_order,
        colors = grDevices::hcl.colors(length(group_order), palette = "Dark 3"),
        bins = 500,
        return_data = TRUE
      )
    })(),
    pattern = map(peak_gene_correlation_top_link_loci_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_primary_cell_ATAC_accessibility_track.WNN,
    description = "Build the primary-cell-group ATAC coverage track from shared per-locus coverage [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.WNN
      if (!isTRUE(locus_tibble$is_analyzable_link[[1]])) {
        return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links"))
      }
      primary_cell_group <- locus_tibble$cell_group[[1]]
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      primary_coverage_tibble <-
        peak_gene_correlation_top_link_ATAC_coverage_tibble.WNN |>
        dplyr::filter(as.character(.data$group) == primary_cell_group) |>
        dplyr::mutate(
          group = factor(as.character(.data$group), levels = primary_cell_group)
        ) |>
        droplevels()
      make_BPCells_ATAC_coverage_track_from_tibble(
        coverage_tibble = primary_coverage_tibble,
        region = region
      )
    })(),
    pattern = map(
      peak_gene_correlation_top_link_loci_tibble.WNN,
      peak_gene_correlation_top_link_ATAC_coverage_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_ATAC_accessibility_track.WNN,
    description = "Build the per-cell-group ATAC coverage track from shared per-locus coverage [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.WNN
      if (!isTRUE(locus_tibble$is_analyzable_link[[1]])) {
        return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links"))
      }
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      make_BPCells_ATAC_coverage_track_from_tibble(
        coverage_tibble =
          peak_gene_correlation_top_link_ATAC_coverage_tibble.WNN,
        region = region
      )
    })(),
    pattern = map(
      peak_gene_correlation_top_link_loci_tibble.WNN,
      peak_gene_correlation_top_link_ATAC_coverage_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_peak_gene_loop_track.WNN,
    description = "Build the BPCells peak-gene loop track for filtered links at one top-link locus [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.WNN
      if (!isTRUE(locus_tibble$is_analyzable_link[[1]])) {
        return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links"))
      }
      region <- GenomicRanges::GRanges(
        seqnames = locus_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_tibble$locus_start[[1]],
          end = locus_tibble$locus_end[[1]]
        )
      )

      loops_tibble <- peak_gene_correlation_links_tibble.WNN |>
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
    })(),
    pattern = map(peak_gene_correlation_top_link_loci_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = peak_gene_correlation_top_link_ATAC_tracks_plots,
    description = "Save combined genome annotation, ATAC coverage, and peak-gene loop tracks at top-link loci [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = (function() {
      locus_tibble <- peak_gene_correlation_top_link_loci_tibble.WNN
      if (!isTRUE(locus_tibble$is_analyzable_link[[1]])) {
        return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links") |>
          save_plots_structured(
            override_suffix = locus_tibble$scatter_plot_name[[1]],
            dyn_suffix_in_subdir = TRUE
          ))
      }
      plot <- BPCells::trackplot_combine(
        tracks = list(
          peak_gene_correlation_top_link_genome_annotation_track.WNN,
          peak_gene_correlation_top_link_ATAC_accessibility_track.WNN,
          peak_gene_correlation_top_link_peak_gene_loop_track.WNN
        ),
        title = paste(
          peak_gene_correlation_top_link_loci_tibble.WNN$TargetGene[[1]],
          peak_gene_correlation_top_link_loci_tibble.WNN$peak[[1]],
          peak_gene_correlation_top_link_loci_tibble.WNN$cell_group[[1]],
          sep = " - "
        )
      ) + patchwork::plot_annotation(
        subtitle = stringr::str_wrap("Read candidate loops alongside local accessibility; arcs represent statistical associations, not measured chromatin contacts.", width = 100),
        caption = stringr::str_wrap("Loops connect peak midpoints to gene TSSs for candidate links in this window (adjusted r >= 0.15, conditional BH FDR < 0.05, excluding self-promoters). Coverage uses 500 bins normalized by bin width and group depth, with extremes clipped at the 99.9th percentile. Groups use GEX-derived cell types; pooled coverage does not establish donor replication or causality.", width = 110))

      plot |>
        save_plots_structured(
          override_suffix = peak_gene_correlation_top_link_loci_tibble.WNN$scatter_plot_name[[1]],
          dyn_suffix_in_subdir = TRUE
        )
    })(),
    pattern = map(
      peak_gene_correlation_top_link_loci_tibble.WNN,
      peak_gene_correlation_top_link_genome_annotation_track.WNN,
      peak_gene_correlation_top_link_ATAC_accessibility_track.WNN,
      peak_gene_correlation_top_link_peak_gene_loop_track.WNN
    )
  )
)

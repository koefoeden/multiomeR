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
    command = purrr::map_dfr(peak_gene_correlation_donor_state_records.WNN, "aggregates"),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_donor_state_diagnostics_tibble.WNN,
    description = "Combine donor-state pseudobulk eligibility diagnostics [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = purrr::map_dfr(peak_gene_correlation_donor_state_records.WNN, "diagnostics"),
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
    name = peak_gene_correlation_filter_records.WNN,
    description = "Filter candidate hypotheses by counts and donor support and cache retention diagnostics [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = filter_peak_gene_candidate_pairs(
      aggregate_matrices = peak_gene_correlation_aggregate_matrices.WNN,
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
      candidate_pairs_tibble = peak_gene_correlation_candidate_pairs_tibble.WNN,
      filter = peak_gene_correlation_filter
    ),
    pattern = map(peak_gene_correlation_aggregate_matrices.WNN,
      peak_gene_correlation_normalized_aggregate_matrices.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_filter_diagnostics_tibble.WNN,
    description = "Combine retention and overlapping exclusion counts for all support presets [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = purrr::map_dfr(peak_gene_correlation_filter_records.WNN, "diagnostics"),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = filter_retention_plot,
    description = "Compare candidate retention across measurement-support filters by cell type [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = plot_peak_gene_filter_retention(
      peak_gene_correlation_filter_diagnostics_tibble.WNN, peak_gene_correlation_filter
    ) |> save_plots_structured(width = 13, height = 9)
  ),
  targets::tar_target(
    name = peak_gene_correlation_branch_diagnostics_tibbles.WNN,
    description = "Diagnose skipped or retained peak-gene correlation branches [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = diagnose_peak_gene_correlation_branch(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
      candidate_pairs_tibble = peak_gene_correlation_filter_records.WNN$candidate_pairs
    ),
    pattern = map(peak_gene_correlation_normalized_aggregate_matrices.WNN,
      peak_gene_correlation_filter_records.WNN),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_links_tibble.WNN,
    description = "Select positive, reliable nonpromoter pairs at hierarchical FDR < 0.05 [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_correlation_links(peak_gene_correlation_hierarchical_results_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_finemapping_reference.WNN,
    description = "Summarize within-cell-group FDR breakpoints for chromosome-level SuSiE inputs [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = make_peak_gene_finemapping_reference(peak_gene_correlation_hierarchical_results_tibble.WNN),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = peak_gene_correlation_finemapped_links_tibbles.WNN,
    description = "Conditionally prioritize peaks for candidate enhancer genes with SuSiE [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = finemap_peak_gene_correlations_for_branch(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
      branch_results = restore_peak_gene_correlation_FDR(
        peak_gene_correlation_hierarchical_results_tibbles.WNN,
        peak_gene_correlation_finemapping_reference.WNN
      )
    ),
    pattern = map(
      peak_gene_correlation_normalized_aggregate_matrices.WNN,
      peak_gene_correlation_hierarchical_results_tibbles.WNN
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
    name = correlation_histogram_plot,
    description = "Save a facetted histogram of peak-gene correlations by cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      correlation_bin_width <- 0.025
      correlation_plot_tibble <- summarize_peak_gene_correlation_histogram(
        results_tibble = peak_gene_correlation_hierarchical_results_tibble.WNN,
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
    name = support_counts_plot,
    description = "Save peak-gene correlation tested, significant, and linked pair counts by cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      support_plot_tibble <- summarize_peak_gene_correlation_support_counts(
        peak_gene_correlation_hierarchical_results_tibble.WNN
      )
      plot_peak_gene_correlation_support_counts(support_plot_tibble) |>
        save_plots_structured()
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = significant_pairs_technical_features_tibble.WNN,
    description = "Summarize peak-gene discovery counts and technical features across WNN cell types [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = prepare_peak_gene_support_technical_features(
      results_tibble = peak_gene_correlation_hierarchical_results_tibble.WNN,
      donor_state_aggregates_tibble = peak_gene_correlation_donor_state_aggregates_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = significant_pairs_vs_technical_features_plot,
    description = "Compare FDR-significant peak-gene pairs with sampling and measurement features across WNN cell types [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = plot_peak_gene_significant_pairs_vs_technical_features(
      significant_pairs_technical_features_tibble.WNN
    ) |>
      save_plots_structured(width = 24, height = 17)
  ),
  tarchetypes::tar_file(
    name = distance_correlation_plot,
    description = "Save median peak-gene correlation and -log10(p) by absolute TSS distance and cell group [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = {
      distance_plot_tibble <- summarize_peak_gene_correlation_by_distance(
        peak_gene_correlation_hierarchical_results_tibble.WNN
      )
      plot_peak_gene_correlation_by_distance(distance_plot_tibble) |>
        save_plots_structured(width = 14, height = 15)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = diagnostics_plot,
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
  tarchetypes::tar_file(
    name = peak_gene_KR_native_source_file,
    description = "Track the batched Kenward-Roger C++ source [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = "src/peak_gene_KR.cpp",
    deployment = "main"
  ),
  tarchetypes::tar_file(
    name = peak_gene_REML_native_source_file,
    description = "Track the profiled donor-slope REML C++ source [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = "src/peak_gene_REML.cpp",
    deployment = "main"
  ),
  targets::tar_target(
    name = peak_gene_correlation_hierarchical_results_tibbles.WNN,
    description = "Scan all eligible pairs with donor-varying slopes and Kenward-Roger inference [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = score_peak_gene_hierarchical_associations(
      normalized_aggregate_matrices = peak_gene_correlation_normalized_aggregate_matrices.WNN,
      candidate_pairs_tibble = peak_gene_correlation_filter_records.WNN$candidate_pairs,
      REML_source_file = peak_gene_REML_native_source_file,
      KR_source_file = peak_gene_KR_native_source_file
    ),
    pattern = map(peak_gene_correlation_normalized_aggregate_matrices.WNN,
      peak_gene_correlation_filter_records.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_hierarchical_results_tibble.WNN,
    description = "Combine full-scan hierarchical fits and adjust p-values within cell type [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = finalize_peak_gene_hierarchical_results(peak_gene_correlation_hierarchical_results_tibbles.WNN),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_links_tibble.WNN,
    description = "Select lowest reliable hierarchical p-values across the full scan for QC plots [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = select_peak_gene_hierarchical_top_links(
      hierarchical_results = peak_gene_correlation_hierarchical_results_tibble.WNN,
      n_per_cell_group = peak_gene_correlation_top_links_per_cell_group
    ),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 60)
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
  targets::tar_target(
    name = peak_gene_correlation_top_gene_context_tibbles.WNN,
    description = "Extract selected-gene context from each hierarchical scan branch [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = dplyr::semi_join(
      peak_gene_correlation_hierarchical_results_tibbles.WNN,
      peak_gene_correlation_top_links_tibble.WNN, by = "TargetGeneID"
    ),
    pattern = map(peak_gene_correlation_hierarchical_results_tibbles.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_gene_context_tibble.WNN,
    description = "Combine compact hierarchical gene context with peak coordinates [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = dplyr::bind_rows(peak_gene_correlation_top_gene_context_tibbles.WNN),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_plot_data.WNN,
    description = "Prepare full-window genes, peaks and hierarchical evidence for one focal link [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = prepare_peak_gene_top_link_plot_data(
      peak_gene_correlation_top_links_tibble.WNN,
      peak_gene_correlation_top_gene_context_tibble.WNN,
      marker_validated_Ensembl_annotations_GRanges_list$genes,
      consensus_peak_GRanges.ATAC
    ),
    pattern = map(peak_gene_correlation_top_links_tibble.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = peak_gene_correlation_top_link_focal_ATAC_coverage_tibble.WNN,
    description = "Cache focal-cell-type insertion coverage over the full candidate window [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = compute_peak_gene_focal_coverage(
      peak_gene_correlation_top_link_plot_data.WNN,
      metadata_w_cell_types_tibble.WNN, combined_BPCells_fragment_obj.ATAC
    ),
    pattern = map(peak_gene_correlation_top_link_plot_data.WNN),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = top_link_aggregate_scatter_plots,
    description = "Save hierarchical genomic evidence and donor-residual scatter for each top link [checkpoint:peak_gene_correlation] [part_of_graph:peak_gene_correlation]",
    command = plot_peak_gene_correlation_top_link(
      peak_gene_correlation_top_link_plot_data.WNN,
      peak_gene_correlation_top_link_focal_ATAC_coverage_tibble.WNN,
      peak_gene_correlation_top_link_aggregate_scatter_tibble.WNN,
      sort(unique(metadata_w_cell_types_tibble.WNN$WNN_harmony_SNN_cluster_cell_type))
    ) |> save_plots_structured(
      override_suffix = peak_gene_correlation_top_link_plot_data.WNN$link$scatter_plot_name[[1]],
      dyn_suffix_in_subdir = TRUE, width = 15, height = 19.5
    ),
    pattern = map(peak_gene_correlation_top_link_plot_data.WNN,
      peak_gene_correlation_top_link_focal_ATAC_coverage_tibble.WNN,
      peak_gene_correlation_top_link_aggregate_scatter_tibble.WNN)
  )
)

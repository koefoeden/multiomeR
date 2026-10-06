# Manuscript figures and tables from configuration_dev's public-data aggregations:
# Figure 1 from mixed_human_31x, Supplementary Figure S2 comparing multiomeR with a
# conventional Seurat/Signac workflow on comparison_1x to comparison_20x, and
# Supplementary Table S4 comparing the algorithm reimplementations with their
# references on comparison_5x. See README.md.
source("pipeline_manuscript_figures/helpers.R")
source("pipeline_manuscript_figures/algorithm_parity_helpers.R")
source("pipeline_manuscript_figures/seurat_signac_helpers.R")
base_pipeline <- source("_targets.R")$value
output_dir <- "pipeline_manuscript_figures/results"

figure_pipeline <- list(
  targets::tar_target(
    figure_1_panel_A.mixed_human_31x,
    description = "Figure 1A: strongest ranked hierarchical link; compact focal genomic tracks from top_link_aggregate_scatter_plots data",
    command = local({
      link <- select_figure_1_link(peak_gene_correlation_top_links_tibble.WNN.peak_gene_correlation.mixed_human_31x)
      records <- peak_gene_correlation_top_link_plot_data.WNN.peak_gene_correlation.mixed_human_31x
      index <- which(vapply(records, function(x) identical(x$link$scatter_plot_name, link$scatter_plot_name), logical(1)))
      stopifnot(length(index) == 1L)
      plot_figure_1_link_tracks(records[[index]],
        peak_gene_correlation_top_link_focal_ATAC_coverage_tibble.WNN.peak_gene_correlation.mixed_human_31x[[index]]) |>
        compact_figure_1_plot()
    }),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_B.mixed_human_31x,
    description = "Figure 1B: same aggregate residuals as top_link_aggregate_scatter_plots, with inset hierarchical P and compact donor colours",
    command = local({
      link <- select_figure_1_link(peak_gene_correlation_top_links_tibble.WNN.peak_gene_correlation.mixed_human_31x)
      data <- peak_gene_correlation_top_link_aggregate_values_tibble.WNN.peak_gene_correlation.mixed_human_31x |>
        dplyr::filter(.data$scatter_plot_name == link$scatter_plot_name[[1]],
          .data$cell_group == .data$primary_cell_group)
      stopifnot(nrow(data) > 0L)
      plot_figure_1_scatter(data) |> compact_figure_1_plot()
    }),
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_C.mixed_human_31x,
    description = "Figure 1C: raw_deviation_unscaled data, up to six raw Z >= qnorm(0.95) traits, all cell types, no metadata tracks",
    command = local({
      scores <- chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x
      plot_figure_1_heatmap(scores, select_figure_1_traits(scores)) |> compact_figure_1_plot()
    }),
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_D.mixed_human_31x,
    description = "Figure 1D: variant-detail data for the top absolute-contribution locus in the two strongest raw Z combinations, left to right; protein-coding genes, focal coverage and PIP only",
    command = local({
      selections <- select_figure_1_combinations(chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x)
      locus_records <- lapply(seq_len(nrow(selections)), function(i) {
        selection <- selections[i, ]
        locus <- chromVAR_locus_contribution_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x |>
          dplyr::filter(.data$GWAS_ID == selection$GWAS_ID[[1]], .data$cluster == selection$cluster[[1]]) |>
          dplyr::arrange(dplyr::desc(abs(.data$relative_deviation_contribution)), .data$studyLocusId) |>
          dplyr::slice_head(n = 1)
        records <- purrr::flatten(chromVAR_variant_contribution_detail_plot_records.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x)
        records <- purrr::keep(records, function(x) {
          x$GWAS_ID == selection$GWAS_ID[[1]] &&
            x$variant_tibble$cluster[[1]] == selection$cluster[[1]] &&
            x$variant_tibble$studyLocusId[[1]] == locus$studyLocusId[[1]]
        })
        stopifnot(length(records) == 1L)
        record <- records[[1]]
        gene_ranges <- marker_validated_Ensembl_annotations_GRanges_list.mixed_human_31x$genes
        genes <- prepare_genomic_gene_bodies(
          gene_ranges[gene_ranges$gene_biotype == "protein_coding"], record$region)
        list(record = record, genes = genes, selection = selection)
      })
      plot_figure_1_loci(locus_records) |> compact_figure_1_plot()
    }),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_E.mixed_human_31x,
    description = "Figure 1E: rebuild TRS_UMAPs for the two traits in D, in matching order, with randomized cell draw order, cell-type centroids and a shared TRS scale",
    command = local({
      selections <- select_figure_1_combinations(chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x)
      selected_scores <- TRS_tibble.WNN_harmony_SNN.SCAVENGE.single_nucleus.genetic_enrichment.mixed_human_31x |>
        dplyr::filter(.data$GWAS_ID %in% selections$GWAS_ID)
      score_limits <- range(selected_scores$score, finite = TRUE)
      plots <- lapply(selections$GWAS_ID, function(trait) {
        scores <- selected_scores |>
          dplyr::filter(.data$GWAS_ID == trait) |>
          dplyr::select("barcode_w_prefix", "score")
        stopifnot(nrow(scores) > 0L)
        data <- metadata_w_cell_types_tibble.WNN.mixed_human_31x |>
          dplyr::select("barcode_w_prefix", "WNN_UMAP_1", "WNN_UMAP_2", "PCA_harmony_SNN_cluster_cell_type") |>
          dplyr::left_join(scores, by = "barcode_w_prefix", relationship = "one-to-one") |>
          dplyr::filter(is.finite(.data$WNN_UMAP_1), is.finite(.data$WNN_UMAP_2))
        # Match production random drawing; use the same permutation for both traits.
        data <- withr::with_seed(1L, dplyr::slice_sample(data, prop = 1))
        plot_figure_1_TRS(data, trait, score_limits) |> compact_figure_1_plot()
      })
      ((patchwork::wrap_plots(plots, nrow = 1, guides = "collect") +
        patchwork::plot_annotation(title = "SCAVENGE trait relevance",
          theme = figure_1_theme())) & ggplot2::theme(legend.position = "bottom")) |>
        compact_figure_1_plot()

    }),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_plot.mixed_human_31x,
    description = "Combine manuscript Figure 1 ggplot panels in four rows: AB, C, D, E",
    command = compose_figure_1(
      figure_1_panel_A.mixed_human_31x,
      figure_1_panel_B.mixed_human_31x,
      figure_1_panel_C.mixed_human_31x,
      figure_1_panel_D.mixed_human_31x,
      figure_1_panel_E.mixed_human_31x) |>
      compact_figure_1_plot(),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    figure_1.mixed_human_31x,
    description = "Export full-page manuscript Figure 1 from its combined ggplot object",
    command = {
      path <- file.path(output_dir, "1.png")
      ggplot2::ggsave(path, figure_1_plot.mixed_human_31x,
        device = ragg::agg_png, width = 8, height = 11, dpi = 300, create.dir = TRUE)
      path
    },
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    seurat_signac_comparison_resources,
    description = "Supplementary Figure S2 data: multiomeR and Seurat/Signac resource use from recorded targets metadata and Slurm job monitoring, one row per aggregation and workflow",
    command = {
      path <- file.path(output_dir, "data", "seurat_signac_comparison_resources.tsv")
      fs::dir_create(dirname(path))
      readr::write_tsv(seurat_signac_comparison_resource_tibble(
        paste0("comparison_", c(1, 2, 5, 10, 20), "x")), path)
      path
    },
    # targets lets only file targets read the store metadata, which no upstream dependency tracks.
    # A worker rebuilds the graph and reads the monitor history, so dispatching continues meanwhile.
    cue = targets::tar_cue(mode = "always"), resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    seurat_signac_comparison_resource_plot,
    description = "Supplementary Figure S2: input nuclei against each resource measure, faceted, coloured by workflow",
    command = readr::read_tsv(seurat_signac_comparison_resources, show_col_types = FALSE) |>
      plot_seurat_signac_comparison_resources() |> compact_figure_1_plot(),
    resources = get_tar_resources(RAM_GB_req = 1)
  ),
  tarchetypes::tar_file(
    supplementary_figure_S2,
    description = "Export Supplementary Figure S2 from its ggplot object",
    command = {
      path <- file.path(output_dir, "S2.png")
      ggplot2::ggsave(path, seurat_signac_comparison_resource_plot,
        device = ragg::agg_png, width = 6.5, height = 4.5, dpi = 300, create.dir = TRUE)
      path
    },
    resources = get_tar_resources(RAM_GB_req = 1)
  )
)

# Supplementary Table S4: native and reference implementations of each
# reimplemented algorithm on comparison_5x production objects, one target per
# implementation so the comparisons run in parallel.
algorithm_parity_run <- function(algorithm, implementation, runner, fixture, extra = list()) {
  targets::tar_target_raw(
    paste0("algorithm_parity_run_", algorithm, "_", implementation, ".comparison_5x"),
    command = as.call(c(quote(measure_algorithm_parity_case), list(algorithm, implementation, as.symbol(runner),
      as.symbol(paste0("algorithm_parity_fixture_", fixture, ".comparison_5x")),
      cores = 6L, package_versions = quote(algorithm_parity_package_versions)), extra)),
    description = paste("Supplementary Table S4:", implementation, algorithm,
      "on comparison_5x objects in a fresh process, with wall time and peak process-tree RSS"),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
  )
}
algorithm_parity_cases <- list(
  list("ucell", "UCell"), list("wnn", "WNN", list(native_source_file = quote(WNN_native_source_file))),
  list("scdblfinder_atac", "scDblFinder_ATAC"), list("scdblfinder_atac_predoublet", "scDblFinder_ATAC_predoublet"),
  list("amulet", "amulet", list(native_source_file = quote(amulet_BPCells_native_source_file))),
  list("scavenge", "SCAVENGE", list(), list(reference_file = quote(algorithm_parity_SCAVENGE_reference_file))))
algorithm_parity_runner <- c(ucell = "UCell", wnn = "WNN", scdblfinder_atac = "scDblFinder_ATAC",
  scdblfinder_atac_predoublet = "scDblFinder_ATAC", amulet = "amulet", scavenge = "SCAVENGE")

algorithm_parity_pipeline <- c(
  list(
    targets::tar_target(
      algorithm_parity_package_versions,
      description = "Installed versions of the packages the S4 comparisons run; rerun every time so the comparisons follow the lockfile",
      command = get_algorithm_parity_package_versions(),
      cue = targets::tar_cue(mode = "always"), deployment = "main"
    ),
    tarchetypes::tar_file(
      algorithm_parity_SCAVENGE_reference_file,
      description = "Pinned SCAVENGE reference functions shared with the SCAVENGE parity tests",
      command = "tests/testthat/fixtures/algorithm-deviations.R"
    ),
    targets::tar_target(
      algorithm_parity_fixture_UCell.comparison_5x,
      description = "S4 UCell fixture: pre-doublet GEX counts and the configured marker sets",
      command = get_UCell_parity_fixture(aggregated_counts_BPCells_matrix.GEX.comparison_5x,
        metadata_w_clusters_tibble.GEX.comparison_5x, UCell_GEX_marker_genes_list.comparison_5x),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      algorithm_parity_fixture_WNN.comparison_5x,
      description = "S4 WNN fixture: aligned GEX and ATAC Harmony embeddings",
      command = get_WNN_parity_fixture(embedding_matrices.WNN.comparison_5x),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      algorithm_parity_fixture_scDblFinder_ATAC.comparison_5x,
      description = "S4 ATAC scDblFinder fixture: one well's peak counts, global LSI feature groups and clusters",
      command = get_scDblFinder_ATAC_parity_fixture(scDblFinder_GEM_well_tibble.ATAC.comparison_5x,
        peak_QC_filtered_BPCells_matrix.ATAC.comparison_5x, scDblFinder_feature_groups.ATAC.comparison_5x,
        LSI_loadings_tibble.ATAC.comparison_5x, "comparison_pbmc_10k_chromium_x"),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      algorithm_parity_scDblFinder_ATAC_predoublet_matrix_dir.comparison_5x,
      description = "S4 pre-doublet ATAC fixture: insertion-counted consensus peaks of one well after ordinary QC",
      command = write_scDblFinder_ATAC_predoublet_parity_matrix(
        cellranger_kept_metadata_tibble.comparison_pbmc_10k_chromium_x,
        excluded_cellranger_only_barcodes_by_type_list.comparison_pbmc_10k_chromium_x,
        consensus_peak_GRanges.ATAC.comparison_5x, fragments_w_prefix_bpcells.comparison_pbmc_10k_chromium_x,
        dir = file.path(targets::tar_path_store(), "files", "pipeline_manuscript_figures",
          "algorithm_parity_scDblFinder_ATAC_predoublet_peak_matrix")),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      algorithm_parity_fixture_scDblFinder_ATAC_predoublet.comparison_5x,
      description = "S4 pre-doublet ATAC fixture with the global LSI feature groups and no clusters",
      command = get_scDblFinder_ATAC_predoublet_parity_fixture(
        open_BPCells_dir(algorithm_parity_scDblFinder_ATAC_predoublet_matrix_dir.comparison_5x),
        scDblFinder_feature_groups.ATAC.comparison_5x, "comparison_pbmc_10k_chromium_x"),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      algorithm_parity_fixture_amulet.comparison_5x,
      description = "S4 AMULET fixture: one well's prefixed BPCells fragments and original Cell Ranger fragment file",
      command = get_amulet_parity_fixture(fragments_w_prefix_bpcells.comparison_pbmc_unsorted_10k,
        cellranger_barcodes_tsv.comparison_pbmc_unsorted_10k, cellranger_summary_file.comparison_pbmc_unsorted_10k,
        "comparison_pbmc_unsorted_10k"),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      algorithm_parity_fixture_SCAVENGE.comparison_5x,
      description = "S4 SCAVENGE fixture: final WNN SNN graph and motif-family chromVAR Z-scores as seed signal",
      command = get_SCAVENGE_parity_fixture(WNN_results.comparison_5x, motif_family_chromVAR_results.ATAC.comparison_5x,
        "cluster_001"),
      resources = get_tar_resources(RAM_GB_req = 16)
    )
  ),
  unlist(lapply(algorithm_parity_cases, \(case) {
    runner <- algorithm_parity_runner[[case[[1]]]]
    list(
      algorithm_parity_run(case[[1]], "native", paste0("run_", runner, "_parity_native"), case[[2]],
        if (length(case) >= 3) case[[3]] else list()),
      algorithm_parity_run(case[[1]], "reference", paste0("run_", runner, "_parity_reference"), case[[2]],
        if (length(case) >= 4) case[[4]] else list()))
  }), recursive = FALSE),
  list(
    targets::tar_target(
      algorithm_parity_checks.comparison_5x,
      description = "S4 result agreement: exact-parity or similarity metrics for each algorithm against its threshold",
      command = dplyr::bind_rows(
        compare_UCell_parity(algorithm_parity_run_ucell_native.comparison_5x$result,
          algorithm_parity_run_ucell_reference.comparison_5x$result),
        compare_WNN_parity(algorithm_parity_run_wnn_native.comparison_5x$result,
          algorithm_parity_run_wnn_reference.comparison_5x$result),
        compare_scDblFinder_ATAC_parity(algorithm_parity_run_scdblfinder_atac_native.comparison_5x$result,
          algorithm_parity_run_scdblfinder_atac_reference.comparison_5x$result, "scdblfinder_atac"),
        compare_scDblFinder_ATAC_parity(
          algorithm_parity_run_scdblfinder_atac_predoublet_native.comparison_5x$result,
          algorithm_parity_run_scdblfinder_atac_predoublet_reference.comparison_5x$result,
          "scdblfinder_atac_predoublet"),
        compare_amulet_parity(algorithm_parity_run_amulet_native.comparison_5x$result,
          algorithm_parity_run_amulet_reference.comparison_5x$result),
        compare_SCAVENGE_parity(algorithm_parity_run_scavenge_native.comparison_5x$result,
          algorithm_parity_run_scavenge_reference.comparison_5x$result)),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      algorithm_parity_scDblFinder_agreement.comparison_5x,
      description = "S4 ATAC scDblFinder agreement: shared cells, doublet calls, score-rank correlation and class agreement",
      command = dplyr::bind_rows(
        summarize_scDblFinder_ATAC_parity(algorithm_parity_run_scdblfinder_atac_native.comparison_5x$result,
          algorithm_parity_run_scdblfinder_atac_reference.comparison_5x$result, "scdblfinder_atac"),
        summarize_scDblFinder_ATAC_parity(
          algorithm_parity_run_scdblfinder_atac_predoublet_native.comparison_5x$result,
          algorithm_parity_run_scdblfinder_atac_predoublet_reference.comparison_5x$result,
          "scdblfinder_atac_predoublet")),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      supplementary_table_S4.comparison_5x,
      description = "Supplementary Table S4: measurements, parity checks and agreement as source TSVs and the Typst table the supplement includes; stops if a parity threshold fails",
      command = {
        runs <- list(
          algorithm_parity_run_ucell_native.comparison_5x,
          algorithm_parity_run_ucell_reference.comparison_5x,
          algorithm_parity_run_wnn_native.comparison_5x,
          algorithm_parity_run_wnn_reference.comparison_5x,
          algorithm_parity_run_scdblfinder_atac_native.comparison_5x,
          algorithm_parity_run_scdblfinder_atac_reference.comparison_5x,
          algorithm_parity_run_scdblfinder_atac_predoublet_native.comparison_5x,
          algorithm_parity_run_scdblfinder_atac_predoublet_reference.comparison_5x,
          algorithm_parity_run_amulet_native.comparison_5x,
          algorithm_parity_run_amulet_reference.comparison_5x,
          algorithm_parity_run_scavenge_native.comparison_5x,
          algorithm_parity_run_scavenge_reference.comparison_5x)
        checks <- algorithm_parity_checks.comparison_5x
        agreement <- algorithm_parity_scDblFinder_agreement.comparison_5x
        summary <- summarize_algorithm_parity(runs, checks)
        dir <- file.path(output_dir, "data", "algorithm_parity")
        fs::dir_create(dir)
        paths <- file.path(dir, c("measurements.tsv", "parity.tsv", "scdblfinder_atac_agreement.tsv", "summary.tsv"))
        readr::write_tsv(dplyr::bind_rows(lapply(runs, `[[`, "measurement")), paths[[1]])
        readr::write_tsv(checks, paths[[2]])
        readr::write_tsv(agreement, paths[[3]])
        readr::write_tsv(summary, paths[[4]])
        if (!all(checks$passed)) {
          stop("S4 parity thresholds failed: ", paste(checks$algorithm[!checks$passed], checks$metric[!checks$passed],
            collapse = ", "), call. = FALSE)
        }
        table_path <- file.path(output_dir, "S4.md")
        writeLines(format_algorithm_parity_table(summary, checks, agreement), table_path)
        c(paths, table_path)
      },
      resources = get_tar_resources(RAM_GB_req = 4)
    )
  )
)

# The conventional Seurat/Signac workflow. Its default packages and resources are
# set after the other targets are defined, so they apply to this chain only.
comparison_resources <- get_tar_resources(cores_req = 6, RAM_GB_req = 60)

targets::tar_option_set(
  packages = c("Seurat", "SeuratObject", "Signac"),
  resources = comparison_resources
)

comparison_aggregation_tibble <- tibble::tibble(
  aggregation = comparison_aggregation_names(),
  comparison_peak_calling_fragment_files = rlang::syms(
    paste0("fragments_per_peak_calling_cluster_discovery.fragments.ATAC.", aggregation)
  ),
  comparison_peak_calling_cluster_names = rlang::syms(paste0("peak_calling_cluster_names.ATAC.", aggregation)),
  comparison_blacklist_GRanges = rlang::syms(paste0("blacklist_GRanges.ATAC.", aggregation))
)

comparison_seurat_signac_targets <- function(values, object_resources, peak_resources) {
  tarchetypes::tar_map(
    values = values,
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    targets::tar_target(
      comparison_seurat_signac_config,
      {
        comparison_seurat_signac_config_files
        comparison_build_config(aggregation)
      }
    ),
    targets::tar_target(
      comparison_seurat_signac_cellranger_input_files,
      comparison_input_files(comparison_seurat_signac_config),
      format = "file"
    ),
    targets::tar_target(
      comparison_seurat_signac_input_multimodal_object,
      {
        comparison_seurat_signac_cellranger_input_files
        comparison_read_inputs(comparison_seurat_signac_config)
      },
      resources = object_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_RNA_preprocessed_object,
      comparison_preprocess_RNA(
        comparison_seurat_signac_input_multimodal_object,
        comparison_seurat_signac_config
      ),
      resources = object_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_RNA_clustered_object,
      comparison_cluster_RNA(
        comparison_seurat_signac_RNA_preprocessed_object,
        comparison_seurat_signac_config
      ),
      resources = object_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_consensus_peak_GRanges,
      comparison_call_consensus_peaks(
        comparison_peak_calling_fragment_files,
        comparison_peak_calling_cluster_names,
        comparison_seurat_signac_config,
        comparison_blacklist_GRanges
      ),
      resources = peak_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_requantified_multimodal_object,
      comparison_rebuild_ATAC_assay(
        comparison_seurat_signac_RNA_clustered_object,
        comparison_seurat_signac_consensus_peak_GRanges
      ),
      resources = peak_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_ATAC_preprocessed_object,
      {
        ATAC_object <- comparison_preprocess_ATAC(
          comparison_seurat_signac_requantified_multimodal_object,
          comparison_seurat_signac_config
        )
        if (is.null(ATAC_object[["lsi"]])) {
          stop("ATAC preprocessing did not create the expected lsi reduction.", call. = FALSE)
        }
        ATAC_object
      },
      resources = object_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_object,
      comparison_run_WNN(
        comparison_seurat_signac_ATAC_preprocessed_object,
        comparison_seurat_signac_config
      ),
      resources = object_resources
    ),
    targets::tar_target(
      comparison_seurat_signac_summary,
      comparison_summarize(
        comparison_seurat_signac_object,
        comparison_seurat_signac_config,
        comparison_seurat_signac_consensus_peak_GRanges
      )
    )
  )
}

comparison_seurat_signac_config_paths <- c(
  configuration_path("cfg_aggregations.yaml"), configuration_path("cfg_GEM_wells.tsv"),
  "cfg_pipeline_parameters.tsv"
)

comparison_seurat_signac_pipeline <- rlang::list2(
  targets::tar_target(
    comparison_seurat_signac_config_files,
    comparison_seurat_signac_config_paths,
    format = "file"
  ),
  comparison_seurat_signac_targets(
    values = dplyr::filter(comparison_aggregation_tibble, .data$aggregation %in% c("comparison_1x", "comparison_2x")),
    object_resources = comparison_resources,
    peak_resources = comparison_resources
  ),
  comparison_seurat_signac_targets(
    values = dplyr::filter(comparison_aggregation_tibble, .data$aggregation == "comparison_5x"),
    object_resources = comparison_resources,
    peak_resources = get_tar_resources(cores_req = 6, RAM_GB_req = 200)
  ),
  comparison_seurat_signac_targets(
    values = dplyr::filter(comparison_aggregation_tibble,
      .data$aggregation %in% c("comparison_10x", "comparison_20x")),
    object_resources = get_tar_resources(cores_req = 6, RAM_GB_req = 500),
    peak_resources = get_tar_resources(cores_req = 6, RAM_GB_req = 500)
  )
)

c(base_pipeline, figure_pipeline, algorithm_parity_pipeline, comparison_seurat_signac_pipeline)

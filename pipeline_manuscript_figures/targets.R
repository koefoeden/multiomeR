# Manuscript figures and tables from configuration_dev's public-data aggregations:
# Figure 1, Supplementary Figure S3 and Supplementary Table S3 from mixed_human_31x,
# Supplementary Figure S2 and Table S6 on what multiomeR costs and saves: against a
# conventional Seurat/Signac workflow on comparison_1x to comparison_31x, with its defaults
# on sizing_1x to sizing_91x and for revisions of revision_31x, and Supplementary Table S4 comparing
# the algorithm reimplementations with their references on comparison_5x. See README.md.
source("pipeline_manuscript_figures/helpers.R")
source("pipeline_manuscript_figures/algorithm_parity_helpers.R")
source("pipeline_manuscript_figures/seurat_signac_helpers.R")
source("pipeline_manuscript_figures/benchmark_helpers.R")
base_pipeline <- source("_targets.R")$value
output_dir <- "pipeline_manuscript_figures/results"

figure_pipeline <- list(
  targets::tar_target(
    figure_1_panel_A.mixed_human_31x,
    description = "Figure 1A: per WNN cell type, the configured marker with the most specific expression and promoter accessibility in donor-summed cell-type pseudobulks",
    command = figure_1_marker_profiles(UCell_WNN_marker_genes_list.mixed_human_31x,
      pseudobulk_counts_matrix.GEX.mixed_human_31x, pseudobulk_counts_matrix.ATAC.mixed_human_31x,
      consensus_peak_annotated_GRanges.ATAC.mixed_human_31x) |>
      plot_figure_1_markers() |> compact_figure_1_plot(),
    packages = w_def("BPCells"), resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    figure_1_panel_B.mixed_human_31x,
    description = "Figure 1B: strongest ranked hierarchical link; compact focal genomic tracks from top_link_aggregate_scatter_plots data",
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
    figure_1_panel_C.mixed_human_31x,
    description = "Figure 1C: same aggregate residuals as top_link_aggregate_scatter_plots, with inset hierarchical P and compact donor colours",
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
    figure_1_traits.mixed_human_31x,
    description = "Figure 1D-F traits: broadly supported cell-type enrichments (q < 0.05, top locus below half of the signal), one per study and per cell type and trait category",
    command = select_figure_1_traits(chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x,
      chromVAR_locus_contribution_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x),
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_D.mixed_human_31x,
    description = "Figure 1D: raw_deviation_unscaled data of the selected traits, all assigned cell types in configured order, traits as a staircase, no metadata tracks",
    command = plot_figure_1_heatmap(chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x,
      figure_1_traits.mixed_human_31x, names(UCell_WNN_marker_genes_list.mixed_human_31x)) |> compact_figure_1_plot(),
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_E.mixed_human_31x,
    description = "Figure 1E: variant-detail data for the top absolute-contribution locus of the showcase trait; protein-coding genes, focal coverage and PIP only",
    command = local({
      selections <- select_figure_1_showcase(figure_1_traits.mixed_human_31x)
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
        list(record = record, genes = genes, selection = selection, locus = locus)
      })
      plot_figure_1_loci(locus_records) |> compact_figure_1_plot()
    }),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_panel_F.mixed_human_31x,
    description = "Figure 1F: rebuild the TRS_UMAP of the trait in E, with randomized cell draw order and the cell types enriched for it labelled",
    command = local({
      trait <- select_figure_1_showcase(figure_1_traits.mixed_human_31x)$GWAS_ID
      scores <- TRS_tibble.WNN_harmony_SNN.SCAVENGE.single_nucleus.genetic_enrichment.mixed_human_31x |>
        dplyr::filter(.data$GWAS_ID == trait) |>
        dplyr::select("barcode_w_prefix", "score")
      stopifnot(nrow(scores) > 0L)
      data <- metadata_w_cell_types_tibble.WNN.mixed_human_31x |>
        dplyr::select("barcode_w_prefix", "WNN_UMAP_1", "WNN_UMAP_2", "WNN_harmony_SNN_cluster_cell_type") |>
        dplyr::left_join(scores, by = "barcode_w_prefix", relationship = "one-to-one") |>
        dplyr::filter(is.finite(.data$WNN_UMAP_1), is.finite(.data$WNN_UMAP_2))
      # Match production random drawing.
      data <- withr::with_seed(1L, dplyr::slice_sample(data, prop = 1))
      enriched <- chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x |>
        dplyr::filter(.data$GWAS_ID == trait, .data$z_q < .05, is_figure_1_cell_type(.data$cluster))
      plot_figure_1_TRS(data, trait, gsub("-", "_", enriched$cluster)) |> compact_figure_1_plot()
    }),
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    figure_1_plot.mixed_human_31x,
    description = "Combine manuscript Figure 1 ggplot panels in four rows: A, BC, D, EF",
    command = compose_figure_1(
      figure_1_panel_A.mixed_human_31x,
      figure_1_panel_B.mixed_human_31x,
      figure_1_panel_C.mixed_human_31x,
      figure_1_panel_D.mixed_human_31x,
      figure_1_panel_E.mixed_human_31x,
      figure_1_panel_F.mixed_human_31x) |>
      compact_figure_1_plot(),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    figure_1.mixed_human_31x,
    description = "Export full-page manuscript Figure 1 from its combined ggplot object",
    command = {
      path <- file.path(output_dir, "1.png")
      ggplot2::ggsave(path, figure_1_plot.mixed_human_31x,
        device = ragg::agg_png, width = 8, height = 9.5, dpi = 300, create.dir = TRUE)
      path
    },
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    figure_1_data.mixed_human_31x,
    description = "Figure 1 values that the manuscript quotes: the aggregation's retained nuclei and assigned cell types, its enriched (q < 0.05) cell type-study pairs, studies and cell types, the selected peak-gene link and, for each heatmap trait's selected cell type, its z, q, contributing loci and the largest locus share of their summed contribution",
    command = {
      paths <- file.path(fs::dir_create(file.path(output_dir, "data")),
        c("figure_1_aggregation.tsv", "figure_1_link.tsv", "figure_1_traits.tsv"))
      cell_types <- metadata_w_cell_types_tibble.WNN.mixed_human_31x$WNN_harmony_SNN_cluster_cell_type
      enriched <- figure_1_enrichments(chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x)
      readr::write_tsv(tibble::tibble(retained_nuclei = length(cell_types),
        assigned_cell_types = dplyr::n_distinct(cell_types[is_figure_1_cell_type(cell_types)]),
        enriched_pairs = nrow(enriched), enriched_studies = dplyr::n_distinct(enriched$GWAS_ID),
        enriched_cell_types = dplyr::n_distinct(enriched$cluster)), paths[[1]])
      select_figure_1_link(peak_gene_correlation_top_links_tibble.WNN.peak_gene_correlation.mixed_human_31x) |>
        dplyr::select("cell_group", "TargetGene", "peak", "distance", "n_donors", "n_aggregates",
          "correlation", "hierarchical_pvalue") |>
        readr::write_tsv(paths[[2]])
      figure_1_traits.mixed_human_31x |>
        dplyr::select("GWAS_ID", "cluster", "z", "z_q", "n_loci", "top_locus_share") |>
        readr::write_tsv(paths[[3]])
      paths
    },
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    supplementary_figure_S3.mixed_human_31x,
    description = "Supplementary Figure S3: the module's unscaled cell-type deviation heatmap for every GWAS study and assigned cell type, with trait-category and loci tracks, on Figure 1D's colour scale with its q < 0.05 stars, and its values",
    command = {
      scores <- chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.mixed_human_31x |>
        dplyr::filter(is_figure_1_cell_type(.data$cluster))
      paths <- file.path(output_dir, c("S3.png", "data/cell_type_trait_enrichment.tsv"))
      ggplot2::ggsave(paths[[1]], plot_figure_S3(scores, GWAS_metadata_tracks_plot.genetic_enrichment.mixed_human_31x,
        figure_1_traits.mixed_human_31x), device = ragg::agg_png, width = 9, height = 10.5, dpi = 300)
      readr::write_tsv(dplyr::select(scores, "GWAS_ID", "Category", "cluster", "deviation", "z", "z_q"), paths[[2]])
      paths
    },
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    supplementary_table_S3.mixed_human_31x,
    description = "Supplementary Table S3: Open Targets metadata and credible sets of each configured GWAS study, with the credible-set loci and variants that overlap consensus peaks and the peaks they weight",
    command = {
      path <- file.path(output_dir, "data", "table_S3.tsv")
      GWAS_study_tibble(GWAS_inputs_tibble.genetic_enrichment.mixed_human_31x,
        GWAS_peak_variant_weight_records.genetic_enrichment.mixed_human_31x, open_targets_study_dataset_path) |>
        readr::write_tsv(path)
      path
    },
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  tarchetypes::tar_file(
    supplementary_table_S6,
    description = "Supplementary Table S6, the data of Figure S2: multiomeR against Seurat/Signac like for like, sizing with defaults by scope and stage, revisions, and the simulated wall time against measured runs, from recorded runs, the targets graph and Slurm job monitoring",
    command = {
      context <- benchmark_context()
      sizing <- benchmark_sizing_tibbles(context)
      tables <- list(seurat_signac_comparison_resources = seurat_signac_comparison_resource_tibble(context),
        benchmark_sizes = benchmark_cellranger_inputs(BENCHMARK_SIZING_AGGREGATIONS),
        benchmark_sizing = sizing$sizing, benchmark_sizing_stages = sizing$stages,
        benchmark_revisions = benchmark_revision_tibble(context),
        benchmark_simulation_check = benchmark_simulation_check(context),
        benchmark_worker_attribution = benchmark_worker_attribution_tibble(context))
      paths <- file.path(fs::dir_create(file.path(output_dir, "data")), paste0(names(tables), ".tsv"))
      purrr::walk2(tables, paths, readr::write_tsv)
      paths
    },
    # targets lets only file targets read the store metadata, which no upstream dependency tracks.
    # A worker rebuilds the graph and reads the monitor history, so dispatching continues meanwhile.
    cue = targets::tar_cue(mode = "always"), resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    supplementary_figure_S2_plot,
    description = "Supplementary Figure S2: like for like, sizing with defaults and revisions from the Table S6 data",
    command = stats::setNames(supplementary_table_S6, tools::file_path_sans_ext(basename(supplementary_table_S6))) |>
      lapply(readr::read_tsv, show_col_types = FALSE) |>
      plot_benchmark_figure() |> compact_figure_1_plot(),
    packages = w_def("patchwork"), resources = get_tar_resources(RAM_GB_req = 1)
  ),
  tarchetypes::tar_file(
    supplementary_figure_S2,
    description = "Export Supplementary Figure S2 from its ggplot object",
    command = {
      path <- file.path(output_dir, "S2.png")
      ggplot2::ggsave(path, supplementary_figure_S2_plot,
        device = ragg::agg_png, width = 8, height = 10.5, dpi = 300, create.dir = TRUE)
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
        LSI_BPCells.ATAC.comparison_5x, "comparison_pbmc_10k_chromium_x"),
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
        # The recorded host names an internal compute node.
        measurements <- dplyr::bind_rows(lapply(runs, `[[`, "measurement")) |> dplyr::select(-"host")
        summary <- summarize_algorithm_parity(measurements, checks)
        dir <- file.path(output_dir, "data", "algorithm_parity")
        fs::dir_create(dir)
        paths <- file.path(dir, c("measurements.tsv", "parity.tsv", "scdblfinder_atac_agreement.tsv", "summary.tsv"))
        readr::write_tsv(measurements, paths[[1]])
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
# set after the other targets are defined, so they apply to this chain only; only
# its peak calling, which runs MACS3 on six clusters at once, needs more RAM.
comparison_resources <- get_tar_resources(cores_req = 6, RAM_GB_req = 60)

targets::tar_option_set(
  packages = c("Seurat", "SeuratObject", "Signac"),
  resources = comparison_resources
)

comparison_aggregation_tibble <- tibble::tibble(
  aggregation = comparison_aggregation_names(),
  comparison_peak_calling_cut_site_dirs = rlang::syms(
    paste0("peak_calling_cut_sites_per_GEM_well.fragments.ATAC.", aggregation)
  ),
  comparison_peak_calling_discovery_tibble = rlang::syms(paste0("peak_calling_cluster_discovery_tibble.ATAC.", aggregation)),
  comparison_blacklist_GRanges = rlang::syms(paste0("blacklist_GRanges.ATAC.", aggregation))
)

comparison_seurat_signac_targets <- function(values, peak_resources) {
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
    tarchetypes::tar_file(
      comparison_seurat_signac_RNA_BPCells_dirs,
      {
        comparison_seurat_signac_cellranger_input_files
        comparison_write_RNA_dirs(comparison_seurat_signac_config, get_structured_file_path())
      }
    ),
    targets::tar_target(
      comparison_seurat_signac_RNA_object,
      comparison_process_RNA(comparison_seurat_signac_RNA_BPCells_dirs, comparison_seurat_signac_config)
    ),
    targets::tar_target(
      comparison_seurat_signac_consensus_peak_GRanges,
      comparison_call_consensus_peaks(
        comparison_peak_calling_cut_site_dirs,
        comparison_peak_calling_discovery_tibble,
        comparison_seurat_signac_config,
        comparison_blacklist_GRanges
      ),
      resources = peak_resources
    ),
    tarchetypes::tar_file(
      comparison_seurat_signac_ATAC_BPCells_dirs,
      comparison_write_ATAC_dirs(
        comparison_seurat_signac_RNA_BPCells_dirs,
        comparison_seurat_signac_config,
        comparison_seurat_signac_consensus_peak_GRanges,
        get_structured_file_path()
      )
    ),
    targets::tar_target(
      comparison_seurat_signac_object,
      comparison_build_object(
        comparison_seurat_signac_RNA_object,
        comparison_seurat_signac_ATAC_BPCells_dirs,
        comparison_seurat_signac_config
      )
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
  "manifests/cfg_pipeline_parameters.tsv"
)

comparison_seurat_signac_pipeline <- rlang::list2(
  targets::tar_target(
    comparison_seurat_signac_config_files,
    comparison_seurat_signac_config_paths,
    format = "file"
  ),
  comparison_seurat_signac_targets(
    values = dplyr::filter(comparison_aggregation_tibble, .data$aggregation %in% c("comparison_1x", "comparison_2x")),
    peak_resources = comparison_resources
  ),
  comparison_seurat_signac_targets(
    values = dplyr::filter(comparison_aggregation_tibble, .data$aggregation == "comparison_5x"),
    peak_resources = get_tar_resources(cores_req = 6, RAM_GB_req = 200)
  ),
  comparison_seurat_signac_targets(
    values = dplyr::filter(comparison_aggregation_tibble,
      .data$aggregation %in% c("comparison_10x", "comparison_20x", "comparison_31x")),
    peak_resources = get_tar_resources(cores_req = 6, RAM_GB_req = 500)
  )
)

c(base_pipeline, figure_pipeline, algorithm_parity_pipeline, comparison_seurat_signac_pipeline)

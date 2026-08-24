rlang::list2(
  ATAC_cluster_specific_fragment_file_generation = rlang::list2(
    targets::tar_target(
      name = blacklist_GRanges.ATAC,
      description = "Load genome blacklist regions for the configured genome",
      command = with_annotation_hub_cache_lock({
        AnnotationHub::setAnnotationHubOption("CACHE", file.path(targets::tar_config_get("store"), "files", "AnnotationHub"))
        get_blacklist_GRanges(
          genome = aggregated_cellranger_ref_list$genomes[[1]],
          annot_hub_interface = AnnotationHub::AnnotationHub(ask = FALSE)
        )
      }),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = combined_BPCells_fragment_obj.ATAC,
      description = "Merge per GEM well BPCells fragment objects and filter to standard chromosomes [part_of_graph:ATAC] [part_of_graph:parallel] [part_of_graph:seurat_export]",
      command = {
        standard_chroms <- switch(
          organism_chr,
          "Homo_sapiens" = paste0("chr", c(1:22, "X", "Y")),
          "Mus_musculus" = paste0("chr", c(1:19, "X", "Y")),
          stop("Invalid organism.")
        )

        combined_fragments <- purrr::reduce(aggregation_fragments_w_prefix_bpcells_syms, c)
        BPCells::select_chromosomes(combined_fragments, standard_chroms)
      },
      packages = w_def("BPCells")
    ),
    targets::tar_target(
      name = BCs_per_peak_cluster_list.ATAC,
      description = "Split cell barcodes into per-cluster lists for peak calling",
      command = {
        cluster_metadata <- group_split_by(metadata_w_cell_types_tibble.GEX, aggregation_call_peaks_by_cluster_col)
        purrr::map(cluster_metadata, \(cluster_df) dplyr::pull(cluster_df, barcode_w_prefix))
      }
    ),
    targets::tar_target(
      name = peak_calling_cluster_names.ATAC,
      description = "Get sorted unique cluster names used for peak calling",
      command = {
        cluster_names <- as.character(metadata_w_cell_types_tibble.GEX[[aggregation_call_peaks_by_cluster_col]])
        gtools::mixedsort(unique(cluster_names))
      }
    ),
    targets::tar_target(
      name = peak_calling_cluster_discovery_tibble.ATAC,
      description = "Select full-genome ATAC peak-calling discovery cells per cluster with optional downsampling",
      command = build_peak_calling_cluster_discovery_tibble(
        BCs_per_peak_cluster_list = BCs_per_peak_cluster_list.ATAC,
        peak_calling_cluster_names = peak_calling_cluster_names.ATAC
      ),
      iteration = "vector"
    ),
    tarchetypes::tar_file(
      name = fragments_per_cluster.fragments.ATAC,
      description = "Export full-genome ATAC fragments for each peak-calling cluster",
      command = write_ATAC_fragments_for_peak_calling_cluster(
        ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
        BCs_per_peak_cluster = peak_calling_cluster_discovery_tibble.ATAC$BCs_per_peak_cluster[[1]],
        peak_calling_cluster_name = peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]]
      ),
      pattern = map(peak_calling_cluster_discovery_tibble.ATAC)
    ),
    tarchetypes::tar_file(
      name = fragments_per_peak_calling_cluster_discovery.fragments.ATAC,
      description = "Export full-genome ATAC fragments used for per-cluster peak discovery",
      command = write_ATAC_fragments_for_peak_calling_cluster(
        ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
        BCs_per_peak_cluster = peak_calling_cluster_discovery_tibble.ATAC$BCs_for_peak_discovery[[1]],
        peak_calling_cluster_name = peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]],
        output_suffix = paste0(peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]], "__discovery")
      ),
      pattern = map(peak_calling_cluster_discovery_tibble.ATAC)
    )
  ),
  ATAC_peak_calling_per_cluster_w_iterative_collapsing_targets = rlang::list2(
    tarchetypes::tar_file(
      name = peaks_per_cluster_narrowPeaks.peaks.ATAC,
      description = "Call full-genome ATAC peaks per cluster using capped discovery cells [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = {
        peak_calling_cluster_name <- peak_calling_cluster_discovery_tibble.ATAC$peak_calling_cluster_name[[1]]
        if (identical(aggregation_ATAC_peak_calling_method, "macs3")) {
          call_peaks_w_MACS3(
            ATAC_fragments_per_cluster = fragments_per_peak_calling_cluster_discovery.fragments.ATAC,
            ATAC_peak_calling_cluster_names = peak_calling_cluster_name,
            genome = aggregated_cellranger_ref_list$genomes[[1]],
            output_suffix = peak_calling_cluster_name,
            allow_no_peaks = TRUE
          )
        } else if (identical(aggregation_ATAC_peak_calling_method, "bpcells_tile")) {
          call_peaks_w_BPCells_tile(
            ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
            ATAC_BCs_per_peak_cluster = peak_calling_cluster_discovery_tibble.ATAC$BCs_for_peak_discovery[[1]],
            ATAC_peak_calling_cluster_names = peak_calling_cluster_name,
            genome = aggregated_cellranger_ref_list$genomes[[1]],
            output_suffix = peak_calling_cluster_name,
            allow_no_peaks = TRUE
          )
        } else {
          stop("aggregation_ATAC_peak_calling_method must be either 'macs3' or 'bpcells_tile'.")
        }
      },
      pattern = map(fragments_per_peak_calling_cluster_discovery.fragments.ATAC, peak_calling_cluster_discovery_tibble.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = peak_GRanges_per_cluster.ATAC,
      description = "Parse full-genome per-cluster ATAC peak calls as fixed-width blacklist-filtered GRanges",
      command = get_peak_GRanges_w_fixed_width(
        peaks_per_cluster_narrowPeaks.peaks.ATAC,
        genome = aggregated_cellranger_ref_list$genomes[[1]],
        blacklist_GRanges = blacklist_GRanges.ATAC
      ),
      pattern = map(peaks_per_cluster_narrowPeaks.peaks.ATAC)
    ),
    targets::tar_target(
      name = chromHMMs_list_proj.ATAC,
      description = "Subset chromHMM annotation list to the configured Roadmap/EDACC states",
      command = chromHMMs_list_general[aggregation_roadmap_EDACC_names]
    ),
    targets::tar_target(
      name = within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC,
      description = "Collapse per-cluster ATAC peaks into non-overlapping peak sets [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = {
        collapsed_peaks <- combine_collapse_GRanges_ArchR(list(peak_GRanges_per_cluster.ATAC), by = "neg_log10pvalue_summit")
        format_peak_GRanges(collapsed_peaks, cluster_name = peak_calling_cluster_names.ATAC)
      },
      pattern = map(peak_GRanges_per_cluster.ATAC, peak_calling_cluster_names.ATAC)
    ),
    tarchetypes::tar_file(
      name = peaks_similarity_tiles_plot.ATAC,
      description = "Plot pairwise peak-set similarity matrix across clusters. [checkpoint:ATAC]",
      command = {
        plot <- plot_similarity_matrix_from_GRanges_list(within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC)
        save_plots_structured(plot)
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    )
  ),
  ATAC_between_clusters_iterative_peak_collapsing_targets = rlang::list2(
    targets::tar_target(
      name = consensus_peak_GRanges.ATAC,
      description = "Iteratively collapse per-cluster peak sets into a single consensus peak set [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = combine_collapse_GRanges_ArchR(GRanges_list = within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC, by = "fold_change")
    ),
    targets::tar_target(
      name = consensus_peak_annotated_GRanges.ATAC,
      description = "Annotate the final consensus peak set with transcript, gene, and chromHMM metadata",
      command = {
        peak_GRanges <- CAGEfightR::assignTxType(
          consensus_peak_GRanges.ATAC,
          marker_validated_Ensembl_annotations_GRanges_list$transcripts
        ) |>
          CAGEfightR::assignGeneID(
            marker_validated_Ensembl_annotations_GRanges_list$genes,
            outputColumn = "geneName"
          )

        if (length(chromHMMs_list_proj.ATAC) > 0) {
          peak_GRanges <- purrr::reduce2(
            .x = chromHMMs_list_proj.ATAC,
            .y = names(chromHMMs_list_proj.ATAC),
            .init = peak_GRanges,
            .f = \(accumulated_GRanges, chromHMM_txModel, chromHMM_name) {
              CAGEfightR::assignTxType(
                object = accumulated_GRanges,
                txModels = chromHMM_txModel,
                outputColumn = chromHMM_name
              )
            }
          )
        }

        peak_GRanges
      }
    ),
    # tar_target(
    #   within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC |>
    #     map(\(x) mutate(as_tibble(as.data.frame(x)), chr = seqnames)) |>
    #     bind_rows() |>
    #     arrange(desc(fold_change)) |>
    #     BPCells::merge_peaks_iterative() # TODO: doesn't give the same result at all? 150 k peaks for muscle test, vs 250 k for traditional approach. Bug in BPCells?
    # ),
    tarchetypes::tar_file(
      name = consensus_peak_BPCells_matrix_dir.ATAC,
      description = "Compute the consensus peak-by-cell count matrix using BPCells [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = {
        out_dir <- get_structured_file_path()
        frags <- BPCells::select_cells(
          combined_BPCells_fragment_obj.ATAC,
          metadata_w_cell_types_tibble.GEX$barcode_w_prefix
        )

        peaks_df <- GenomicRanges::as.data.frame(consensus_peak_GRanges.ATAC) |>
          dplyr::transmute(chr = as.character(seqnames), start, end)

        consensus_mat <- BPCells::peak_matrix(frags, peaks_df, mode = aggregation_ATAC_peak_matrix_mode)
        rownames(consensus_mat) <- names(consensus_peak_GRanges.ATAC)
        BPCells::write_matrix_dir(consensus_mat, out_dir, overwrite = TRUE)
        out_dir
      },
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = consensus_peak_BPCells_matrix.ATAC,
      description = "Open the consensus ATAC peak BPCells matrix",
      command = BPCells::open_matrix_dir(consensus_peak_BPCells_matrix_dir.ATAC),
      resources = get_tar_resources(RAM_GB_req = 8)
    )
  ),
  ATAC_peak_diff_accesibility_testing_targets = rlang::list2(
    targets::tar_target(
      name = diff_peak_acc_results_tibble.ATAC,
      description = "Test differential peak accessibility across clusters with a Wilcoxon test",
      command = {
        groups <- metadata_w_cell_types_tibble.GEX[[aggregation_call_peaks_by_cluster_col]] |>
          as.character() |>
          purrr::set_names(metadata_w_cell_types_tibble.GEX$barcode_w_prefix)
        groups <- groups[colnames(consensus_peak_BPCells_matrix.ATAC)]

        BPCells::marker_features(consensus_peak_BPCells_matrix.ATAC, groups, method = "wilcoxon") |>
          dplyr::mutate(
            avg_log2FC = log2(foreground_mean / background_mean),
            p_val_adj = stats::p.adjust(p_val_raw, method = "BH"),
            cluster = foreground,
            gene = feature
          ) |>
          dplyr::mutate(peak = gene) |>
          tidyr::separate_wider_delim(cols = peak, names = c("chr", "start", "end"), delim = "-")
      },
      packages = w_def("BPCells"),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    tarchetypes::tar_file(
      name = diff_acc_peaks_per_cluster_BEDs.peaks.ATAC,
      description = "Export significant differentially accessible peaks per cluster as BED files",
      command = {
        BED_out_path <- get_structured_file_path(
          filetype = "BED",
          override_suffix = peak_calling_cluster_names.ATAC
        )

        cluster_peaks <- dplyr::filter(
          diff_peak_acc_results_tibble.ATAC,
          cluster == peak_calling_cluster_names.ATAC &
            avg_log2FC > 0.5 &
            p_val_adj < 0.05 &
            foreground_mean >= 0.10
        )
        cluster_peak_GRanges <- GenomicRanges::makeGRangesFromDataFrame(cluster_peaks, keep.extra.columns = TRUE)
        rtracklayer::export(cluster_peak_GRanges, con = BED_out_path, format = "bed")

        BED_out_path
      },
      pattern = map(peak_calling_cluster_names.ATAC)
    )
  ),

  ATAC_peak_objects_targets = rlang::list2(
    targets::tar_target(
      name = consensus_peak_tibble.ATAC,
      description = "Annotated tibble representation of the consensus peak set",
      command = tibble::as_tibble(as.data.frame(consensus_peak_annotated_GRanges.ATAC))
    )
  ),

  ATAC_coverage_plots_targets = rlang::list2(
    targets::tar_target(
      name = coverage_regions_tibble.ATAC,
      description = "Build genomic regions to display in ATAC coverage tracks",
      command = get_coverage_regions_tibble_BPCells(
        formatted_peak_tibble = consensus_peak_tibble.ATAC,
        gene_GRanges = marker_validated_Ensembl_annotations_GRanges_list$genes,
        marker_genes = ATAC_coverage_marker_genes
      ),
      iteration = "vector"
    ),
    tarchetypes::tar_file(
      name = coverage_tracks_plots.ATAC,
      description = "Plot ATAC coverage tracks at marker gene and DA peak loci. [checkpoint:ATAC]",
      command = plot_coverage_at_region_BPCells(
        fragments = combined_BPCells_fragment_obj.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        region_id = coverage_regions_tibble.ATAC$region,
        collapsed_peak_tibble = consensus_peak_tibble.ATAC,
        gene_GRanges = marker_validated_Ensembl_annotations_GRanges_list$genes
      ) |>
        save_plots_structured(override_suffix = coverage_regions_tibble.ATAC$name, dyn_suffix_in_subdir = TRUE),
      pattern = map(coverage_regions_tibble.ATAC)
    )
  ),

  ATAC_dim_reduc_and_plots_targets = rlang::list2(
    targets::tar_target(
      name = metadata_w_QC_tibble.ATAC,
      description = "Compute ATAC QC metadata from BPCells peak and fragment objects",
      command = get_ATAC_QC_metadata_from_BPCells(
        metadata_df = metadata_w_cell_types_tibble.GEX,
        ATAC_peak_BPCells_matrix = consensus_peak_BPCells_matrix.ATAC,
        ATAC_combined_BPCells_fragment_obj = combined_BPCells_fragment_obj.ATAC,
        blacklist_GRanges = blacklist_GRanges.ATAC,
        ATAC_peak_GRanges = consensus_peak_GRanges.ATAC,
        genome = aggregated_cellranger_ref_list$genomes[[1]],
        peak_matrix_mode = aggregation_ATAC_peak_matrix_mode
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = QC_excluded_BCs_list.ATAC,
      description = "Identify barcodes failing ATAC QC thresholds",
      command = get_excluded_BCs(metadata_w_QC_tibble.ATAC, QC_exclude_vector = aggregation_QC_exclude_list_combined_object)
    ),
    targets::tar_target(
      name = QC_filtered_BCs.ATAC,
      description = "Barcodes passing BPCells-native ATAC QC filters",
      command = setdiff(metadata_w_QC_tibble.ATAC$barcode_w_prefix, unlist(QC_excluded_BCs_list.ATAC))
    ),
    targets::tar_target(
      name = metadata_filtered_tibble.ATAC,
      description = "BPCells-native ATAC metadata after QC filtering",
      command = dplyr::filter(metadata_w_QC_tibble.ATAC, .data$barcode_w_prefix %in% QC_filtered_BCs.ATAC)
    ),
    targets::tar_target(
      name = peak_QC_filtered_BPCells_matrix.ATAC,
      description = "Lazy QC-filtered ATAC peak BPCells matrix",
      command = consensus_peak_BPCells_matrix.ATAC[, QC_filtered_BCs.ATAC, drop = FALSE],
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    targets::tar_target(
      name = LSI_BPCells.ATAC,
      description = "Run BPCells-native TF-IDF and SVD on the QC-filtered ATAC peak matrix [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = run_ATAC_LSI_BPCells(
        ATAC_peak_BPCells_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
        n_components = utils::tail(aggregation_ATAC_data_PCs, 1),
        threads = 6
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = LSI_embeddings_tibble.ATAC,
      description = "Cell embeddings from BPCells-native ATAC LSI",
      command = LSI_BPCells.ATAC$cell_embeddings |>
        tibble::as_tibble(rownames = "barcode_w_prefix")
    ),
    targets::tar_target(
      name = LSI_loadings_tibble.ATAC,
      description = "Peak loadings from BPCells-native ATAC LSI",
      command = LSI_BPCells.ATAC$feature_loadings |>
        tibble::as_tibble(rownames = "peak")
    ),
    targets::tar_target(
      name = harmony_embeddings_matrix.ATAC,
      description = "Harmony-corrected BPCells-native ATAC LSI embeddings [part_of_graph:ATAC] [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = run_harmony_on_embedding_matrix(
        embedding_matrix = LSI_BPCells.ATAC$cell_embeddings,
        metadata_tibble = metadata_filtered_tibble.ATAC,
        harmony_correction_metadata_col_names = c(aggregation_harmony_correction_metadata_col_names, aggregation_extra_harmony_covars_ATAC),
        dims = aggregation_ATAC_data_PCs,
        cores = 6
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = UMAP_embeddings_tibble.ATAC,
      description = "UMAP coordinates from BPCells-native Harmony-corrected ATAC LSI embeddings [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = run_UMAP_from_embedding_matrix(
        embedding_matrix = harmony_embeddings_matrix.ATAC,
        dims = aggregation_UMAP_ATAC_PCs,
        n_neighbors = aggregation_UMAP_nNNs,
        min_dist = aggregation_UMAP_min_dist
      ) |>
        tibble::as_tibble(rownames = "barcode_w_prefix"),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = LSI_clusters.ATAC,
      description = "Leiden clusters from a BPCells nearest-neighbour graph on Harmony-corrected ATAC LSI embeddings [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = cluster_embedding_matrix_BPCells(
        embedding_matrix = harmony_embeddings_matrix.ATAC,
        dims = aggregation_ATAC_data_PCs,
        k = aggregation_data_nNNs,
        resolution = aggregation_ATAC_cluster_res,
        threads = 6,
        min_barcodes = aggregation_cluster_min_barcodes
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = metadata_w_clusters_tibble.ATAC,
      description = "BPCells-native ATAC metadata with LSI clusters and UMAP coordinates",
      command = metadata_filtered_tibble.ATAC |>
        dplyr::filter(.data$barcode_w_prefix %in% names(LSI_clusters.ATAC)) |>
        dplyr::left_join(UMAP_embeddings_tibble.ATAC, by = "barcode_w_prefix") |>
        dplyr::mutate(LSI_harmony_SNN_cluster = LSI_clusters.ATAC[.data$barcode_w_prefix])
    ),
    targets::tar_target(
      name = metadata_w_cell_types_unfiltered_tibble.ATAC,
      description = "Assign cell type labels to BPCells-native ATAC clusters from GEX UCell marker scores in metadata",
      command = add_cell_types_to_metadata_from_module_scores(
        metadata_tibble = metadata_w_clusters_tibble.ATAC,
        named_marker_genes_list = UCell_GEX_marker_genes_list,
        allow_multiple_cell_types = aggregation_allow_multiple_cell_types,
        cluster_column = "LSI_harmony_SNN_cluster"
      )
    ),
    tarchetypes::tar_file(
      name = VizDimLoadings_plots.ATAC,
      description = "Plot top feature loadings for each LSI dimension. [checkpoint:ATAC]",
      command = plot_LSI_loadings_from_tibble(
        LSI_loadings_tibble = LSI_loadings_tibble.ATAC,
        dims = aggregation_ATAC_data_PCs,
        nfeatures = 50
      ) |>
        save_plots_structured()
    ),
    tarchetypes::tar_file(
      name = LSI_singular_values_elbow_plot.ATAC,
      description = "Elbow plot of native ATAC LSI singular values. [checkpoint:ATAC]",
      command = {
        plot <- plot_embedding_singular_values(
          singular_values = LSI_BPCells.ATAC$singular_values,
          dims = aggregation_ATAC_data_PCs
        )
        save_plots_structured(plot)
      }
    ),
    tarchetypes::tar_file(
      name = LSI_embedding_sdev_plot.ATAC,
      description = "Non-Harmony and Harmony ATAC LSI embedding coordinate SD plot. [checkpoint:ATAC]",
      command = {
        plot <- plot_embedding_sdev(
          embedding_matrix = LSI_BPCells.ATAC$cell_embeddings,
          dims = aggregation_ATAC_data_PCs,
          harmony_embedding_matrix = if (
            length(c(aggregation_harmony_correction_metadata_col_names, aggregation_extra_harmony_covars_ATAC) %||% character()) > 0
          ) {
            harmony_embeddings_matrix.ATAC
          } else {
            NULL
          },
          dim_prefix = "LSI_"
        )
        save_plots_structured(plot)
      }
    ),
    tarchetypes::tar_file(
      name = LSI_metadata_association_barplots.ATAC,
      description = "Non-Harmony and Harmony ATAC LSI metadata association bar plots. [checkpoint:ATAC]",
      command = {
        plots <- plot_embedding_metadata_association_barplots(
          embedding_matrix = LSI_BPCells.ATAC$cell_embeddings,
          harmony_embedding_matrix = if (
            length(c(aggregation_harmony_correction_metadata_col_names, aggregation_extra_harmony_covars_ATAC) %||% character()) > 0
          ) {
            harmony_embeddings_matrix.ATAC
          } else {
            NULL
          },
          metadata_tibble = metadata_filtered_tibble.ATAC,
          dims = aggregation_ATAC_data_PCs,
          dim_prefix = "LSI_",
          continuous_technical_cols = c(
            "log10_nCount_ATAC",
            "atac_fragments",
            "atac_peak_counts_frac",
            "atac_peak_counts_blacklist_frac",
            "TSS.enrichment",
            "nucleosome_signal",
            "amulet_q.value",
            "log10_nCount_RNA",
            "vireo_max_prob_singlet"
          ),
          categorical_technical_cols = c(
            "GEM_well_ID",
            "multiplex_batch",
            "multiplex_pool",
            "run_harmony",
            "run_harmony_batched",
            "vireo_type",
            "discarded_by_cellranger",
            "discarded_by_amulet"
          ),
          continuous_biological_cols = aggregation_continuous_vars %||% character(),
          categorical_biological_cols = aggregation_categorical_vars %||% character()
        )
        save_plots_structured(
          plots,
          height = purrr::map_dbl(plots, \(plot) max(5, 0.75 * get_num_facet_rows(plot)))
        )
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    )
  ),

  ATAC_graph_and_cluster_targets = rlang::list2(
    targets::tar_target(
      name = scDblFinder_GEM_well_tibble.ATAC,
      description = "Prepare per GEM well ATAC barcodes and cluster labels for scDblFinder",
      command = prepare_scDblFinder_GEM_well_tibble(
        metadata_tibble = metadata_w_cell_types_unfiltered_tibble.ATAC |>
          dplyr::filter(.data$barcode_w_prefix %in% QC_filtered_BCs.ATAC),
        cluster_collapse_list = aggregation_scDblFinder_GEX_cell_type_collapse_list,
        cluster_col = "LSI_harmony_SNN_cluster_cell_type"
      ),
      iteration = "vector"
    ),
    targets::tar_target(
      name = scDblFinder_feature_groups.ATAC,
      description = "Cluster consensus ATAC peaks from global LSI loadings for memory-light scDblFinder input",
      command = get_feature_groups_from_LSI_loadings(
        LSI_loadings_tibble = LSI_loadings_tibble.ATAC,
        dims = intersect(2:20, aggregation_ATAC_data_PCs),
        n_groups = 50,
        seed = 1
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = scDblFinder_feature_matrix.ATAC,
      description = "Aggregate ATAC peak counts into global LSI-derived feature groups for scDblFinder",
      command = aggregate_BPCells_rows_by_group(
        feature_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
        feature_groups = scDblFinder_feature_groups.ATAC,
        threads = 6
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = scDblFinder_results_by_GEM_well_tibble.ATAC,
      description = "Run ATAC scDblFinder independently for each 10x Genomics GEM well from pre-aggregated peak-group counts",
      command = run_scDblFinder_BPCells_GEM_well(
        feature_matrix = scDblFinder_feature_matrix.ATAC,
        scDblFinder_GEM_well_tibble = scDblFinder_GEM_well_tibble.ATAC,
        output_suffix = "ATAC",
        dbr.sd = 1.0,
        aggregateFeatures = FALSE,
        nfeatures = 50,
        processing = "normFeatures"
      ),
      pattern = map(scDblFinder_GEM_well_tibble.ATAC),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    # Previous full-peak ATAC scDblFinder path. Kept as fallback, but disabled because
    # scDblFinder::aggregateFeatures() recomputes TF-IDF and can OOM before aggregation.
    # tar_target(
    #   name = scDblFinder_results_by_GEM_well_tibble.ATAC,
    #   description = "Run ATAC scDblFinder independently for each 10x Genomics GEM well from BPCells peak-count slices",
    #   command = run_scDblFinder_BPCells_GEM_well(
    #     feature_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
    #     scDblFinder_GEM_well_tibble = scDblFinder_GEM_well_tibble.ATAC,
    #     output_suffix = "ATAC",
    #     dbr.sd = 1.0,
    #     aggregateFeatures = TRUE,
    #     nfeatures = 50,
    #     processing = "normFeatures"
    #   ),
    #   pattern = map(scDblFinder_GEM_well_tibble.ATAC),
    #   resources = get_tar_resources(RAM_GB_req = 16)
    # ),
    targets::tar_target(
      name = scDblFinder_results_df.ATAC,
      description = "Combine per GEM well ATAC scDblFinder classifications [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = scDblFinder_results_by_GEM_well_tibble.ATAC |>
        dplyr::bind_rows() |>
        dplyr::select(-dplyr::any_of("GEM_well_ID")),
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    tarchetypes::tar_file(
      name = scDblFinder_score_violins_plot.ATAC,
      description = "Violin plots of scDblFinder doublet scores by cluster. [checkpoint:ATAC]",
      command = {
        scDblFinder_metadata <- metadata_w_cell_types_unfiltered_tibble.ATAC |>
          dplyr::left_join(scDblFinder_results_df.ATAC, by = "barcode_w_prefix")
        plot_data <- scDblFinder_metadata |>
          dplyr::select(LSI_harmony_SNN_cluster, LSI_harmony_SNN_cluster_cell_type, scDblFinder.score_ATAC, scDblFinder.class_ATAC, GEM_well_ID) |>
          tidyr::pivot_longer(
            cols = dplyr::all_of(c("LSI_harmony_SNN_cluster", "LSI_harmony_SNN_cluster_cell_type")),
            names_to = "cluster_type",
            values_to = "cluster_id"
          )

        plot <- plot_data |>
          ggplot2::ggplot(ggplot2::aes(x = cluster_id, y = scDblFinder.score_ATAC)) +
          ggplot2::geom_violin(scale = "width") +
          ggplot2::geom_jitter(
            data = \(d) dplyr::filter(d, scDblFinder.class_ATAC == "doublet"),
            ggplot2::aes(color = GEM_well_ID),
            size = 0.3,
            alpha = 0.5,
            width = 0.2
          ) +
          ggplot2::facet_wrap(~cluster_type, scales = "free_x") +
          ggplot2::labs(subtitle = "Points are cells classified as doublet by scDblFinder, colored by GEM well.") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::theme(legend.position = "none")

        save_plots_structured(plot)
      }
    ),
    targets::tar_target(
      name = metadata_w_cell_types_tibble.ATAC,
      description = "Annotate ATAC metadata with scDblFinder QC and apply configured cell- and cluster-level doublet filters [part_of_graph:ATAC] [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = filter_metadata_by_scDblFinder(
        metadata_tibble = metadata_w_cell_types_unfiltered_tibble.ATAC,
        scDblFinder_results_df = scDblFinder_results_df.ATAC,
        class_col = "scDblFinder.class_ATAC",
        cluster_col = "LSI_harmony_SNN_cluster",
        remove_called_doublets = aggregation_scDblFinder_ATAC_remove_called_doublets,
        max_doublet_fraction_per_cluster = aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster
      )
    )
  ),

  ATAC_QC_plots_targets = rlang::list2(
    tarchetypes::tar_file(
      name = peaks_QC_violins_plot.ATAC,
      description = "Violin plots of peak-based ATAC QC metrics per GEM well. [checkpoint:ATAC]",
      command = {
        plot_data <- metadata_w_QC_tibble.ATAC |>
          dplyr::select(GEM_well_ID, "dataset", dplyr::any_of(aggregation_peak_based_continuous_QC_vars)) |>
          tidyr::pivot_longer(cols = dplyr::any_of(aggregation_peak_based_continuous_QC_vars), names_to = "feature", values_to = "value")

        plot <- plot_data |>
          ggplot2::ggplot(ggplot2::aes(x = GEM_well_ID, y = value, fill = dataset)) +
          ggplot2::geom_violin(scale = "width") +
          ggplot2::facet_wrap(~feature, scales = "free", ncol = 1) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(x = ggplot2::element_blank())
        save_plots_structured(plot)
      }
    ),
    tarchetypes::tar_file(
      name = categorical_bars_plots.ATAC,
      description = "Bar plots of categorical metadata proportions by cell-type cluster. [checkpoint:ATAC]",
      command = {
        plot <- plot_categorical_bars_plot(
          metadata_w_cell_types_tibble.ATAC,
          metadata_cols = aggregation_ATAC_categorical_vars,
          cluster_col = "LSI_harmony_SNN_cluster_cell_type"
        )
        save_plots_structured(plot)
      }
    ),
    tarchetypes::tar_file(
      name = QC_excluded_upset_plot.ATAC,
      description = "UpSet plot of overlapping ATAC QC exclusion reasons. [checkpoint:ATAC]",
      command = {
        plot <- plot_upset_from_excluded_BCs_list(QC_excluded_BCs_list.ATAC, n_total = nrow(metadata_w_QC_tibble.ATAC))
        save_plots_structured(plot)
      }
    ),
    targets::tar_target(
      name = categorical_UMAP_var.ATAC,
      description = "Categorical ATAC UMAP variables to plot one at a time",
      command = aggregation_ATAC_categorical_vars,
      iteration = "vector"
    ),
    tarchetypes::tar_file(
      name = categorical.UMAPs.ATAC,
      description = "UMAPs colored by categorical metadata variables. [checkpoint:ATAC]",
      command = metadata_w_cell_types_tibble.ATAC |>
        plot_UMAP_from_metadata(variable = categorical_UMAP_var.ATAC) |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = stringr::str_replace_all(categorical_UMAP_var.ATAC, "[/\\\\]", "_")
        ),
      pattern = map(categorical_UMAP_var.ATAC)
    ),
    targets::tar_target(
      name = UMAP_n_dims_seq.ATAC,
      description = "Sequence of LSI dimension counts for cross-parameter UMAP sweep",
      command = round(seq(5, length(aggregation_ATAC_data_PCs), length.out = 3))
    ),
    tarchetypes::tar_file(
      name = cross.UMAPs.ATAC,
      description = "Compute ATAC UMAPs across a sweep of LSI dimension counts and neighbour counts. [checkpoint:ATAC]",
      command = {
        sweep_umap <- run_UMAP_from_embedding_matrix(
          embedding_matrix = harmony_embeddings_matrix.ATAC[metadata_w_cell_types_tibble.ATAC$barcode_w_prefix, , drop = FALSE],
          dims = 1:UMAP_n_dims_seq.ATAC,
          n_neighbors = UMAP_neighbors_seq,
          min_dist = aggregation_UMAP_min_dist
        ) |>
          tibble::as_tibble(rownames = "barcode_w_prefix")

        plot <- metadata_w_cell_types_tibble.ATAC |>
          dplyr::select(-dplyr::any_of(c("LSI_UMAP_1", "LSI_UMAP_2"))) |>
          dplyr::left_join(sweep_umap, by = "barcode_w_prefix") |>
          plot_UMAP_from_metadata(variable = "LSI_harmony_SNN_cluster_cell_type")

        save_plots_structured(plot, dyn_suffix_in_subdir = TRUE, override_suffix = paste0(UMAP_n_dims_seq.ATAC, "_", UMAP_neighbors_seq))
      },
      pattern = cross(UMAP_n_dims_seq.ATAC, UMAP_neighbors_seq),
      resources = get_tar_resources(RAM_GB_req = 60)
    )
  ),

  ATAC_w_motif_family_accessibility_targets = rlang::list2(
    targets::tar_target(
      name = ATAC_marker_motif_families_list,
      description = "Resolve configured ATAC marker TFs to JASPAR2026 motif-similarity families [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = resolve_marker_motif_families(
        marker_TFs_list = aggregation_ATAC_marker_TFs,
        motif_family_members_tibble = JASPAR_motif_family_members_tibble
      )
    ),
    targets::tar_target(
      name = ATAC_marker_motif_families_vec,
      description = "Flatten configured ATAC marker motif families to plain feature names",
      command = {
        ATAC_marker_motif_families_list |>
          unlist(use.names = FALSE) |>
          stringr::str_remove_all("[+-]$") |>
          unique()
      }
    ),
    targets::tar_target(
      name = peak_TF_motif_family_matrix.ATAC,
      description = "Match JASPAR2026 familial root motifs directly to ATAC peaks [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:differential_analyses]",
      command = get_motif_matrix_from_ATAC_peak_names(
        ATAC_peak_names = rownames(consensus_peak_BPCells_matrix.ATAC),
        ATAC_peak_GRanges = consensus_peak_GRanges.ATAC,
        motif_matrix_list = JASPAR_familial_root_motif_matrix_list,
        genome = aggregated_cellranger_ref_list$genomes[[1]]
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = chromVAR_obj.ATAC,
      description = "Build the lightweight BPCells-backed chromVAR RSE for motif-family accessibility and genetic enrichment [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:genetic_enrichment_single_nucleus]",
      command = get_chromVAR_obj_from_peak_matrix(
        ATAC_peak_matrix = consensus_peak_BPCells_matrix.ATAC[, metadata_w_cell_types_tibble.ATAC$barcode_w_prefix],
        ATAC_peak_GRanges = consensus_peak_GRanges.ATAC,
        genome = aggregated_cellranger_ref_list$genomes[[1]]
      ),
      resources = get_tar_resources(RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = chromVAR_peak_expectation.ATAC,
      description = "Compute shared peak-level accessibility expectations for chromVAR analyses [part_of_graph:genetic_enrichment_single_nucleus]",
      command = get_chromVAR_peak_expectation(chromVAR_obj.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = chromVAR_background_bins.ATAC,
      description = "Compute shared peak-level betterChromVAR background bins for chromVAR analyses [part_of_graph:genetic_enrichment_single_nucleus]",
      command = get_chromVAR_background_bins(
        chromVAR_obj = chromVAR_obj.ATAC,
        expectation = chromVAR_peak_expectation.ATAC
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = chromVAR_cell_chunks.ATAC,
      description = "Split ATAC cells into reusable chromVAR chunks bounded by nonzero counts",
      command = get_chromVAR_cell_chunk_tibble(
        chromVAR_obj = chromVAR_obj.ATAC,
        chunk_nonzero_limit = 2^27
      ) |>
        dplyr::group_by(chunk_id) |>
        targets::tar_group(),
      iteration = "group",
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = chromVAR_chunk_context_records.ATAC,
      description = "Materialize one ATAC cell chunk and compute its reusable betterChromVAR background",
      command = get_chromVAR_chunk_context_record(
        chunk_record = chromVAR_cell_chunks.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC,
        background_bins = chromVAR_background_bins.ATAC,
        expectation = chromVAR_peak_expectation.ATAC
      ),
      pattern = map(chromVAR_cell_chunks.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = chromVAR_TF_motif_family_matrix.ATAC,
      description = "Align JASPAR2026 familial root-motif annotations to the filtered chromVAR peak order [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:differential_analyses]",
      command = align_chromVAR_annotations_to_obj(
        annotations = peak_TF_motif_family_matrix.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC
      ),
      resources = get_tar_resources(RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = motif_family_chromVAR_chunk_results.ATAC,
      description = "Compute per-cell motif-family accessibility scores for one reusable ATAC chunk",
      command = {
        compute_chromVAR_annotation_chunk_result(
          annotations = chromVAR_TF_motif_family_matrix.ATAC,
          chromVAR_obj = chromVAR_obj.ATAC,
          chunk_context_record = chromVAR_chunk_context_records.ATAC,
          compute = c("deviations", "z")
        )
      },
      pattern = map(chromVAR_chunk_context_records.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = motif_family_chromVAR_results.ATAC,
      description = "Combine chunk-level motif-family accessibility scores into the standard chromVAR result object",
      command = combine_chromVAR_chunk_results(
        chunk_results = motif_family_chromVAR_chunk_results.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC,
        annotations = chromVAR_TF_motif_family_matrix.ATAC
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    tarchetypes::tar_file(
      name = motif_family_accessibility_BPCells_matrix_dir.ATAC,
      description = "Write betterChromVAR motif-family accessibility Z-scores to a BPCells matrix directory [part_of_graph:ATAC] [part_of_graph:seurat_export]",
      command = {
        out_dir <- get_structured_file_path()
        column_major_dir <- tempfile(pattern = "motif_family_accessibility_column_major_", tmpdir = dirname(out_dir))
        on.exit(if (fs::dir_exists(column_major_dir)) fs::dir_delete(column_major_dir), add = TRUE)
        if (fs::dir_exists(out_dir)) {
          fs::dir_delete(out_dir)
        }

        motif_family_chromVAR_results.ATAC$chromVAR_z_scores |>
          Matrix::Matrix(sparse = TRUE) |>
          methods::as("dgCMatrix") |>
          BPCells::write_matrix_dir(column_major_dir, compress = FALSE, overwrite = TRUE)
        BPCells::transpose_storage_order(BPCells::open_matrix_dir(column_major_dir), outdir = out_dir)
        out_dir
      },
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = motif_family_accessibility_BPCells_matrix.ATAC,
      description = "Open the disk-backed betterChromVAR motif-family accessibility Z-score matrix",
      command = BPCells::open_matrix_dir(motif_family_accessibility_BPCells_matrix_dir.ATAC),
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    targets::tar_target(
      name = motif_family_accessibility_markers.ATAC,
      description = "Find marker motif-family accessibility scores per ATAC cell type",
      command = get_marker_motif_family_accessibility_from_chromVAR_BPCells_z_scores(
        chromVAR_z_scores_BPCells_matrix = motif_family_accessibility_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        group_col = "LSI_harmony_SNN_cluster_cell_type"
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    tarchetypes::tar_file(
      name = motif_family_accessibility_marker_volcano_plots.ATAC,
      description = "Facetted volcano plot of marker motif-family accessibility per ATAC cell type. [checkpoint:ATAC]",
      command = {
        plot <- motif_family_accessibility_markers.ATAC |>
          dplyr::mutate(avg_log2FC = .data$avg_diff) |>
          plot_markers_volcano_simple() +
          ggplot2::labs(x = "Mean motif-family accessibility difference")
        save_plots_structured(plot)
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = motif_family_accessibility_heatmap.ATAC,
      description = "Heatmap of configured marker motif-family accessibility per ATAC cell type. [checkpoint:ATAC]",
      command = plot_feature_scores_heatmap_from_matrix(
        feature_matrix = motif_family_accessibility_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        features = ATAC_marker_motif_families_vec,
        group_col = "LSI_harmony_SNN_cluster_cell_type"
      ) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = continuous_UMAP_spec.ATAC,
      description = "Continuous ATAC UMAP variables to plot one at a time",
      command = {
        feature_vars <- intersect(ATAC_marker_motif_families_vec, rownames(motif_family_accessibility_BPCells_matrix.ATAC))
        tibble::tibble(
          variable = c(aggregation_w_peaks_continuous_vars, feature_vars),
          value_source = c(
            rep("metadata", length(aggregation_w_peaks_continuous_vars)),
            rep("feature", length(feature_vars))
          )
        )
      },
      iteration = "vector"
    ),
    tarchetypes::tar_file(
      name = continuous.UMAPs.ATAC,
      description = "UMAPs colored by continuous motif-family accessibility and peak accessibility metrics. [checkpoint:ATAC]",
      command = plot_UMAP_from_metadata(
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        variable = continuous_UMAP_spec.ATAC$variable,
        value_source = continuous_UMAP_spec.ATAC$value_source,
        feature_matrix = motif_family_accessibility_BPCells_matrix.ATAC
      ) |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = stringr::str_replace_all(continuous_UMAP_spec.ATAC$variable, "[/\\\\]", "_")
        ),
      pattern = map(continuous_UMAP_spec.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    )
  ),
  ATAC_w_gene_activity_targets = rlang::list2(
    tarchetypes::tar_file(
      name = gene_score_archr_BPCells_matrix_dir.ATAC,
      description = "Compute ArchR-style ATAC gene activity scores with BPCells and write them to disk",
      command = {
        out_dir <- get_structured_file_path()
        tile_matrix_dir <- tempfile(pattern = "gene_score_tiles_", tmpdir = dirname(out_dir))
        on.exit(if (fs::dir_exists(tile_matrix_dir)) fs::dir_delete(tile_matrix_dir), add = TRUE)

        fragments <- BPCells::select_cells(
          combined_BPCells_fragment_obj.ATAC,
          metadata_w_cell_types_tibble.ATAC$barcode_w_prefix
        )
        genes <- marker_validated_Ensembl_annotations_GRanges_list$genes
        chromosome_sizes <- get_chrom_sizes_for_BPCells_tile_calling(aggregated_cellranger_ref_list$genomes[[1]])
        genes <- filter_GRanges_to_chrom_sizes(genes, chromosome_sizes)
        blacklist <- filter_GRanges_to_chrom_sizes(blacklist_GRanges.ATAC, chromosome_sizes, allow_empty = TRUE)

        if (!"gene_name" %in% colnames(S4Vectors::mcols(genes))) {
          genes$gene_name <- names(genes)
        }
        # Put genes with no retained ArchR tiles first so BPCells can keep a
        # full zero row for them after blacklist filtering.
        genes <- order_genes_for_BPCells_gene_score_archr(
          genes = genes,
          chromosome_sizes = chromosome_sizes,
          blacklist = blacklist
        )
        gene_scores <- BPCells::gene_score_archr(
          fragments = fragments,
          genes = genes,
          chromosome_sizes = chromosome_sizes,
          blacklist = blacklist,
          gene_name_column = "gene_name",
          tile_matrix_path = tile_matrix_dir
        )

        BPCells::write_matrix_dir(gene_scores, out_dir, overwrite = TRUE)
        out_dir
      },
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = gene_score_archr_BPCells_matrix.ATAC,
      description = "Open the BPCells ArchR-style ATAC gene activity score matrix",
      command = BPCells::open_matrix_dir(gene_score_archr_BPCells_matrix_dir.ATAC),
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    tarchetypes::tar_file(
      name = marker_gene_activity_dot_plot.ATAC,
      description = "Dot plot of marker gene activity scores per ATAC cell type. [checkpoint:ATAC]",
      command = {
        plot <- plot_marker_gene_activity_dot_BPCells(
          feature_matrix = gene_score_archr_BPCells_matrix.ATAC,
          metadata_tibble = metadata_w_cell_types_tibble.ATAC,
          features = GEX_marker_genes_vec,
          group_col = "LSI_harmony_SNN_cluster_cell_type"
        )
        save_plots_structured(plot)
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    )
  )
)

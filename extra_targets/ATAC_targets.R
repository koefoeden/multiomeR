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
      resources = get_tar_resources(cores_req = 1, RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = combined_BPCells_fragment_obj.ATAC,
      description = "Merge per-reaction BPCells fragment objects and filter to standard chromosomes",
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
      description = "Call full-genome ATAC peaks per cluster using capped discovery cells",
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
      description = "Collapse per-cluster ATAC peaks into non-overlapping peak sets",
      command = {
        collapsed_peaks <- combine_collapse_GRanges_ArchR(list(peak_GRanges_per_cluster.ATAC), by = "neg_log10pvalue_summit")
        format_peak_GRanges(collapsed_peaks, cluster_name = peak_calling_cluster_names.ATAC)
      },
      pattern = map(peak_GRanges_per_cluster.ATAC, peak_calling_cluster_names.ATAC)
    ),
    tarchetypes::tar_file(
      name = peaks_similarity_tiles_plot.ATAC,
      description = "Plot pairwise peak-set similarity matrix across clusters. [checkpoint:multimodal]",
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
      description = "Iteratively collapse per-cluster peak sets into a single consensus peak set",
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
      description = "Compute the consensus peak-by-cell count matrix using BPCells",
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

  ATAC_consensus_peak_files_targets = rlang::list2(
    tarchetypes::tar_file(
      name = consensus_peak_BED.ATAC,
      description = "Write consensus peak set as a BED file",
      command = write_peak_bed_files_from_GRanges(consensus_peak_GRanges.ATAC, cluster_name = "")
    ),
    tarchetypes::tar_file(
      name = consensus_peak_bigBed.ATAC,
      description = "Convert consensus BED to indexed bigBed format",
      command = write_bigBed_from_BED_file(consensus_peak_BED.ATAC, cluster_name = "", genome = aggregated_cellranger_ref_list$genomes[[1]])
    ),
    # tar_target( # disabled while testing.
    #   name = consensus_peak_bigBed_galaxy.ATAC,
    #   description = "Upload consensus bigBed track to Galaxy",
    #   command = copy_to_galaxy(consensus_peak_bigBed.ATAC, aggregation_GALAXY_track_upload_HISTORY_ID, aggregation_GALAXY_track_upload_API_KEY)
    # )
  ),
  ATAC_coverage_plots_targets = rlang::list2(
    targets::tar_target(
      name = coverage_regions_tibble.ATAC,
      description = "Build genomic regions to display in ATAC coverage tracks",
      command = get_coverage_regions_tibble_BPCells(
        formatted_peak_tibble = consensus_peak_tibble.ATAC,
        gene_GRanges = marker_validated_Ensembl_annotations_GRanges_list$genes,
        marker_genes = GEX_marker_genes_vec
      ),
      iteration = "vector"
    ),
    tarchetypes::tar_file(
      name = coverage_tracks_plots.ATAC,
      description = "Plot ATAC coverage tracks at marker gene and DA peak loci. [checkpoint:multimodal]",
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
      description = "Run BPCells-native TF-IDF and SVD on the QC-filtered ATAC peak matrix",
      command = run_ATAC_LSI_BPCells(
        ATAC_peak_BPCells_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
        n_components = utils::tail(aggregation_ATAC_data_PCs, 1),
        threads = 15
      ),
      resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
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
      description = "Harmony-corrected BPCells-native ATAC LSI embeddings",
      command = run_harmony_on_embedding_matrix(
        embedding_matrix = LSI_BPCells.ATAC$cell_embeddings,
        metadata_tibble = metadata_filtered_tibble.ATAC,
        harmony_correction_metadata_col_names = c(aggregation_harmony_correction_metadata_col_names, aggregation_extra_harmony_covars_ATAC),
        dims = aggregation_ATAC_data_PCs,
        cores = 15
      ),
      resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = UMAP_embeddings_tibble.ATAC,
      description = "UMAP coordinates from BPCells-native Harmony-corrected ATAC LSI embeddings",
      command = run_UMAP_from_embedding_matrix(
        embedding_matrix = harmony_embeddings_matrix.ATAC,
        dims = aggregation_UMAP_ATAC_PCs,
        n_neighbors = aggregation_UMAP_nNNs,
        min_dist = aggregation_UMAP_min_dist
      ) |>
        tibble::as_tibble(rownames = "barcode_w_prefix"),
      resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = LSI_clusters.ATAC,
      description = "Leiden clusters from a BPCells nearest-neighbour graph on Harmony-corrected ATAC LSI embeddings",
      command = cluster_embedding_matrix_BPCells(
        embedding_matrix = harmony_embeddings_matrix.ATAC,
        dims = aggregation_ATAC_data_PCs,
        k = aggregation_data_nNNs,
        resolution = aggregation_ATAC_cluster_res,
        threads = 15,
        min_barcodes = aggregation_cluster_min_barcodes
      ),
      resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
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
      description = "Assign cell type labels to BPCells-native ATAC clusters from metadata module scores",
      command = add_cell_types_to_metadata_from_module_scores(
        metadata_tibble = metadata_w_clusters_tibble.ATAC,
        named_marker_genes_list = UCell_GEX_marker_genes_list,
        allow_multiple_cell_types = aggregation_allow_multiple_cell_types,
        cluster_column = "LSI_harmony_SNN_cluster"
      )
    ),
    tarchetypes::tar_file(
      name = VizDimLoadings_plots.ATAC,
      description = "Plot top feature loadings for each LSI dimension. [checkpoint:multimodal]",
      command = plot_LSI_loadings_from_tibble(
        LSI_loadings_tibble = LSI_loadings_tibble.ATAC,
        dims = aggregation_ATAC_data_PCs,
        nfeatures = 50
      ) |>
        save_plots_structured()
    ),
    tarchetypes::tar_file(
      name = LSI_singular_values_elbow_plot.ATAC,
      description = "Elbow plot of native ATAC LSI singular values. [checkpoint:multimodal]",
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
      description = "Non-Harmony and Harmony ATAC LSI embedding coordinate SD plot. [checkpoint:multimodal]",
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
      description = "Non-Harmony and Harmony ATAC LSI metadata association bar plots. [checkpoint:multimodal]",
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
            "TENX_reaction_ID",
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
      name = scDblFinder_reaction_tibble.ATAC,
      description = "Prepare per-reaction ATAC barcodes and cluster labels for scDblFinder",
      command = prepare_scDblFinder_reaction_tibble(
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
        threads = 15
      ),
      resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = scDblFinder_results_by_reaction_tibble.ATAC,
      description = "Run ATAC scDblFinder independently per 10x reaction from pre-aggregated peak-group counts",
      command = run_scDblFinder_BPCells_reaction(
        feature_matrix = scDblFinder_feature_matrix.ATAC,
        scDblFinder_reaction_tibble = scDblFinder_reaction_tibble.ATAC,
        output_suffix = "ATAC",
        dbr.sd = 1.0,
        aggregateFeatures = FALSE,
        nfeatures = 50,
        processing = "normFeatures"
      ),
      pattern = map(scDblFinder_reaction_tibble.ATAC),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    # Previous full-peak ATAC scDblFinder path. Kept as fallback, but disabled because
    # scDblFinder::aggregateFeatures() recomputes TF-IDF and can OOM before aggregation.
    # tar_target(
    #   name = scDblFinder_results_by_reaction_tibble.ATAC,
    #   description = "Run ATAC scDblFinder independently per 10x reaction from BPCells peak-count slices",
    #   command = run_scDblFinder_BPCells_reaction(
    #     feature_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
    #     scDblFinder_reaction_tibble = scDblFinder_reaction_tibble.ATAC,
    #     output_suffix = "ATAC",
    #     dbr.sd = 1.0,
    #     aggregateFeatures = TRUE,
    #     nfeatures = 50,
    #     processing = "normFeatures"
    #   ),
    #   pattern = map(scDblFinder_reaction_tibble.ATAC),
    #   resources = get_tar_resources(RAM_GB_req = 16)
    # ),
    targets::tar_target(
      name = scDblFinder_results_df.ATAC,
      description = "Combine per-reaction ATAC scDblFinder classifications",
      command = scDblFinder_results_by_reaction_tibble.ATAC |>
        dplyr::bind_rows() |>
        dplyr::select(-dplyr::any_of("TENX_reaction_ID")),
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    tarchetypes::tar_file(
      name = scDblFinder_score_violins_plot.ATAC,
      description = "Violin plots of scDblFinder doublet scores by cluster. [checkpoint:multimodal]",
      command = {
        scDblFinder_metadata <- metadata_w_cell_types_unfiltered_tibble.ATAC |>
          dplyr::left_join(scDblFinder_results_df.ATAC, by = "barcode_w_prefix")
        plot_data <- scDblFinder_metadata |>
          dplyr::select(LSI_harmony_SNN_cluster, LSI_harmony_SNN_cluster_cell_type, scDblFinder.score_ATAC, scDblFinder.class_ATAC, TENX_reaction_ID) |>
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
            ggplot2::aes(color = TENX_reaction_ID),
            size = 0.3,
            alpha = 0.5,
            width = 0.2
          ) +
          ggplot2::facet_wrap(~cluster_type, scales = "free_x") +
          ggplot2::labs(subtitle = "Points are cells classified as doublet by scDblFinder, colored by reaction.") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::theme(legend.position = "none")

        save_plots_structured(plot)
      }
    ),
    targets::tar_target(
      name = metadata_w_cell_types_tibble.ATAC,
      description = "Remove doublets and high-doublet clusters from BPCells-native ATAC metadata",
      command = {
        excluded_doublets_barcodes <- dplyr::filter(
          scDblFinder_results_df.ATAC,
          .data$scDblFinder.class_ATAC == "doublet"
        ) |>
          dplyr::pull(.data$barcode_w_prefix)

        high_doublet_clusters_vec <- metadata_w_cell_types_unfiltered_tibble.ATAC |>
          dplyr::mutate(is_excluded = .data$barcode_w_prefix %in% excluded_doublets_barcodes) |>
          dplyr::group_by(.data$LSI_harmony_SNN_cluster) |>
          dplyr::summarise(frac_excluded = mean(is_excluded), .groups = "drop") |>
          dplyr::filter(frac_excluded > 0.5) |>
          dplyr::pull(.data$LSI_harmony_SNN_cluster)

        excluded_high_doublet_cluster_barcodes_vec <- metadata_w_cell_types_unfiltered_tibble.ATAC |>
          dplyr::filter(.data$LSI_harmony_SNN_cluster %in% high_doublet_clusters_vec) |>
          dplyr::pull(.data$barcode_w_prefix)

        metadata_w_cell_types_unfiltered_tibble.ATAC |>
          dplyr::filter(!.data$barcode_w_prefix %in% base::union(excluded_doublets_barcodes, excluded_high_doublet_cluster_barcodes_vec)) |>
          dplyr::left_join(scDblFinder_results_df.ATAC, by = "barcode_w_prefix")
      }
    )
  ),

  ATAC_QC_plots_targets = rlang::list2(
    tarchetypes::tar_file(
      name = peaks_QC_violins_plot.ATAC,
      description = "Violin plots of peak-based ATAC QC metrics per reaction. [checkpoint:multimodal]",
      command = {
        plot_data <- metadata_w_cell_types_tibble.ATAC |>
          dplyr::select(TENX_reaction_ID, "dataset", dplyr::any_of(aggregation_peak_based_continuous_QC_vars)) |>
          tidyr::pivot_longer(cols = dplyr::any_of(aggregation_peak_based_continuous_QC_vars), names_to = "feature", values_to = "value")

        plot <- plot_data |>
          ggplot2::ggplot(ggplot2::aes(x = TENX_reaction_ID, y = value, fill = dataset)) +
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
      description = "Bar plots of categorical metadata proportions by cell-type cluster. [checkpoint:multimodal]",
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
      description = "UpSet plot of overlapping ATAC QC exclusion reasons. [checkpoint:multimodal]",
      command = {
        plot <- plot_upset_from_excluded_BCs_list(QC_excluded_BCs_list.ATAC, n_total = nrow(metadata_w_QC_tibble.ATAC))
        save_plots_structured(plot)
      }
    ),
    tarchetypes::tar_file(
      name = categorical.UMAPs.ATAC,
      description = "UMAPs colored by categorical metadata variables. [checkpoint:multimodal]",
      command = {
        plots <- plot_UMAP_from_metadata(metadata_w_cell_types_tibble.ATAC, metadata_cols = aggregation_ATAC_categorical_vars)
        save_plots_structured(plots)
      }
    ),
    targets::tar_target(
      name = UMAP_n_dims_seq.ATAC,
      description = "Sequence of LSI dimension counts for cross-parameter UMAP sweep",
      command = round(seq(5, length(aggregation_ATAC_data_PCs), length.out = 3))
    ),
    tarchetypes::tar_file(
      name = cross.UMAPs.ATAC,
      description = "Compute ATAC UMAPs across a sweep of LSI dimension counts and neighbour counts. [checkpoint:multimodal]",
      command = {
        sweep_umap <- run_UMAP_from_embedding_matrix(
          embedding_matrix = harmony_embeddings_matrix.ATAC[metadata_w_cell_types_tibble.ATAC$barcode_w_prefix, , drop = FALSE],
          dims = 1:UMAP_n_dims_seq.ATAC,
          n_neighbors = UMAP_neighbors_seq,
          min_dist = aggregation_UMAP_min_dist
        ) |>
          tibble::as_tibble(rownames = "barcode_w_prefix")

        plots <- metadata_w_cell_types_tibble.ATAC |>
          dplyr::select(-dplyr::any_of(c("LSI_UMAP_1", "LSI_UMAP_2"))) |>
          dplyr::left_join(sweep_umap, by = "barcode_w_prefix") |>
          plot_UMAP_from_metadata(metadata_cols = "LSI_harmony_SNN_cluster_cell_type")

        save_plots_structured(plots, dyn_suffix_in_subdir = TRUE, override_suffix = paste0(UMAP_n_dims_seq.ATAC, "_", UMAP_neighbors_seq))
      },
      pattern = cross(UMAP_n_dims_seq.ATAC, UMAP_neighbors_seq),
      resources = get_tar_resources(RAM_GB_req = 60)
    )
  ),

  ATAC_w_TF_activity_targets = rlang::list2(
    targets::tar_target(
      name = ATAC_marker_TFs_list,
      description = "Validate configured ATAC marker TFs against organism motif names",
      command = {
        motif_list <- TF_motif_matrix_list
        motif_names <- names(motif_list)
        motif_symbols <- stringr::str_split_i(motif_names, "__", 1)
        resolve_marker_TFs <- function(marker_TFs) {
          purrr::map_chr(
            marker_TFs,
            \(marker_TF) {
              marker_TF_suffix <- stringr::str_extract(marker_TF, "[+-]$")
              marker_TF_name <- stringr::str_remove_all(marker_TF, "[+-]$")
              full_match <- motif_names[stringr::str_to_upper(motif_names) == stringr::str_to_upper(marker_TF_name)]
              if (length(full_match) == 1) {
                return(paste0(full_match, dplyr::coalesce(marker_TF_suffix, "")))
              }

              symbol_matches <- motif_names[stringr::str_to_upper(motif_symbols) == stringr::str_to_upper(marker_TF_name)]
              if (length(symbol_matches) == 0) {
                stop(
                  "Configured ATAC marker TF is missing from vertebrate motif names in ",
                  targets::tar_name(),
                  ": ",
                  marker_TF_name,
                  ". Use a motif symbol or ID-qualified motif feature from TF_motif_matrix_list.",
                  call. = FALSE
                )
              }
              if (length(symbol_matches) > 1) {
                stop(
                  "Configured ATAC marker TF is ambiguous in ",
                  targets::tar_name(),
                  ": ",
                  marker_TF_name,
                  ". Available ID-qualified motif features: ",
                  paste(symbol_matches, collapse = ", "),
                  ". Update aggregation_ATAC_marker_TFs in cfg_aggregations.yaml to one of these exact feature names.",
                  call. = FALSE
                )
              }
              paste0(symbol_matches, dplyr::coalesce(marker_TF_suffix, ""))
            }
          )
        }
        purrr::map(aggregation_ATAC_marker_TFs, resolve_marker_TFs)
      }
    ),
    targets::tar_target(
      name = ATAC_marker_TFs_vec,
      description = "Flatten validated ATAC marker TFs to plain feature names",
      command = {
        ATAC_marker_TFs_list |>
          unlist(use.names = FALSE) |>
          stringr::str_remove_all("[+-]$") |>
          unique()
      }
    ),
    targets::tar_target(
      name = peak_TF_motif_matrix.ATAC, # TODO: this runs quite slowly for large peak sets - consider looking into faster alternatives.
      description = "Match TF motifs to ATAC peaks for reuse by betterChromVAR",
      command = {
        ATAC_marker_TFs_list
        get_motif_matrix_from_ATAC_peak_names(
          ATAC_peak_names = rownames(consensus_peak_BPCells_matrix.ATAC),
          ATAC_peak_GRanges = consensus_peak_GRanges.ATAC,
          TF_motif_matrix_list = TF_motif_matrix_list,
          genome = aggregated_cellranger_ref_list$genomes[[1]]
        )
      },
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = chromVAR_obj.ATAC,
      description = "Build the lightweight BPCells-backed chromVAR RSE for TF activity and genetic enrichment",
      command = get_chromVAR_obj_from_peak_matrix(
        ATAC_peak_matrix = consensus_peak_BPCells_matrix.ATAC[, metadata_w_cell_types_tibble.ATAC$barcode_w_prefix],
        ATAC_peak_GRanges = consensus_peak_GRanges.ATAC,
        genome = aggregated_cellranger_ref_list$genomes[[1]]
      ),
      resources = get_tar_resources(RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = chromVAR_peak_expectation.ATAC,
      description = "Compute shared peak-level accessibility expectations for chromVAR analyses",
      command = get_chromVAR_peak_expectation(chromVAR_obj.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = chromVAR_background_bins.ATAC,
      description = "Compute shared peak-level betterChromVAR background bins for chromVAR analyses",
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
      name = chromVAR_TF_motif_matrix.ATAC,
      description = "Align the TF motif annotation matrix to the filtered chromVAR peak order",
      command = align_chromVAR_annotations_to_obj(
        annotations = peak_TF_motif_matrix.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC
      ),
      resources = get_tar_resources(RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = motif_chromVAR_chunk_results.ATAC,
      description = "Compute per-cell TF motif accessibility scores for one reusable ATAC chunk",
      command = {
        ATAC_marker_TFs_list
        compute_chromVAR_annotation_chunk_result(
          annotations = chromVAR_TF_motif_matrix.ATAC,
          chromVAR_obj = chromVAR_obj.ATAC,
          chunk_context_record = chromVAR_chunk_context_records.ATAC,
          compute = c("deviations", "z")
        )
      },
      pattern = map(chromVAR_chunk_context_records.ATAC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = motif_chromVAR_results.ATAC,
      description = "Combine chunk-level TF motif accessibility scores into the standard chromVAR result object",
      command = combine_chromVAR_chunk_results(
        chunk_results = motif_chromVAR_chunk_results.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC,
        annotations = chromVAR_TF_motif_matrix.ATAC
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    tarchetypes::tar_file(
      name = TF_activity_BPCells_matrix_dir.ATAC,
      description = "Write betterChromVAR TF activity Z-scores to a BPCells matrix directory",
      command = {
        out_dir <- get_structured_file_path()
        column_major_dir <- tempfile(pattern = "TF_activity_column_major_", tmpdir = dirname(out_dir))
        on.exit(if (fs::dir_exists(column_major_dir)) fs::dir_delete(column_major_dir), add = TRUE)
        if (fs::dir_exists(out_dir)) {
          fs::dir_delete(out_dir)
        }

        motif_chromVAR_results.ATAC$chromVAR_z_scores |>
          Matrix::Matrix(sparse = TRUE) |>
          methods::as("dgCMatrix") |>
          BPCells::write_matrix_dir(column_major_dir, compress = FALSE, overwrite = TRUE)
        BPCells::transpose_storage_order(BPCells::open_matrix_dir(column_major_dir), outdir = out_dir)
        out_dir
      },
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = TF_activity_BPCells_matrix.ATAC,
      description = "Open the disk-backed betterChromVAR TF activity Z-score matrix",
      command = BPCells::open_matrix_dir(TF_activity_BPCells_matrix_dir.ATAC),
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    targets::tar_target(
      name = TF_activity_markers.ATAC,
      description = "Find marker TF activities per ATAC cell type from betterChromVAR Z-scores",
      command = get_marker_TF_activities_from_chromVAR_BPCells_z_scores(
        chromVAR_z_scores_BPCells_matrix = TF_activity_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        group_col = "LSI_harmony_SNN_cluster_cell_type"
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    tarchetypes::tar_file(
      name = TF_activity_marker_volcano_plots.ATAC,
      description = "Facetted volcano plot of marker TF activity scores per ATAC cell type. [checkpoint:multimodal]",
      command = {
        plot <- TF_activity_markers.ATAC |>
          dplyr::mutate(avg_log2FC = .data$avg_diff) |>
          plot_markers_volcano_simple() +
          ggplot2::labs(x = "Mean TF activity difference")
        save_plots_structured(plot)
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = TF_activity_dot_plot.ATAC,
      description = "Dot plot of configured marker TF activity scores per ATAC cell type. [checkpoint:multimodal]",
      command = plot_feature_scores_dot_from_matrix(
        feature_matrix = TF_activity_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        features = ATAC_marker_TFs_vec,
        group_col = "LSI_harmony_SNN_cluster_cell_type"
      ) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = continuous.UMAPs.ATAC,
      description = "UMAPs colored by continuous TF activity and peak accessibility metrics. [checkpoint:multimodal]",
      command = plot_UMAP_from_metadata(
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        metadata_cols = aggregation_w_peaks_continuous_vars,
        feature_matrix = TF_activity_BPCells_matrix.ATAC,
        feature_rows = ATAC_marker_TFs_vec
      ) |>
        save_plots_structured(),
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
      description = "Dot plot of marker gene activity scores per ATAC cell type. [checkpoint:multimodal]",
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
  ),
  ATAC_tracks_targets = rlang::list2(
    tarchetypes::tar_file(
      name = peaks_per_cluster_BEDs.peaks.ATAC,
      description = "Write per-cluster collapsed peak sets as BED files",
      command = write_peak_bed_files_from_GRanges(
        within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC,
        cluster_name = peak_calling_cluster_names.ATAC
      ),
      pattern = map(within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC, peak_calling_cluster_names.ATAC)
    ),
    tarchetypes::tar_file(
      name = peaks_per_cluster_bigBeds.peaks.ATAC,
      description = "Convert per-cluster peak BEDs to indexed bigBed format",
      command = write_bigBed_from_BED_file(
        peaks_per_cluster_BEDs.peaks.ATAC,
        cluster_name = peak_calling_cluster_names.ATAC,
        genome = aggregated_cellranger_ref_list$genomes[[1]]
      ),
      pattern = map(peaks_per_cluster_BEDs.peaks.ATAC, peak_calling_cluster_names.ATAC)
    ),
    # tar_target( # disabled while testing.
    #   name = peaks_per_cluster_bigBeds_galaxy.ATAC,
    #   description = "Upload per-cluster bigBed peak tracks to Galaxy",
    #   command = copy_to_galaxy(peaks_per_cluster_bigBeds.peaks.ATAC, aggregation_GALAXY_track_upload_HISTORY_ID, aggregation_GALAXY_track_upload_API_KEY),
    #   pattern = map(peaks_per_cluster_bigBeds.peaks.ATAC)
    # ),
    tarchetypes::tar_file(
      name = coverage_per_cluster_bedgraphs.coverage.ATAC,
      description = "Compute CPM-normalized Tn5 insertion bedGraphs per peak-calling cluster with BPCells",
      command = {
        out_paths <- purrr::map_chr(
          peak_calling_cluster_names.ATAC,
          \(cluster_name) {
            get_structured_file_path(
              filetype = "bedgraph",
              override_suffix = cluster_name
            )
          }
        )
        names(out_paths) <- peak_calling_cluster_names.ATAC

        if (!identical(names(BCs_per_peak_cluster_list.ATAC), peak_calling_cluster_names.ATAC)) {
          stop("Peak-calling barcode groups and cluster names are not in the same order.")
        }

        cell_group_lookup <- stats::setNames(
          rep(peak_calling_cluster_names.ATAC, lengths(BCs_per_peak_cluster_list.ATAC)),
          unlist(BCs_per_peak_cluster_list.ATAC, use.names = FALSE)
        )
        fragments <- BPCells::select_cells(combined_BPCells_fragment_obj.ATAC, names(cell_group_lookup))
        cell_groups <- cell_group_lookup[BPCells::cellNames(fragments)]

        if (anyNA(cell_groups)) {
          stop("Missing peak-calling cluster labels for one or more BPCells fragment barcodes.")
        }

        BPCells::write_insertion_bedgraph(
          fragments = fragments,
          path = out_paths,
          cell_groups = factor(cell_groups, levels = peak_calling_cluster_names.ATAC),
          normalization_method = "cpm",
          chrom_sizes = get_chrom_sizes_for_BPCells_tile_calling(aggregated_cellranger_ref_list$genomes[[1]])
        )
        out_paths
      }
    ),
    tarchetypes::tar_files(
      name = coverage_per_cluster_bedgraphs_files.coverage.ATAC,
      description = "Convert per-cluster insertion bedgraphs to bigWig format",
      command = coverage_per_cluster_bedgraphs.coverage.ATAC,
      deployment = "main"
    )
  ),
  GALAXY_ATAC_tracks_targets = rlang::list2(
    targets::tar_target(
      name = coverage_per_cluster_BWs_galaxy.ATAC,
      description = "Upload per-cluster bigWig coverage tracks to Galaxy",
      command = copy_to_galaxy(
        coverage_per_cluster_bedgraphs_files.coverage.ATAC,
        aggregation_GALAXY_track_upload_HISTORY_ID,
        aggregation_GALAXY_track_upload_API_KEY
      ),
      pattern = map(coverage_per_cluster_bedgraphs_files.coverage.ATAC)
    ),
    tarchetypes::tar_file(
      name = tracks_hub_file_galaxy.ATAC,
      description = "Write UCSC track hub file pointing to Galaxy-hosted ATAC tracks",
      command = write_hub_txt(
        bigWig_files = coverage_per_cluster_BWs_galaxy.ATAC,
        peak_files = peaks_per_cluster_bigBeds_galaxy.ATAC,
        ATAC_peak_calling_cluster_names = peak_calling_cluster_names.ATAC,
        consensus_peak_file = consensus_peak_bigBed_galaxy.ATAC
      )
    ),
    targets::tar_target(
      name = tracks_hub_galaxy.ATAC,
      description = "Upload track hub file to Galaxy",
      command = copy_to_galaxy(tracks_hub_file_galaxy.ATAC, aggregation_GALAXY_track_upload_HISTORY_ID, aggregation_GALAXY_track_upload_API_KEY)
    ),
    tarchetypes::tar_file(
      name = tracks_hub_galaxy_link.ATAC,
      description = "Save Galaxy track hub URL to a text file",
      command = {
        file_out <- get_structured_file_path(filetype = "link.txt")
        readr::write_lines(tracks_hub_galaxy.ATAC, file = file_out)
        file_out
      }
    )
  ),

  ATAC_tracks_w_scE2G_targets = rlang::list2(
    tarchetypes::tar_file(
      name = tracks_w_scE2G_tracks_hub_file.ATAC,
      description = "Write UCSC track hub including scE2G enhancer-gene interaction tracks",
      command = write_hub_txt(
        bigWig_files = coverage_per_cluster_BWs_galaxy.ATAC,
        peak_files = peaks_per_cluster_bigBeds_galaxy.ATAC,
        ATAC_peak_calling_cluster_names = peak_calling_cluster_names.ATAC,
        bigInteract_files = scE2G_bigInteract_galaxy,
        consensus_peak_file = consensus_peak_bigBed_galaxy.ATAC
      )
    ),
    targets::tar_target(
      name = tracks_w_scE2G_tracks_hub_galaxy.ATAC,
      description = "Upload scE2G-inclusive track hub to Galaxy",
      command = copy_to_galaxy(tracks_w_scE2G_tracks_hub_file.ATAC, aggregation_GALAXY_track_upload_HISTORY_ID, aggregation_GALAXY_track_upload_API_KEY)
    ),
    tarchetypes::tar_file(
      name = tracks_w_scE2G_tracks_hub_galaxy_link.ATAC,
      description = "Save Galaxy scE2G hub URL to a text file",
      command = {
        file_out <- get_structured_file_path(filetype = "link.txt")
        readr::write_lines(tracks_w_scE2G_tracks_hub_galaxy.ATAC, file = file_out)
        file_out
      }
    )
  )
)

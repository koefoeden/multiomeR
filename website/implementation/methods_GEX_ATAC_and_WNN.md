# GEX, ATAC, batch correction and WNN

This chapter covers normalization and dimensional reduction of both modalities, peak definition, batch correction, graph construction and clustering, and weighted nearest-neighbour (WNN) integration. The target structure is shown in the [primary-module graph](implementation_main.md), and the WNN implementation is compared with Seurat in [Algorithmic implementations](algorithm_validation.md#native-weighted-nearest-neighbors). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

## GEX normalization and PCA

{{< include _shared_methods/GEX_normalization.md >}}

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Gene filter | Minimum total count per gene, both backends | Fixed | more than 50 counts over the PCA nuclei (`min_feature_count = 50`) | `run_GEX_PCA_BPCells()` |
| Backend | Normalization and PCA backend | Configurable |  | [`aggregation_GEX_PCA_backend`](../parameters.html#aggregation_GEX_PCA_backend) |
| Regression | Cell-level covariates regressed | Configurable |  | [`aggregation_SCT_regress_vars`](../parameters.html#aggregation_SCT_regress_vars) |
| Components | Number of PCs computed (last element of the list) | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs) |
| Variable genes | Genes retained by residual variance, both backends | Fixed | 3,000 (`n_variable_features = 3000`) | `run_GEX_PCA_BPCells()` |
| BPCells branch | Pearson residuals | Fixed | method-of-moments theta clamped to \[1e-6, 1e6\]; `clip_range = c(-10, 10)`; `min_var = 0` | `run_BPCells_native_GEX_PCA()` |
| BPCells branch | Regression and SVD | Fixed | `BPCells::regress_out(prediction_axis = "row")` on the residuals; `BPCells::svds()` without centring, embeddings = right singular vectors × singular values | `run_BPCells_native_GEX_PCA()` |
| Seurat branch | SCTransform | Fixed | `conserve.memory = TRUE`, `do.correct.umi = FALSE`; package defaults otherwise, including `vst.flavor = "v2"` and `ncells = 5000` | `run_Seurat_SCT_for_PCA()`, `Seurat::SCTransform()` |
| Seurat branch | PCA | Fixed | exact eigendecomposition of the gene × gene residual Gram matrix | `run_dense_feature_gram_PCA()` |
| Cell cycle | Gene sets and scoring | Fixed | `Seurat::cc.genes.updated.2019` for both organisms; log-normalization with scale factor 10,000; binned control-gene scores with seed 1; `CC.Difference = S − G2M` | `cell_cycle_gene_sets()`, `add_cell_cycle_scores_to_cell_attr()`, `calculate_BPCells_cell_cycle_scores_from_matrix()` |

## ATAC peak definition and dimensional reduction

{{< include _shared_methods/ATAC_peaks_and_LSI.md >}}

### ATAC peak calling and consensus peaks

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Fragments | Chromosomes retained | Fixed | autosomes, X and Y | `combined_BPCells_fragment_obj.ATAC` |
| Groups | Peak-calling grouping column | Configurable |  | [`aggregation_call_peaks_by_cluster_col`](../parameters.html#aggregation_call_peaks_by_cluster_col) |
| Groups | Discovery cap per group | Fixed | 50,000 nuclei; larger groups randomly downsampled with seed 1 + group index | `build_peak_calling_cluster_discovery_tibble()` |
| Caller | Peak-calling method | Configurable |  | [`aggregation_ATAC_peak_calling_method`](../parameters.html#aggregation_ATAC_peak_calling_method) |
| MACS3 | Arguments | Fixed | `-f BED --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all`; MACS3 default q-value cutoff 0.05 | `call_peaks_w_MACS3()` |
| Genome | Effective genome size, used by both callers and the peak-count enrichment | Fixed | GRCh38 2.913e9; mm10 and GRCm39 2.65e9 | `get_effective_genome_size()` |
| Tile caller | BPCells settings | Fixed | `peak_width = 500`, `peak_tiling = 3`, `fdr_cutoff = 0.01`, `merge_peaks = "none"` | `call_peaks_w_BPCells_tile()` |
| Peak shape | Summit extension | Fixed | 250 bp either side of the summit, giving 500 bp peaks | `get_peak_GRanges_w_fixed_width()` |
| Blacklist | Source and rule | Fixed | GRCh38 `hg38.Kundaje.GRCh38_unified_Excludable` and mm10 `mm10.Boyle.mm10-Excludable.v2` from AnnotationHub excluderanges; GRCm39 `resources/mm39.excluderanges.bed`; any overlap removes a peak | `get_blacklist_GRanges()`, `get_peak_GRanges_w_fixed_width()` |
| Consensus | Overlap ranking | Fixed | within groups `neg_log10pvalue_summit`, across groups `fold_change`, both decreasing | `within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC`, `consensus_peak_GRanges.ATAC` |
| Peak matrix | Counting mode, also used for blacklist QC counts | Configurable |  | [`aggregation_ATAC_peak_matrix_mode`](../parameters.html#aggregation_ATAC_peak_matrix_mode) |
| Peak matrix | Nuclei | Fixed | GEX nuclei after doublet removal | `consensus_peak_BPCells_matrix_dir.ATAC` |

### ATAC TF-IDF and LSI

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| TF-IDF | Transform | Fixed | `log1p(10000 × TF × IDF)`, with TF = count / nucleus total and IDF = number of nuclei / peak total | `run_ATAC_LSI_BPCells()` |
| SVD | Components computed (last element of the list) | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| SVD | Decomposition | Fixed | `BPCells::svds()` without centring; embeddings = right singular vectors × singular values | `run_ATAC_LSI_BPCells()` |
| Dimensions | Components used downstream | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |

## Batch correction, clustering and weighted nearest neighbours

{{< include _shared_methods/batch_correction_clustering_WNN.md >}}

### Harmony batch correction

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Covariates | Shared correction columns | Configurable |  | [`aggregation_harmony_correction_metadata_col_names`](../parameters.html#aggregation_harmony_correction_metadata_col_names) |
| Covariates | Additional ATAC covariates | Configurable |  | [`aggregation_extra_harmony_covars_ATAC`](../parameters.html#aggregation_extra_harmony_covars_ATAC) |
| Dimensions | Corrected dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| Harmony | Arguments | Fixed | `max_iter = 25`, `lambda = 1`; package defaults otherwise | `run_harmony_on_embedding_matrix()`, `harmony::RunHarmony()` |

### Graph construction, Leiden clustering and UMAP

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| kNN | Neighbour count | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| kNN | Search | Fixed | `BPCells::knn_hnsw()` with cosine metric and `ef = 500` | `cluster_embedding_matrix_BPCells()` |
| SNN | Construction and pruning | Fixed | package default: Jaccard weights, `min_val = 1/15`, no self loops | `cluster_knn_snn_leiden()`, `BPCells::knn_to_snn_graph()` |
| Leiden | Resolution | Configurable |  | [`aggregation_GEX_cluster_res`](../parameters.html#aggregation_GEX_cluster_res), [`aggregation_ATAC_cluster_res`](../parameters.html#aggregation_ATAC_cluster_res), [`aggregation_WNN_cluster_res`](../parameters.html#aggregation_WNN_cluster_res) |
| Leiden | Objective, iterations and seed | Fixed | modularity with SNN edge weights, seed 1; package defaults: 2 iterations, `beta = 0.01` | `cluster_knn_snn_leiden()`, `igraph::cluster_leiden()` |
| UMAP | Dimensions | Configurable |  | [`aggregation_UMAP_GEX_PCs`](../parameters.html#aggregation_UMAP_GEX_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |
| UMAP | Neighbours and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs), [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Fixed arguments | Fixed | cosine metric, 2 components, seed 1; package defaults otherwise, including spectral initialisation | `run_UMAP_from_embedding_matrix()`, `uwot::umap()` |
| UMAP sweeps | Grids | Fixed | 3 dimension counts from 5 to the number of data dimensions; 3 neighbour counts from 10 to the configured UMAP neighbour count | `UMAP_n_dims_seq.GEX`, `UMAP_n_dims_seq.ATAC`, `UMAP_neighbors_seq` |

### Weighted nearest neighbours

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Inputs | Embeddings after optional Harmony, and their dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| Candidates | Candidates per modality | Fixed | `candidate_k = 200`, from `BPCells::knn_hnsw()` on the L2-normalized embeddings with Euclidean metric, `k = candidate_k + 1` and `ef = 500` | `WNN_results_raw`, `WNN_results`, `weighted_nearest_neighbors_BPCells()` |
| Final neighbours | Neighbour count, also the bandwidth and imputation neighbourhood | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| Bandwidth | Small-SNN kernel width | Fixed | mean distance, minus the nearest non-self distance, to the k cells with the fewest shared neighbours among those sharing at least one; `sd_scale = 1`; floored at machine epsilon | `weighted_nearest_neighbors_BPCells()`, `calculate_small_SNN_bandwidth()` |
| Weights | Modality weight kernel | Fixed | `exp(−d/σ)`; ratio `within / (cross + 1e-4)` clipped to \[0, 200\]; softmax across modalities | `weighted_nearest_neighbors_BPCells()` |
| Selection | Weighted score and distance | Fixed | `Σ w_m · exp(−d_m/σ_m)` (`kernel_power = 1`); `nn_dist = sqrt((1 − score)/2)` | `weighted_nearest_neighbors_BPCells()` |
| SNN and Leiden | Graph and clustering | Fixed | as for GEX and ATAC: package-default SNN, Leiden modularity, seed 1 | `cluster_WNN_graph()`, `cluster_knn_snn_leiden()` |
| UMAP | Neighbours (capped at the final neighbour count) and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs), [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Fixed arguments | Fixed | precomputed WNN neighbours, 2 components, seed 1; package defaults otherwise | `run_WNN_UMAP()`, `uwot::umap()` |

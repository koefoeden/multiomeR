# GEX, ATAC, batch correction and WNN

This chapter describes normalization and dimensional reduction of both modalities, peak definition, batch correction, weighted nearest-neighbour (WNN) integration, and the shared graph, clustering and UMAP steps. It lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [primary-module graph](implementation_main.md). Library versions are pinned by the Pixi environment: BPCells 0.3.1, igraph 2.3.0, harmony 2.0.2, uwot 0.2.4 and Seurat 5.5.0 at the time of writing.

## GEX normalization and PCA

The combined GEX matrix keeps every gene. At PCA, genes whose total count over the Cell Ranger-kept cells is at or below the minimum are excluded, for both backends. The backend is configurable. The BPCells-native branch computes Pearson residuals with a per-gene method-of-moments theta, clips them, optionally regresses configured cell-level covariates, keeps the genes with the highest residual variance and takes a truncated SVD. The Seurat branch runs SCTransform v2 on the same cells and genes, keeps the same number of variable features and computes the PCA from the dense residual matrix. Regression of the cell-cycle difference is the manifest default; when a cell-cycle column is regressed, S and G2M scores are computed from log-normalized counts with the Seurat 2019 gene sets.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Gene filter | Minimum total count per gene, both backends | Hardcoded: helper default | `min_feature_count = 50`, genes with more than 50 counts kept | `run_GEX_PCA_BPCells()` in `R/processing_GEX_helpers.R` |
| Backend | Normalization and PCA backend | Configurable |  | [`aggregation_GEX_PCA_backend`](../parameters.html#aggregation_GEX_PCA_backend) |
| Regression | Cell-level covariates regressed | Configurable |  | [`aggregation_SCT_regress_vars`](../parameters.html#aggregation_SCT_regress_vars) |
| Components | Number of PCs computed | Configurable |  | last element of [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs) |
| Variable genes | Genes retained by residual variance | Hardcoded: helper default | `n_variable_features = 3000` | `R/processing_GEX_helpers.R` |
| BPCells branch | Residual clip range | Hardcoded: helper default | `clip_range = c(-10, 10)` | `R/processing_GEX_helpers.R` |
| BPCells branch | Minimum variance passed to `sctransform_pearson` | Hardcoded: helper default | `min_var = 0` | `R/processing_GEX_helpers.R` |
| BPCells branch | Theta estimate | Hardcoded: inline literal | method of moments, clamped to \[1e-6, 1e6\] | `R/processing_GEX_helpers.R` |
| BPCells branch | Regression | Hardcoded: inline literal | `BPCells::regress_out(prediction_axis = "row")` on residuals | `R/processing_GEX_helpers.R` |
| BPCells branch | SVD | Hardcoded: inline literal | `BPCells::svds()`, embeddings = right singular vectors × singular values, no centring | `R/processing_GEX_helpers.R` |
| Seurat branch | SCTransform arguments | Hardcoded: inline literal | `conserve.memory = TRUE`, `do.correct.umi = FALSE`, otherwise Seurat 5.5.0 defaults (v2, 5,000 model cells) | `run_Seurat_SCT_for_PCA()` in `R/processing_GEX_helpers.R` |
| Seurat branch | PCA | Hardcoded: inline literal | eigendecomposition of the residual gram matrix | `R/processing_GEX_helpers.R` |
| Cell cycle | Gene sets and scoring | Hardcoded: inline literal | `Seurat::cc.genes.updated.2019`, both organisms; log-normalization with scale factor 10,000; `CC.Difference = S − G2M`; seed 1 | `R/processing_GEX_helpers.R` |

## ATAC peak calling and consensus peaks

Fragments are restricted to the standard chromosomes. Peak-calling groups are defined by a configurable metadata column of the GEX-retained nuclei; groups above the discovery cap are randomly downsampled and fragments are exported per group. MACS3 is run with ATAC-style shift and extension and summit calling, or the BPCells tile caller is used. Summits are extended to fixed-width peaks, peaks overlapping the reference blacklist are removed, and an ArchR-style iterative overlap removal produces one non-overlapping consensus set: within a group, overlapping peaks are ranked by summit significance; across groups, by MACS3 fold enrichment. The peak matrix counts insertions or fragment overlaps, as configured, for the GEX-retained nuclei.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Fragments | Chromosomes retained | Hardcoded: target literal | autosomes, X and Y | `extra_targets/ATAC_targets.R` |
| Groups | Peak-calling grouping column | Configurable |  | [`aggregation_call_peaks_by_cluster_col`](../parameters.html#aggregation_call_peaks_by_cluster_col) |
| Groups | Discovery cap per group | Hardcoded: helper default | `max_cells_per_cluster = 50000`, seed 1 + group index | `R/processing_ATAC_helpers.R` |
| Caller | Peak-calling method | Configurable |  | [`aggregation_ATAC_peak_calling_method`](../parameters.html#aggregation_ATAC_peak_calling_method) |
| MACS3 | Command-line arguments | Hardcoded: inline literal | `-f BED --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all`; default q = 0.05 | `R/processing_ATAC_helpers.R` |
| MACS3 | Effective genome size | Hardcoded: inline literal | GRCh38 2.913e9; mm10 and GRCm39 2.65e9 | `R/processing_ATAC_helpers.R` |
| Tile caller | BPCells settings | Hardcoded: helper default | `peak_width = 500`, `peak_tiling = 3`, `fdr_cutoff = 0.01`, `merge_peaks = "none"` | `R/processing_ATAC_helpers.R` |
| Peak shape | Summit extension | Hardcoded: helper default | `extend_summits = 250`, giving 500 bp peaks | `R/processing_ATAC_helpers.R` |
| Blacklist | Source per genome | Hardcoded: inline literal | GRCh38 Kundaje unified; mm10 Boyle v2; `resources/mm39.excluderanges.bed` | `R/processing_ATAC_helpers.R` |
| Blacklist | Rule | Hardcoded: inline literal | any overlap removes the peak | `R/processing_ATAC_helpers.R` |
| Consensus | Within-group ranking | Hardcoded: target literal | `neg_log10pvalue_summit`, decreasing | `extra_targets/ATAC_targets.R` |
| Consensus | Cross-group ranking | Hardcoded: target literal | `fold_change`, decreasing | `extra_targets/ATAC_targets.R` |
| Consensus | Overlap removal | Hardcoded: inline literal | iterative: reduce, keep best per cluster, drop overlaps, repeat | `R/processing_ATAC_helpers.R` |
| Peak matrix | Counting mode, also used for blacklist QC counts | Configurable |  | [`aggregation_ATAC_peak_matrix_mode`](../parameters.html#aggregation_ATAC_peak_matrix_mode) |
| Peak matrix | Cells | Hardcoded: target literal | post-doublet-filter GEX metadata barcodes | `extra_targets/ATAC_targets.R` |

## ATAC TF-IDF and LSI

The QC-filtered peak matrix is transformed with Signac's TF-IDF method 1 and decomposed with BPCells SVD to obtain latent semantic indexing (LSI) embeddings. No code rule drops the first component; the manifest defaults for the data and UMAP dimensions start at the second component.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| TF-IDF | Method and scale factor | Hardcoded: helper default | term frequency × inverse document frequency, `log1p(scale_factor × TF·IDF)`, `scale_factor = 10000` | `R/processing_ATAC_helpers.R` |
| SVD | Components computed | Configurable |  | last element of [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| SVD | Function | Hardcoded: inline literal | `BPCells::svds()`, embeddings = right singular vectors × singular values | `R/processing_ATAC_helpers.R` |
| Dimensions | Components used downstream | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |

## Harmony batch correction

Harmony can be applied separately to the selected GEX PCs and ATAC LSI dimensions. When no covariates are configured, the uncorrected embedding is returned unchanged. Several covariates are collapsed into one interaction batch factor, and nuclei with missing covariate values are dropped with a warning.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Covariates | Shared correction columns | Configurable |  | [`aggregation_harmony_correction_metadata_col_names`](../parameters.html#aggregation_harmony_correction_metadata_col_names) |
| Covariates | Additional ATAC covariates | Configurable |  | [`aggregation_extra_harmony_covars_ATAC`](../parameters.html#aggregation_extra_harmony_covars_ATAC) |
| Dimensions | Corrected dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| Function | Call | Hardcoded: inline literal | `harmony::RunHarmony(max_iter = 25, lambda = 1)`, other arguments harmony 2.0.2 defaults | `R/processing_ATAC_helpers.R` |
| Covariates | Combination rule | Hardcoded: inline literal | interaction of all covariates as one batch factor | `R/processing_ATAC_helpers.R` |
| Covariates | Missing values | Hardcoded: inline literal | nuclei removed with a warning | `R/processing_ATAC_helpers.R` |
| Resources | Cores | Hardcoded: target literal | 6 | `extra_targets/GEX_merge_and_dim_reduc_targets.R`, `extra_targets/ATAC_targets.R` |

## Graph construction, Leiden clustering and UMAP

For GEX and ATAC, approximate nearest neighbours are found with BPCells HNSW on the selected dimensions, converted to a shared-nearest-neighbour (SNN) graph with Jaccard weights and clustered with Leiden. Because the query cell is its own first neighbour, a neighbour count of k yields k − 1 non-self neighbours. Clusters are renumbered by size. UMAP is computed with uwot on the selected dimensions. QC parameter sweeps additionally render UMAPs over grids of dimensions and neighbour counts.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| kNN | Neighbour count | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| kNN | Function, metric and search effort | Hardcoded: inline literal | `BPCells::knn_hnsw(metric = "cosine", ef = 500)` | `R/processing_ATAC_helpers.R` |
| SNN | Construction and pruning | Hardcoded: environment pin | `BPCells::knn_to_snn_graph()` defaults: Jaccard weights, `min_val = 1/15`, no self loops | `R/processing_ATAC_helpers.R` |
| Leiden | Resolution | Configurable |  | [`aggregation_GEX_cluster_res`](../parameters.html#aggregation_GEX_cluster_res), [`aggregation_ATAC_cluster_res`](../parameters.html#aggregation_ATAC_cluster_res), [`aggregation_WNN_cluster_res`](../parameters.html#aggregation_WNN_cluster_res) |
| Leiden | Objective, weights and seed | Hardcoded: inline literal | `igraph::cluster_leiden(objective_function = "modularity")` with SNN weights; seed 1 | `R/processing_ATAC_helpers.R` |
| Leiden | Iterations and randomness | Hardcoded: environment pin | igraph 2.3.0 defaults: 2 iterations, beta 0.01 | `R/processing_ATAC_helpers.R` |
| Leiden | Relabelling | Hardcoded: inline literal | renumbered by size, largest first | `R/processing_ATAC_helpers.R` |
| UMAP | Dimensions | Configurable |  | [`aggregation_UMAP_GEX_PCs`](../parameters.html#aggregation_UMAP_GEX_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |
| UMAP | Neighbours and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs), [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Function and fixed arguments | Hardcoded: inline literal | `uwot::umap(metric = "cosine", n_components = 2, n_sgd_threads = 0)`, seed 1 | `R/processing_ATAC_helpers.R` |
| UMAP | Other arguments | Hardcoded: environment pin | uwot 0.2.4 defaults, spectral initialisation | `R/processing_ATAC_helpers.R` |
| UMAP sweeps | Grids | Hardcoded: target literal | dimensions from 5 to the configured count in 3 steps; neighbours from 10 to the configured UMAP count in 3 steps | `extra_targets/GEX_graph_and_cluster_targets.R`, `extra_targets/ATAC_targets.R` |
| Resources | Threads | Hardcoded: target literal | 6 | `extra_targets/GEX_graph_and_cluster_targets.R`, `extra_targets/ATAC_targets.R` |

## Weighted nearest neighbours

The BPCells-native WNN implementation is described, with its deviations from Seurat and its validation, in [Algorithmic implementations](algorithm_validation.md#native-weighted-nearest-neighbors). The selected GEX and ATAC dimensions after optional Harmony correction are aligned by barcode and L2-normalized. HNSW finds a large candidate set per modality; Seurat's small-SNN bandwidth strategy sets the kernel width per cell from the configured neighbour count; per-cell modality weights come from within- versus cross-modality prediction kernels; and the union of candidates is ranked by the weighted kernel score to select the final neighbours. BPCells builds the SNN graph, which is clustered with Leiden and embedded with UMAP using the precomputed neighbours.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Inputs | Embeddings and dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) after Harmony |
| Inputs | Normalization | Hardcoded: helper default | row-wise L2 | `weighted_nearest_neighbors_BPCells()` in `R/processing_multimodal_helpers.R` |
| Candidates | Candidates per modality | Hardcoded: target literal | `candidate_k = 200` | `extra_targets/WNN_targets.R` |
| Candidates | Per-modality search | Hardcoded: inline literal | `BPCells::knn_hnsw(k = candidate_k + 1, metric = "euclidean", ef = 500)` | `R/processing_multimodal_helpers.R` |
| Final neighbours | Neighbour count, also the bandwidth and imputation neighbourhood | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| Bandwidth | Small-SNN kernel width | Hardcoded: inline literal | mean distance to the k lowest-shared-neighbour cells after subtracting the nearest non-self distance; `sd_scale = 1`; floor at machine epsilon | `src/wnn_snn_bandwidth.cpp`, `R/processing_multimodal_helpers.R` |
| Weights | Modality weight kernel | Hardcoded: inline literal | `exp(−d/σ)`; ratio `within / (cross + 1e-4)` clipped to \[0, 200\]; softmax-normalized | `R/processing_multimodal_helpers.R` |
| Selection | Weighted score and distance | Hardcoded: inline literal | `Σ w_m · exp(−(d_m/σ_m))`, `kernel_power = 1`; `nn_dist = sqrt((1 − score)/2)` | `R/processing_multimodal_helpers.R` |
| SNN and Leiden | Graph and clustering | Hardcoded: environment pin | `knn_to_snn_graph(min_val = 1/15)`; Leiden modularity, seed 1 | `R/processing_multimodal_helpers.R` |
| UMAP | Neighbours and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs) capped at the neighbour count; [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Fixed arguments | Hardcoded: inline literal | precomputed neighbours, 2 components, seed 1 | `R/processing_multimodal_helpers.R` |
| Resources | Threads | Hardcoded: target literal | 6 | `extra_targets/WNN_targets.R` |
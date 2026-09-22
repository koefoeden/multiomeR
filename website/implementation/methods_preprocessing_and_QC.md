# Preprocessing and nucleus QC

This chapter describes per-GEM-well preprocessing and the successive nucleus filters up to the final WNN cell set, and lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [main pipeline graph](implementation_main.md). GEM-well-level settings are columns of `cfg_GEM_wells.tsv`, described in [GEM well table](../reference_GEM_wells.html); aggregation-level settings are manifest parameters.

## Aggregation inputs and operational settings

These manifest parameters select inputs, plot variables and execution behaviour rather than algorithm settings. They are listed here so that every manifest row has a home in this book.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Inputs | GEM wells combined | Configurable |  | [`aggregation_GEM_well_IDs`](../parameters.html#aggregation_GEM_well_IDs) |
| Inputs | Donor metadata table | Configurable |  | [`aggregation_donor_id_metadata_tsv`](../parameters.html#aggregation_donor_id_metadata_tsv) |
| Inputs | Aggregation active | Configurable |  | [`is_active`](../parameters.html#is_active) |
| Inputs | Optional modules enabled | Configurable |  | [`modules`](../parameters.html#modules) |
| Plots | Categorical and continuous metadata plotted | Configurable |  | [`aggregation_categorical_vars`](../parameters.html#aggregation_categorical_vars), [`aggregation_continuous_vars`](../parameters.html#aggregation_continuous_vars) |
| Plots | Additional genes plotted | Configurable |  | [`aggregation_other_interesting_genes`](../parameters.html#aggregation_other_interesting_genes) |
| Tracks | Roadmap epigenome tracks | Configurable |  | [`aggregation_roadmap_EDACC_names`](../parameters.html#aggregation_roadmap_EDACC_names) |
| Execution | Targets skipped by `tar_make()` helpers | Configurable |  | [`aggregation_tar_make_skip_regex_patterns`](../parameters.html#aggregation_tar_make_skip_regex_patterns) |

## Per-GEM-well inputs and metrics

Each GEM well supplies a Cell Ranger ARC count directory. The GEX matrix is imported into BPCells from the filtered feature-barcode matrix, or from a CellBender output when configured, keeping only gene-expression features. ATAC fragments are imported from the Cell Ranger fragment file. Both carry a GEM-well-prefixed barcode so nuclei stay distinct across wells. The called-cell universe is the set of barcodes that Cell Ranger flagged as cells in its per-barcode metrics. ATAC QC metrics come from BPCells, GEX metrics from the imported matrix. Optional genotype demultiplexing runs cellsnp-lite on the ATAC BAM and Vireo against the configured VCF; the resulting donor label replaces the GEM-well donor identifier, and doublet or unassigned calls are recorded as metadata only. The BPCells-native AMULET implementation, described in [Algorithmic implementations](algorithm_validation.md#bpcells-native-amulet), calculates overlap metrics and q-values on the called cells; no automatic AMULET filter is applied.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| GEX input | CellBender versus Cell Ranger matrix | Configurable |  | `GEM_well_add_cellbender`, `GEM_well_cellbender_h5_file` in `cfg_GEM_wells.tsv` |
| GEX input | Feature type kept | Hardcoded: target literal | `Gene Expression` | `GEX_counts_BPcells_matrix_dir` in `extra_targets/per_GEM_well_targets.R` |
| Barcodes | Prefix on GEX and fragment barcodes | Hardcoded: target literal | `<GEM_well_ID>_` | `extra_targets/per_GEM_well_targets.R` |
| Barcodes | Called-cell universe | Hardcoded: target literal | rows of `per_barcode_metrics.csv` with `is_cell == 1` | `cellranger_kept_metadata_tibble` in `extra_targets/per_GEM_well_targets.R` |
| ATAC metrics | QC function and blacklist | Hardcoded: target literal | `BPCells::qc_scATAC()` with Ensembl genes and an empty blacklist | `ATAC_qc_metrics_tibble` in `extra_targets/per_GEM_well_targets.R` |
| ATAC metrics | TSS enrichment window | Hardcoded: environment pin | BPCells 0.3.1: 101 bp centre window, 100 bp flanks at ±1.9–2 kb, denominator floor 0.1 | BPCells `qc_scATAC()` |
| ATAC metrics | Nucleosome signal | Hardcoded: target literal | mono-nucleosomal / sub-nucleosomal fragment counts | `extra_targets/per_GEM_well_targets.R` |
| GEX metrics | Mitochondrial gene pattern | Hardcoded: target literal | `(?i)^MT-` | `extra_targets/per_GEM_well_targets.R` |
| Demultiplexing | VCF, donor count and donor label | Configurable |  | `GEM_well_donors_VCF_file`, `GEM_well_n_donors`, `GEM_well_donor_id` in `cfg_GEM_wells.tsv` |
| Demultiplexing | cellsnp-lite settings | Hardcoded: helper default | `--minMAF 0.1`, `--minCOUNT 20`, `--UMItag None`, ATAC BAM | `R/parallel_GEM_well_preprocessing_helpers.R` |
| Demultiplexing | Vireo settings | Hardcoded: helper default | genotypes from the VCF (`-t GT`), no genotype learning | `R/parallel_GEM_well_preprocessing_helpers.R` |
| Demultiplexing | Cores | Hardcoded: target literal | cellsnp 6, Vireo 4 | `extra_targets/per_GEM_well_targets.R` |
| Demultiplexing | Use of doublet and unassigned calls | Hardcoded: inline literal | metadata only; removable through the GEM-well exclusion list | `extra_targets/per_GEM_well_targets.R` |
| AMULET | Barcodes scored | Hardcoded: target literal | Cell Ranger called cells | `amulet_metrics_tibble` in `extra_targets/per_GEM_well_targets.R` |
| AMULET | Maximum fragment size | Hardcoded: helper default | 1,000 bp | `calculate_amulet_metrics_BPCells()` in `R/amulet_BPCells_helpers.R` |
| AMULET | Excluded regions | Hardcoded: helper default | chrM, chrX, chrY and their aliases | `R/amulet_BPCells_helpers.R` |
| AMULET | Cell Ranger end-inclusive shift | Hardcoded: target literal | end − 1 | `extra_targets/per_GEM_well_targets.R` |
| AMULET | High-overlap-site removal | Hardcoded: helper default | on; Poisson P \< 0.01 | `remove_high_overlap_amulet_loci()` in `R/amulet_BPCells_helpers.R` |
| AMULET | Per-cell test | Hardcoded: inline literal | upper-tail Poisson on loci covered by more than two fragments; BH q-values | `R/amulet_BPCells_helpers.R` |
| AMULET | Minimum fragments argument | Hardcoded: target literal | 1,000, inert because barcodes are supplied | `extra_targets/per_GEM_well_targets.R` |
| AMULET | Automatic filter | Hardcoded: inline literal | none; q-values usable in the GEM-well exclusion list | `R/processing_and_aggregation_constants.R` |

## GEM-well exclusions and the GEX cell universe

Configured GEM-well exclusions are dplyr filter expressions evaluated on the called cells of each well; matching nuclei are removed from the union of called cells before any aggregation-level analysis. The per-well GEX matrices are then column-bound without further filtering. At PCA, the metadata are intersected with the matrix barcodes, so only barcodes present in the aligned count data continue.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Exclusions | Per-well exclusion expressions | Configurable |  | `GEM_well_QC_exclude_list` in `cfg_GEM_wells.tsv` |
| Exclusions | Where applied | Hardcoded: target literal | union of called cells, before GEX combination | `GEX_cellranger_kept_metadata_tibble` in `extra_targets/general_aggregation_targets.R` |
| Combination | Matrix merge | Hardcoded: target literal | column bind of every per-well matrix, no gene or cell filter | `extra_targets/GEX_merge_and_dim_reduc_targets.R` |
| Retention | Barcodes carried into PCA | Hardcoded: inline literal | intersection of metadata and matrix barcodes | `run_GEX_PCA_BPCells()` in `R/processing_GEX_helpers.R` |

## Cluster-size filter and GEX doublets

After GEX graph construction and Leiden clustering, clusters below the configured minimum size are removed; the same threshold is applied to the ATAC and WNN clusters and to the optional subgroup clusterings. The filter keeps clusters at or above the threshold and is disabled when the value is missing or at most one. scDblFinder is then run per GEM well on the raw GEX counts of the retained nuclei. Its cluster labels are the annotation-derived scDblFinder groups: the cell-type label when the cluster was assigned, otherwise the status and cluster identifier, optionally collapsed with the configured map. Doublet removal is controlled separately at the nucleus level and at the cluster level; the cluster-level rule drops a Leiden cluster whose doublet fraction exceeds the threshold.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cluster filter | Minimum cluster size, applied to GEX, ATAC, WNN and subgroups | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |
| Cluster filter | Comparison and disabling rule | Hardcoded: inline literal | keep size ≥ threshold; disabled when NULL, NA or ≤ 1 | `cluster_embedding_matrix_BPCells()` in `R/processing_ATAC_helpers.R` |
| GEX scDblFinder | Input | Hardcoded: target literal | per-well slice of the raw aggregated GEX counts | `extra_targets/GEX_graph_and_cluster_targets.R` |
| GEX scDblFinder | Cluster labels | Hardcoded: inline literal | annotation-derived scDblFinder groups | `R/cluster_annotation_helpers.R`, `R/processing_GEX_helpers.R` |
| GEX scDblFinder | Label collapse map, reused for ATAC | Configurable |  | [`aggregation_scDblFinder_GEX_cell_type_collapse_list`](../parameters.html#aggregation_scDblFinder_GEX_cell_type_collapse_list) |
| GEX scDblFinder | `dbr.sd` | Hardcoded: target literal | 1.0 | `extra_targets/GEX_graph_and_cluster_targets.R` |
| GEX scDblFinder | Return type and parallelism | Hardcoded: inline literal | scores; `BiocParallel::SerialParam()` | `R/processing_GEX_helpers.R` |
| GEX scDblFinder | Other arguments | Hardcoded: environment pin | scDblFinder 1.24.0 defaults | `extra_targets/GEX_graph_and_cluster_targets.R` |
| GEX doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_GEX_remove_called_doublets) |
| GEX doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster) |
| GEX doublets | Cluster fraction rule | Hardcoded: inline literal | fraction over all calls on the raw Leiden cluster; strictly greater than the threshold | `R/processing_GEX_helpers.R` |

## ATAC QC and doublets

The ATAC branch starts from the GEX-retained nuclei: peak-calling groups, the consensus peak matrix and peak-level QC metrics are all built on that cell set. Peak-level metrics are evaluated against the configured aggregation-level exclusions, ATAC clusters below the minimum size are removed, and a separate scDblFinder run uses the feature-aggregation adaptation described in [Algorithmic implementations](algorithm_validation.md#bpcells-backed-atac-scdblfinder-feature-aggregation): peaks are grouped by k-means on the aggregation-wide LSI loadings, the groups are summed with BPCells, and the compact matrix is passed to the unchanged classifier.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cell set | Nuclei entering the ATAC branch | Hardcoded: target literal | post-doublet-filter GEX metadata | `extra_targets/ATAC_targets.R` |
| Peak QC | Metrics | Hardcoded: inline literal | ATAC counts in peaks, blacklist counts and fraction, peak-count fraction of fragments, peak-count enrichment | `R/processing_ATAC_helpers.R` |
| Peak QC | Exclusion expressions | Configurable |  | [`aggregation_QC_exclude_list_combined_object`](../parameters.html#aggregation_QC_exclude_list_combined_object) |
| ATAC scDblFinder | Nuclei scored | Hardcoded: target literal | post-cluster-filter ATAC cells passing peak QC | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | LSI dimensions for feature groups | Hardcoded: target literal | `intersect(2:20, aggregation_ATAC_data_PCs)` | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | Number of feature groups and seed | Hardcoded: target literal | 50 groups; seed 1 | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | Grouping method | Hardcoded: inline literal | `stats::kmeans(iter.max = 50, nstart = 1)` on loadings; sums via `BPCells::pseudobulk_matrix(method = "sum")` | `get_feature_groups_from_LSI_loadings()`, `aggregate_BPCells_rows_by_group()` in `R/processing_GEX_helpers.R` |
| ATAC scDblFinder | Classifier arguments | Hardcoded: target literal | `dbr.sd = 1.0`, `aggregateFeatures = FALSE`, `nfeatures = 50`, `processing = "normFeatures"` | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | Cluster labels | Hardcoded: inline literal | ATAC annotation-derived groups collapsed with the GEX map | `extra_targets/ATAC_targets.R` |
| ATAC doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_ATAC_remove_called_doublets) |
| ATAC doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster) |

## WNN cell set

WNN integration uses the nuclei retained by the ATAC branch that have rows in both corrected embeddings. The joint clusters are filtered with the minimum-cluster-size rule; when nuclei are dropped, the WNN graph and UMAP are recomputed on the retained nuclei while the cluster labels from the first run are kept. The GEX, ATAC and WNN metadata therefore describe successive cell universes.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cell set | Nuclei entering WNN | Hardcoded: target literal | post-doublet-filter ATAC metadata ∩ rows of both embeddings | `extra_targets/WNN_targets.R`, `R/processing_multimodal_helpers.R` |
| Cluster filter | Threshold and recomputation | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |
| Cluster filter | Recompute rule | Hardcoded: target literal | graph and UMAP recomputed when any nucleus is dropped; labels kept | `extra_targets/WNN_targets.R` |
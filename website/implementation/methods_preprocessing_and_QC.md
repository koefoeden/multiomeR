# Preprocessing and nucleus QC

This chapter covers per-GEM-well preprocessing and the successive nucleus filters up to the final WNN cell set. The target structure is shown in the [primary-module graph](implementation_main.md), and GEM-well settings are columns of `cfg_GEM_wells.tsv`, described in [GEM well table](../reference_GEM_wells.html). The AMULET and ATAC scDblFinder adaptations are compared with their references in [Algorithmic implementations](algorithm_validation.md). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

{{< include _shared_methods/quality_control.md >}}

## Aggregation inputs and operational settings

These manifest parameters select inputs, plot variables and execution behaviour rather than algorithm settings.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Inputs | GEM wells combined | Configurable |  | [`aggregation_GEM_well_IDs`](../parameters.html#aggregation_GEM_well_IDs) |
| Inputs | Donor metadata table | Configurable |  | [`aggregation_donor_id_metadata_tsv`](../parameters.html#aggregation_donor_id_metadata_tsv) |
| Inputs | Aggregation active | Configurable |  | [`is_active`](../parameters.html#is_active) |
| Inputs | Optional modules enabled | Configurable |  | [`modules`](../parameters.html#modules) |
| Plots | Categorical and continuous metadata plotted | Configurable |  | [`aggregation_categorical_vars`](../parameters.html#aggregation_categorical_vars), [`aggregation_continuous_vars`](../parameters.html#aggregation_continuous_vars) |
| Plots | Additional genes plotted | Configurable |  | [`aggregation_other_interesting_genes`](../parameters.html#aggregation_other_interesting_genes) |
| Tracks | Roadmap epigenome tracks | Configurable |  | [`aggregation_roadmap_EDACC_names`](../parameters.html#aggregation_roadmap_EDACC_names) |

## Per-GEM-well inputs and metrics

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| GEX input | CellBender versus Cell Ranger matrix | Configurable |  | `GEM_well_add_cellbender`, `GEM_well_cellbender_h5_file` |
| GEX input | Feature type kept | Fixed | `Gene Expression` | `GEX_counts_BPcells_matrix_dir` |
| Barcodes | Called-cell universe | Fixed | barcodes with `is_cell == 1` in the Cell Ranger `per_barcode_metrics.csv` | `per_barcode_metrics_tibble` |
| ATAC metrics | TSS enrichment | Fixed | `BPCells::qc_scATAC()` on Ensembl gene TSSs with an empty blacklist; package default: 101 bp centre window, 100 bp flanks 1.9–2 kb up- and downstream, flank signal floored at 0.1 | `ATAC_qc_metrics_tibble`, `BPCells::qc_scATAC()` |
| ATAC metrics | Nucleosome signal | Fixed | mono-nucleosomal / sub-nucleosomal fragment counts | `ATAC_qc_metrics_tibble` |
| GEX metrics | Mitochondrial genes | Fixed | gene names matching `(?i)^MT-` | `GEX_basic_metadata_tibble` |
| Demultiplexing | VCF, donor count and donor label | Configurable |  | `GEM_well_donors_VCF_file`, `GEM_well_n_donors`, `GEM_well_donor_id` |
| Demultiplexing | cellsnp-lite | Fixed | ATAC BAM, Cell Ranger-called barcodes, `--minMAF 0.1`, `--minCOUNT 20`, `--UMItag None` | `cellsnp_dir`, `get_cellsnp_dir()` |
| Demultiplexing | Vireo | Fixed | donor genotypes from the VCF (`-t GT`); no genotype learning | `get_vireo_donor_ids_tibble()` |
| AMULET | Nuclei scored | Fixed | Cell Ranger-called barcodes; no fragment minimum | `amulet_metrics_tibble` |
| AMULET | Fragment filters | Fixed | fragments ≤ 1,000 bp; chrM, chrX, chrY and their aliases excluded; Cell Ranger end coordinate shifted by −1 | `calculate_amulet_metrics_BPCells()`, `get_amulet_fragment_overlaps_BPCells()` |
| AMULET | High-overlap-site removal | Fixed | on; loci with Poisson P \< 0.01 across nuclei removed | `calculate_amulet_metrics_BPCells()`, `remove_high_overlap_amulet_loci()` |
| AMULET | Per-nucleus test | Fixed | upper-tail Poisson on the number of loci covered by more than two fragments; BH q-values | `calculate_amulet_metrics_BPCells()` |

## GEM-well exclusions and the GEX cell universe

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Exclusions | Per-well exclusion expressions | Configurable |  | `GEM_well_QC_exclude_list` |
| Cell set | Nuclei entering GEX PCA | Fixed | Cell Ranger-called nuclei passing the per-well exclusions, intersected with the combined matrix barcodes | `GEX_cellranger_kept_metadata_tibble`, `run_GEX_PCA_BPCells()` |

## Cluster-size filter and GEX doublets

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cluster filter | Minimum cluster size, applied to GEX, ATAC, WNN and subgroups | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |
| Cluster filter | Comparison and disabling rule | Fixed | keep clusters with size ≥ threshold; disabled when NULL, NA or ≤ 1 | `filter_clusters_by_min_barcodes()` |
| GEX scDblFinder | Nuclei and counts | Fixed | GEX nuclei after the cluster-size filter and before doublet removal; raw counts; one run per GEM well | `scDblFinder_GEM_well_tibble.GEX`, `scDblFinder_results_by_GEM_well_tibble.GEX` |
| GEX scDblFinder | Label collapse map, reused for ATAC | Configurable |  | [`aggregation_scDblFinder_GEX_cell_type_collapse_list`](../parameters.html#aggregation_scDblFinder_GEX_cell_type_collapse_list) |
| GEX scDblFinder | Classifier arguments | Fixed | `dbr.sd = 1.0`; package defaults otherwise | `scDblFinder_results_by_GEM_well_tibble.GEX`, `scDblFinder::scDblFinder()` |
| GEX doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_GEX_remove_called_doublets) |
| GEX doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster) |
| GEX doublets | Cluster fraction rule, also used for ATAC | Fixed | called-doublet fraction of each Leiden cluster, computed before nucleus-level removal; cluster removed when the fraction is strictly greater than the threshold | `filter_metadata_by_scDblFinder()` |

## ATAC QC and doublets

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cell set | Nuclei entering the ATAC branch | Fixed | GEX nuclei after doublet removal | `BCs_per_peak_cluster_list.ATAC`, `consensus_peak_BPCells_matrix_dir.ATAC`, `metadata_w_QC_tibble.ATAC` |
| Peak QC | Metrics | Fixed | `nCount_ATAC` (counts in consensus peaks); blacklist counts / `nCount_ATAC`; `nCount_ATAC` / Cell Ranger ATAC fragments; that fraction divided by the genome fraction covered by peaks | `get_ATAC_QC_metadata_from_BPCells()` |
| Peak QC | Exclusion expressions | Configurable |  | [`aggregation_QC_exclude_list_combined_object`](../parameters.html#aggregation_QC_exclude_list_combined_object) |
| ATAC scDblFinder | Nuclei scored | Fixed | ATAC nuclei after the cluster-size filter that pass peak QC | `scDblFinder_GEM_well_tibble.ATAC` |
| ATAC scDblFinder | Feature groups | Fixed | 50 groups from `stats::kmeans(iter.max = 50, nstart = 1)`, seed 1, on the peak loadings of LSI dimensions `intersect(2:20, aggregation_ATAC_data_PCs)` | `scDblFinder_feature_groups.ATAC`, `get_feature_groups_from_LSI_loadings()` |
| ATAC scDblFinder | Classifier arguments | Fixed | `dbr.sd = 1.0`, `aggregateFeatures = FALSE`, `nfeatures = 50`, `processing = "normFeatures"`; package defaults otherwise | `scDblFinder_results_by_GEM_well_tibble.ATAC`, `scDblFinder::scDblFinder()` |
| ATAC scDblFinder | Cluster labels | Fixed | ATAC annotation-derived scDblFinder groups, collapsed with the GEX map | `scDblFinder_GEM_well_tibble.ATAC` |
| ATAC doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_ATAC_remove_called_doublets) |
| ATAC doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster) |

## WNN cell set

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cell set | Nuclei entering WNN | Fixed | ATAC nuclei after doublet removal that have rows in both corrected embeddings | `embedding_matrices.WNN`, `get_WNN_embedding_matrices()` |
| Cluster filter | Threshold | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |

# Primary-module methods

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

This page describes the analyses of the primary module in the order of the checkpoints in [Run your own analysis](main_running.md). The same descriptions form the Supplementary Methods of the multiomeR manuscript. Configurable settings and their defaults are listed in the [parameter browser](parameters.html#workflow=aggregation), and the reimplemented algorithms are compared with their references in [Algorithm validation](algorithm_validation.md).

Each section ends with its source files and a simplified view of its target graph. The views keep real target names but omit lower-level targets; [Change the pipeline](change_the_pipeline.md#graph-views) explains how they are made.

## Preprocessing {#preprocessing}

**Checkpoint:** [1: Pre-aggregation QC](main_running.md#checkpoint-1)

{{< include _shared_methods/quality_control.md >}}

GEM-well settings are described in [GEM well table](reference_GEM_wells.md), and the AMULET implementation is compared with scDblFinder in [Algorithm validation](algorithm_validation.md#bpcells-native-amulet).

Two checks at checkpoint 1 look at samples rather than nuclei, and neither removes nuclei. The first takes each sample's median TSS enrichment, fraction of fragments in peaks and mitochondrial UMI percentage over its nuclei that pass per-well QC. It flags a sample whose median lies more than three median absolute deviations on the poor side of the other samples in its dataset. The second sums each sample's raw Cell Ranger UMIs of XIST and of all chromosome Y genes, per 10,000 UMIs. It labels a sample that lies nearer the median of the other recorded sex than of its own, on log scales.

Once the GEX clusters are labelled at checkpoint 3, a third check shows which cell types the per-well rules remove. It scores every Cell Ranger-called nucleus on chromatin, because many removed nuclei carry RNA of other cells. A nucleus's share of ATAC fragments in each marker set's gene windows, the gene body and 2 kb upstream, is scaled from the 90th percentile in nuclei with other cluster labels (0) to the median in nuclei with that label (1). The nucleus takes the top set among the cell types labelled in its GEM well when it reaches 0.4, and is otherwise unclassified. The scored types of removed nuclei are then compared with those of kept nuclei, for all rules and for each.

### Source and target graph

**Source:** [`extra_targets/per_GEM_well_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/per_GEM_well_targets.R), [`extra_targets/general_aggregation_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/general_aggregation_targets.R), [`R/parallel_GEM_well_preprocessing_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/parallel_GEM_well_preprocessing_helpers.R), [`R/amulet_BPCells_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/amulet_BPCells_helpers.R), [`R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R), [`R/processing_ATAC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_ATAC_helpers.R), [`R/sample_QC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/sample_QC_helpers.R).

This view covers per-GEM-well processing and QC, including optional ambient RNA correction, donor demultiplexing, doublet detection, barcode filtering, and handoffs into aggregation-level GEX and ATAC objects.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/parallel_v2.mmd")
```

## GEX {#gex}

**Checkpoints:** [2: GEX dimension reduction](main_running.md#checkpoint-2), [3: GEX clusters and cell types](main_running.md#checkpoint-3)

### Normalization and PCA {#gex-normalization}

{{< include _shared_methods/GEX_normalization.md >}}

### Batch correction and clustering {#batch-correction-and-clustering}

{{< include _shared_methods/batch_correction_and_clustering.md >}}

### Cell-type annotation {#cell-type-annotation}

{{< include _shared_methods/cell_type_annotation.md >}}

The BPCells-native UCell scorer is compared with UCell in [Algorithm validation](algorithm_validation.md#bpcells-native-ucell-scoring).

### Doublet removal {#gex-doublets}

{{< include _shared_methods/GEX_doublets.md >}}

### Source and target graph

**Source:** [`extra_targets/GEX_merge_and_dim_reduc_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_merge_and_dim_reduc_targets.R), [`extra_targets/GEX_graph_and_cluster_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_graph_and_cluster_targets.R), [`R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R), [`R/processing_ATAC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_ATAC_helpers.R), [`R/cluster_annotation_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/cluster_annotation_helpers.R), [`src/cluster_UCell_chunk.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/cluster_UCell_chunk.cpp).

This view runs from the aggregated counts through PCA, Harmony correction, clustering and UMAP to cell-type annotation and GEX doublet calls.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/GEX_v2.mmd")
```

## ATAC {#atac}

**Checkpoints:** [4: Peak QC](main_running.md#checkpoint-4), [5: ATAC filtering](main_running.md#checkpoint-5), [6: ATAC dimension reduction](main_running.md#checkpoint-6), [7: ATAC clusters and motifs](main_running.md#checkpoint-7)

### Peaks and LSI {#peaks-and-lsi}

{{< include _shared_methods/ATAC_peaks_and_LSI.md >}}

### Filtering and doublet removal {#atac-filtering}

ATAC clusters are built from the LSI dimensions as described in [Batch correction and clustering](#batch-correction-and-clustering).

{{< include _shared_methods/ATAC_filtering_and_doublets.md >}}

[Algorithm validation](algorithm_validation.md#bpcells-backed-atac-scdblfinder-feature-aggregation) describes how this feature aggregation differs from scDblFinder's own.

### Motif accessibility {#motif-accessibility}

{{< include _shared_methods/motif_accessibility.md >}}

### Gene activity {#gene-activity}

{{< include _shared_methods/ATAC_gene_activity.md >}}

They are computed only when [`aggregation_ATAC_gene_activity_scores`](parameters.html#aggregation_ATAC_gene_activity_scores) is `true`, and they feed `marker_gene_activity_dot_plot.png` at [checkpoint 7](main_running.md#checkpoint-7).

### Source and target graph

**Source:** [`extra_targets/ATAC_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/ATAC_targets.R), [`extra_targets/ATAC_gene_activity_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/ATAC_gene_activity_targets.R), [`R/processing_ATAC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_ATAC_helpers.R), [`R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R), [`R/celltype_labeling_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/celltype_labeling_helpers.R), [`resources/`](https://github.com/koefoeden/multiomeR/tree/main/resources).

This view runs from the combined fragments and the GEX-annotated nuclei through consensus peaks and LSI to ATAC clusters, doublet calls, motif-family accessibility and coverage tracks.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/ATAC_v2.mmd")
```

## WNN {#wnn}

**Checkpoint:** [8: WNN integration](main_running.md#checkpoint-8)

### Weighted nearest neighbors {#weighted-nearest-neighbors}

{{< include _shared_methods/WNN_integration.md >}}

The WNN implementation is compared with Seurat in [Algorithm validation](algorithm_validation.md#native-weighted-nearest-neighbors).

### Final cell set {#final-cell-set}

{{< include _shared_methods/WNN_cell_set.md >}}

### Source and target graph

**Source:** [`extra_targets/WNN_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/WNN_targets.R), [`R/processing_multimodal_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R), [`src/wnn_snn_bandwidth.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/wnn_snn_bandwidth.cpp).

This view combines the Harmony-corrected GEX and ATAC embeddings into the WNN results, clusters, UMAP embeddings and final annotated metadata.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/WNN_v2.mmd")
```

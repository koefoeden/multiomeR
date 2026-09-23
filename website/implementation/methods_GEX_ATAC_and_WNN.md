# GEX, ATAC, batch correction and WNN

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

This page covers normalization and dimensional reduction of both modalities, peak definition, batch correction, graph construction and clustering, and weighted nearest-neighbour (WNN) integration, with the target graph of each stage. Configurable settings and their defaults are listed in the [parameter browser](../parameters.html). The WNN implementation is compared with Seurat in [Algorithmic implementations](algorithm_validation.md#native-weighted-nearest-neighbors).

## GEX normalization and PCA

{{< include _shared_methods/GEX_normalization.md >}}

**Source:** [`extra_targets/GEX_merge_and_dim_reduc_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_merge_and_dim_reduc_targets.R), [`R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R).

The GEX graph also covers clustering, marker detection, cell-type annotation and GEX review outputs.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/GEX_v2.mmd")
```

## ATAC peak definition and dimensional reduction

{{< include _shared_methods/ATAC_peaks_and_LSI.md >}}

**Source:** [`extra_targets/ATAC_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/ATAC_targets.R), [`R/processing_ATAC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_ATAC_helpers.R).

The ATAC graph also covers ATAC QC, chromVAR scoring and coverage tracks.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/ATAC_v2.mmd")
```

## Batch correction, clustering and weighted nearest neighbours

{{< include _shared_methods/batch_correction_clustering_WNN.md >}}

**Source:** [`extra_targets/GEX_graph_and_cluster_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_graph_and_cluster_targets.R), [`extra_targets/WNN_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/WNN_targets.R), [`R/processing_ATAC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_ATAC_helpers.R), [`R/processing_multimodal_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R), [`src/wnn_snn_bandwidth.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/wnn_snn_bandwidth.cpp).

The WNN graph also covers modality weights, cluster comparison, integrated metadata and the Seurat/Signac compatibility export.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/WNN_v2.mmd")
```

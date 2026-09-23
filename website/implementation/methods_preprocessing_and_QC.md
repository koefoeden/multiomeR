# Preprocessing and nucleus QC

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

The root `_targets.R` maps per-GEM-well targets over the active rows of `cfg_GEM_wells.tsv` and aggregation targets over the active aggregations; [Reading the graph views](graph_methodology.md) traces one aggregation through that mapping. This page covers per-GEM-well preprocessing and the successive nucleus filters up to the final WNN cell set. GEM-well settings are described in [GEM well table](../reference_GEM_wells.html). Configurable settings and their defaults are listed in the [parameter browser](../parameters.html). The AMULET and ATAC scDblFinder adaptations are compared with their references in [Algorithmic implementations](algorithm_validation.md).

{{< include _shared_methods/quality_control.md >}}

**Source:** [`extra_targets/per_GEM_well_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/per_GEM_well_targets.R), [`extra_targets/general_aggregation_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/general_aggregation_targets.R), [`R/parallel_GEM_well_preprocessing_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/parallel_GEM_well_preprocessing_helpers.R), [`R/amulet_BPCells_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/amulet_BPCells_helpers.R), [`R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R), [`R/processing_ATAC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_ATAC_helpers.R).

## Parallel preprocessing graph

This view covers per-GEM-well processing and QC, including optional ambient RNA correction, donor demultiplexing, doublet detection, barcode filtering, and handoffs into aggregation-level GEX and ATAC objects.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/parallel_v2.mmd")
```

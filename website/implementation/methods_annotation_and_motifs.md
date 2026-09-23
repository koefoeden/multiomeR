# Cell-type annotation and motif accessibility

This page covers marker-signature cluster annotation and motif-family accessibility; their targets appear in the GEX, ATAC and WNN graphs of [GEX, ATAC, batch correction and WNN](methods_GEX_ATAC_and_WNN.md). Configurable settings and their defaults are listed in the [parameter browser](../parameters.html). The BPCells-native UCell scorer is compared with UCell in [Algorithmic implementations](algorithm_validation.md#bpcells-native-ucell-scoring).

## Cell-type annotation

{{< include _shared_methods/cell_type_annotation.md >}}

**Source:** [`R/cluster_annotation_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/cluster_annotation_helpers.R), [`R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R), [`extra_targets/GEX_graph_and_cluster_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_graph_and_cluster_targets.R).

## Motif families and motif accessibility

{{< include _shared_methods/motif_accessibility.md >}}

**Source:** [`R/celltype_labeling_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/celltype_labeling_helpers.R), [`extra_targets/ATAC_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/ATAC_targets.R), [`resources/`](https://github.com/koefoeden/multiomeR/tree/main/resources).

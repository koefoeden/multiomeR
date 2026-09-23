# Peak–gene correlation

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`module_peak_gene_correlation/targets.R` maps only opted-in aggregations and binds the primary-module inputs they consume. Parameters use the `peak_gene_correlation` manifest scope, targets end in `.peak_gene_correlation.<aggregation>`, and plot checkpoint tags use `peak_gene_correlation`, outside the numbered QC selections. The prerequisites and module selector are described in [Peak–gene correlation](../downstream_peak_gene_correlation.html). Configurable settings and their defaults are listed in the [parameter browser](../parameters.html). The compiled kernels are compared with their references in [Algorithmic implementations](algorithm_validation.md).

{{< include _shared_methods/peak_gene_correlation.md >}}

**Source:** [`module_peak_gene_correlation/`](https://github.com/koefoeden/multiomeR/tree/main/module_peak_gene_correlation), [`R/peak_gene_correlation_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_correlation_helpers.R), [`R/peak_gene_filter_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_filter_helpers.R), [`R/peak_gene_hierarchical_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_hierarchical_helpers.R), [`R/peak_gene_finemapping_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_finemapping_helpers.R), [`src/peak_gene_REML.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/peak_gene_REML.cpp), [`src/peak_gene_KR.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/peak_gene_KR.cpp).

## Target graph

This view covers the TSS table, candidate peak–gene pairs, the WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/peak_gene_correlation_v2.mmd")
```

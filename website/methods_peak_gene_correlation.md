# Peak–gene correlation methods

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

This page describes the models of the optional peak–gene correlation module. [Peak–gene correlation](downstream_peak_gene_correlation.md) describes its prerequisites, configuration and outputs, and its settings and defaults are listed in the [parameter browser](parameters.html#workflow=peak_gene_correlation).

## Pseudobulks and support filters {#pseudobulks-and-support-filters}

{{< include _shared_methods/peak_gene_pseudobulks.md >}}

## Models and candidate links {#models-and-candidate-links}

{{< include _shared_methods/peak_gene_models.md >}}

The compiled kernels are compared with `lme4`, `pbkrtest` and `stats::p.adjust()` in [Algorithm validation](algorithm_validation.md).

## Source and target graph

`module_peak_gene_correlation/targets.R` maps only opted-in aggregations and binds the primary-module inputs they consume. Parameters use the `peak_gene_correlation` manifest scope, targets end in `.peak_gene_correlation.<aggregation>`, and plot checkpoint tags use `peak_gene_correlation`, outside the numbered QC selections.

**Source:** [`module_peak_gene_correlation/`](https://github.com/koefoeden/multiomeR/tree/main/module_peak_gene_correlation), [`R/peak_gene_correlation_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_correlation_helpers.R), [`R/peak_gene_filter_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_filter_helpers.R), [`R/peak_gene_hierarchical_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_hierarchical_helpers.R), [`R/peak_gene_finemapping_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_finemapping_helpers.R), [`src/peak_gene_REML.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/peak_gene_REML.cpp), [`src/peak_gene_KR.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/peak_gene_KR.cpp).

This view covers the TSS table, candidate peak–gene pairs, the WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/peak_gene_correlation_v2.mmd")
```

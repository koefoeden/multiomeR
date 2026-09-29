# Differential-analysis methods

```{r setup, include = FALSE}
pipeline_name <- "differential_analyses"
source("helpers/_setup.R")
```

This page describes the models of the optional differential-analysis module. [Differential analyses](downstream_differential_analyses.md) describes its prerequisites, configuration and outputs, and its settings and defaults are listed in the [parameter browser](parameters.html#workflow=differential_analyses).

## Cell-type composition

{{< include _shared_methods/cell_type_composition.md >}}

## Molecular pseudobulk analyses

{{< include _shared_methods/pseudobulk_differential_analyses.md >}}

## Source and target graph

`module_differential_analyses/targets.R` filters the aggregations that enabled the module, joins their module configuration, attaches symbols for the accepted WNN metadata and pseudobulk inputs, and maps the target files in its directory. One generic pseudobulk model family is instantiated for each tested feature matrix.

**Source:** [`module_differential_analyses/`](https://github.com/koefoeden/multiomeR/tree/main/module_differential_analyses), [`R/differential_analysis_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/differential_analysis_helpers.R), [`R/pseudobulk_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/pseudobulk_helpers.R), [`R/TF_activity_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/TF_activity_helpers.R).

This view covers the donor metadata, the pseudobulk and motif-family inputs, the CollecTRI network, the composition and feature models, and the cross-modality comparison.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/differential_analyses_v2.mmd")
```

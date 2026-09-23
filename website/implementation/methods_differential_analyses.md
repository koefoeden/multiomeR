# Differential analyses

```{r setup, include = FALSE}
pipeline_name <- "differential_analyses"
source("helpers/_setup.R")
```

`module_differential_analyses/targets.R` filters the aggregations that enabled the module, joins their module configuration, attaches symbols for the accepted WNN metadata and pseudobulk inputs, and maps the target files in its directory. One generic pseudobulk model family is instantiated for each tested feature matrix. The prerequisites and module selector are described in [Differential analyses](../downstream_differential_analyses.html). Configurable settings and their defaults are listed in the [parameter browser](../parameters.html).

## Cell-type composition

{{< include _shared_methods/cell_type_composition.md >}}

**Source:** [`module_differential_analyses/setup_and_cell_type_composition_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/module_differential_analyses/setup_and_cell_type_composition_targets.R), [`R/differential_analysis_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/differential_analysis_helpers.R).

## Molecular pseudobulk analyses

{{< include _shared_methods/pseudobulk_differential_analyses.md >}}

**Source:** [`module_differential_analyses/`](https://github.com/koefoeden/multiomeR/tree/main/module_differential_analyses), [`R/pseudobulk_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/pseudobulk_helpers.R), [`R/TF_activity_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/TF_activity_helpers.R), [`R/differential_analysis_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/differential_analysis_helpers.R).

## Target graph

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/differential_analyses_v2.mmd")
```

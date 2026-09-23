# Differential analyses

```{r setup, include = FALSE}
pipeline_name <- "differential_analyses"
source("helpers/_setup.R")
```

`module_differential_analyses/targets.R` filters the aggregations that enabled the module, joins their module configuration, attaches symbols for the accepted WNN metadata and pseudobulk inputs, and maps the target files in `module_differential_analyses/`. One generic pseudobulk model family is instantiated for each tested feature matrix.

The graph below is an orientation view; the target files hold the complete model and plotting commands. The methods, model routes and settings are in [Differential analyses methods](methods_differential_analyses.md), and the prerequisites and module selector in [Differential analyses](../downstream_differential_analyses.html).

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/differential_analyses_v2.mmd")
```

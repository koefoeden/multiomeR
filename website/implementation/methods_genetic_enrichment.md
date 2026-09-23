# Genetic enrichment

```{r setup, include = FALSE}
pipeline_name <- "genetic_enrichment"
source("helpers/_setup.R")
```

`module_genetic_enrichment/targets.R` selects the aggregations whose `modules` include `genetic_enrichment`, resolves one configured Open Targets study set per aggregation, attaches symbols for the primary-module inputs it consumes, and maps the target files in its directory. It does not check the species; the GRCh38 GWAS inputs assume a human aggregation. The release, method-selection, and interpretation contracts are described in [Genetic enrichment](../downstream_genetic_enrichment.html). Configurable settings and their defaults are listed in the [parameter browser](../parameters.html). The sparse SCAVENGE implementation is compared with its reference in [Algorithmic implementations](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

{{< include _shared_methods/genetic_enrichment.md >}}

**Source:** [`module_genetic_enrichment/`](https://github.com/koefoeden/multiomeR/tree/main/module_genetic_enrichment), [`R/open_targets_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/open_targets_helpers.R), [`R/GWAS_chromVAR_input_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_input_helpers.R), [`R/GWAS_chromVAR_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_helpers.R), [`R/GWAS_chromVAR_absolute_effect_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_absolute_effect_helpers.R), [`R/GWAS_chromVAR_contribution_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_contribution_helpers.R), [`R/SCAVENGE_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/SCAVENGE_helpers.R).

## Target graphs

### Single-nucleus and graph-based enrichment

This view covers the configured GWAS inputs, single-nucleus enrichment state, graph propagation, and downstream trait summaries.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd")
```

### Cell-type pseudobulk enrichment

This view covers the annotation-class pseudobulk deviations that are calculated separately from the nucleus-level results.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_cell_type_absolute_effect_v2.mmd")
```

### Cell-type contributions and locus attribution

This view covers the per-class contribution and locus-attribution branches that are pruned from the single-nucleus view.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_cell_type_contributions_v2.mmd")
```

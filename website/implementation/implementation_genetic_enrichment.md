# Genetic enrichment

```{r setup, include = FALSE}
pipeline_name <- "genetic_enrichment"
source("helpers/_setup.R")
```

`module_genetic_enrichment/targets.R` selects the aggregations whose `modules` include `genetic_enrichment`, resolves one configured Open Targets study set per aggregation, attaches symbols for the primary-module inputs it consumes, and maps the target files in `module_genetic_enrichment/`. It does not check the species; the GRCh38 GWAS inputs assume a human aggregation.

The methods and settings are in [Genetic enrichment methods](methods_genetic_enrichment.md), and the release, method-selection, and interpretation contracts in [Genetic enrichment](../downstream_genetic_enrichment.html).

## Single-nucleus and graph-based enrichment

This view covers the configured GWAS inputs, single-nucleus enrichment state, graph propagation, and downstream trait summaries. Additional cell-type contribution and locus-attribution branches may be pruned from this compact orientation view; use the manifest for the complete graph.

The sparse SCAVENGE reimplementation, deliberate graph and permutation differences, and validation evidence are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd")
```

## Cell-type pseudobulk enrichment

This view covers the annotation-class pseudobulk deviations that are calculated separately from the nucleus-level results.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_cell_type_absolute_effect_v2.mmd")
```

## Cell-type contributions and locus attribution

This view covers the per-cell-type contribution and locus-attribution branches that are pruned from the single-nucleus view above.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_cell_type_contributions_v2.mmd")
```
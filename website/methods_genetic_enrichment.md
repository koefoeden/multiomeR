# Genetic-enrichment methods

```{r setup, include = FALSE}
pipeline_name <- "genetic_enrichment"
source("helpers/_setup.R")
```

This page describes the analyses of the optional genetic-enrichment module. [Genetic enrichment](downstream_genetic_enrichment.md) describes its release, method-selection and interpretation contracts, and its settings and defaults are listed in the [parameter browser](parameters.html#workflow=genetic_enrichment). The module does not check the species; the GRCh38 GWAS inputs assume a human aggregation.

## Nucleus-level enrichment {#nucleus-level-enrichment}

{{< include _shared_methods/genetic_enrichment_nuclei.md >}}

## Cell-type enrichment and locus attribution {#cell-type-enrichment}

{{< include _shared_methods/genetic_enrichment_cell_types.md >}}

## SCAVENGE propagation {#scavenge-propagation}

{{< include _shared_methods/genetic_enrichment_SCAVENGE.md >}}

The sparse SCAVENGE implementation is compared with its reference in [Algorithm validation](algorithm_validation.md#sparse-scavenge-propagation).

## Source and target graphs

`extra_targets/module_genetic_enrichment/targets.R` selects the aggregations whose `modules` include `genetic_enrichment`, resolves one configured Open Targets study set per aggregation, attaches symbols for the primary-module inputs it consumes, and maps the target files in its directory.

**Source:** [`extra_targets/module_genetic_enrichment/`](https://github.com/koefoeden/multiomeR/tree/main/extra_targets/module_genetic_enrichment), [`R/open_targets_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/open_targets_helpers.R), [`R/GWAS_chromVAR_input_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_input_helpers.R), [`R/GWAS_chromVAR_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_helpers.R), [`R/GWAS_chromVAR_absolute_effect_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_absolute_effect_helpers.R), [`R/GWAS_chromVAR_contribution_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/GWAS_chromVAR_contribution_helpers.R), [`R/SCAVENGE_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/SCAVENGE_helpers.R).

### Nucleus-level enrichment and SCAVENGE

This view covers the GWAS inputs, the chromVAR background shared with the motif analysis, and SCAVENGE propagation over the WNN graph to trait-relevance summaries.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd")
```

### Cell-type absolute effects

This view covers the absolute-effect branch, which weights variants by posterior probability times effect size.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_cell_type_absolute_effect_v2.mmd")
```

### Cell-type contributions and locus attribution

This view covers the annotation-class pseudobulk deviations and their decomposition into contributing peaks, variants and loci.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/genetic_enrichment_cell_type_contributions_v2.mmd")
```

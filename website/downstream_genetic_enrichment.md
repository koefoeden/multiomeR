# Genetic enrichment

```{r setup, include = FALSE}
pipeline_name <- "genetic_enrichment"
source("helpers/_setup.R")
```

## When to use this module

Use this module to ask which cell types or nuclei have accessible regions overlapping genetic evidence for a human trait. It connects fine-mapped GWAS variants to ATAC peaks, calculates accessibility-based enrichment, and uses [`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE) to summarize trait relevance across related nuclei.

A **credible set** contains candidate causal variants at a GWAS locus, with probabilities from fine-mapping. Enrichment helps prioritize cellular contexts; it does not by itself identify a causal cell type, gene, or mechanism.

See the [example plots](gallery.md#genetic-enrichment) for representative results and the [implementation graph](implementation/implementation_genetic_enrichment.html) for upstream ATAC and WNN dependencies.

## Prerequisites

Before enabling the module, confirm that:

- the aggregation is human and you have reviewed its main results and cell-type annotations;
- WNN metadata and graph results, consensus peaks, ATAC counts, chromVAR objects, and GEX/ATAC embeddings are available;
- each configured `sourceId` represents the intended trait and population;
- the machine can download and retain the Open Targets study and credible-set Parquet datasets; and
- you have chosen whether to interpret individual nuclei, graph-smoothed scores, or cell-type summaries; these answer related but different questions.

## Outputs

| Result | What to inspect |
|---|---|
| Study and variant-to-peak tables | Which studies, variants, and accessible regions contributed |
| Single-nucleus deviations | Accessibility enrichment for each nucleus |
| SCAVENGE plots | Trait-relevance scores propagated through the cell-neighbor graph |
| Cell-type heatmaps and attribution tables | Enrichment by cell type and the loci or variants contributing to it |

## Configure

Add [`modules`](parameters.html#modules) to the existing human aggregation entry, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [genetic_enrichment]
```

Then create a matching row directly in `configuration/cfg_module_genetic_enrichment.yaml`.

```{.yaml filename="configuration/cfg_module_genetic_enrichment.yaml"}
your_aggregation:
  genetic_enrichment_GWAS_studies:
    lymphocyte_count:
      Category: positive_control
      sourceId: GCST90002388
      finemappingMethod: auto
```

`sourceId` values beginning with `GCST` use the pinned Open Targets datasets. Every other value is a local Parquet filename, resolved from the project root and tracked as a file target. Local files must satisfy the schema enforced by `validate_local_finemapped_GWAS_tibble()`; their study ID, fine-mapping method, build, credible-set probability, and provenance are read from the file rather than repeated in YAML.

The root workflow currently pins Open Targets release `26.03`. That release identifier is recorded in downstream metadata and determines the available studies, credible sets, and fine-mapping methods.

For `finemappingMethod: auto`, multiomeR selects the first available supported method in this order: `SuSie`, `SuSiE-inf`, then `PICS`. Specify a method explicitly when the method itself is part of the analysis contract; the workflow fails if that method is unavailable for the study.

## Run

Preview the selected module outputs before running them:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".your_aggregation")
)
```

## Review

Open the study-selection summaries, chromVAR summaries, SCAVENGE heatmaps and locus-contribution plots produced for your configured studies. The [output gallery](gallery.md#genetic-enrichment) shows one example per plot; interpretation belongs to each plot.

Runtime and disk use grow with studies, cells, graph representations, permutations, and attributed loci.

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=genetic_enrichment).

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_yaml_entry(module_config_file, "immune_human_2x")
```

</details>

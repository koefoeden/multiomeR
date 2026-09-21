# Aggregation configuration

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`cfg_aggregations.yaml` has one top-level entry per aggregation: a joint GEX, ATAC, and WNN analysis of one or more GEM wells. The committed file enables the two human GEM wells in `immune_human_2x`, with optional modules disabled. The mouse and ENCODE validation examples are inactive by default. Edit the file directly; the demo entries can stay as worked examples. This page describes the entry structure and lists every parameter. When to set each parameter, and how to review the effect, is given step by step in [Run your own analysis](main_running.md#steps).

## Minimal entry

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  aggregation_GEM_well_IDs: [your_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
  is_active: true
```

Every other parameter has a default from `cfg_pipeline_parameters.tsv`, listed in the [parameter reference](#parameter-reference) below. Add a parameter to the entry only when you want to change its default.

## Required keys

- [`aggregation_GEM_well_IDs`](parameters.html#aggregation_GEM_well_IDs): the `GEM_well_ID` values to combine. Each must be an active row of the [GEM well table](reference_GEM_wells.md), and all must use the same Cell Ranger reference.
- [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv): the [donor metadata table](reference_donor_metadata.md) for these GEM wells.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes): a named list of marker genes per expected cell type, used for cluster annotation and marker plots.
- [`is_active`](parameters.html#is_active): whether targets are constructed for the aggregation. Deactivate aggregations you are not ready to run before an unqualified `targets::tar_make()`.

## Marker genes and transcription factors

Replace the placeholder genes with symbols appropriate for the tissue and reference. A gene listed without a suffix or with a `+` suffix is a positive marker; a `-` suffix marks a gene that should be absent. The optional [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) list names transcription factors per cell type for the motif-activity plots. How the annotation uses these lists is described in [Output files and metadata](review_outputs.md#cluster-annotation).

## QC filters after peak calling

[`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object) lists dplyr filter expressions applied to the peak-based ATAC metrics of the combined object, for example:

``` {.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

Omit it or set it to `null` until step 4 of [Run your own analysis](main_running.md#steps) has shown the distributions.

## Optional modules

Omit [`modules`](parameters.html#modules) for the first run. After reviewing the main results, enable an optional analysis by listing its name and adding a matching entry for the aggregation in the module's own configuration file:

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

See [Differential analyses](downstream_differential_analyses.md) and [Genetic enrichment](downstream_genetic_enrichment.md) for the module entries and their parameters.

## Parameter reference {#parameter-reference}

The [standalone parameter browser](parameters.html) is generated from `cfg_pipeline_parameters.tsv`, using a shared snapshot of the public runtime defaults and validation schema. Choose the main workflow or an optional module, then search by name or purpose. Cards are grouped by whether a value is required, defaulted, or optional. Defaults are visible beside each parameter; open a row for its type and example. See the [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_aggregations.yaml) for a complete configuration.

<details>

<summary>Show the public <code>immune_human_2x</code> example</summary>

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_yaml_entry(aggregations_config_file, "immune_human_2x")
```

</details>
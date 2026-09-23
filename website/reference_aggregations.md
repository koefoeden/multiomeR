# Aggregation configuration

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`cfg_aggregations.yaml` in the [selected configuration directory](main_overview.md#configuration-directory) has one top-level entry per aggregation: a joint GEX, ATAC and WNN analysis of one or more GEM wells. This page describes the entry structure; the [running guide](main_running.md#steps) explains when to set each parameter and how to review its effect.

The public configuration contains these entries; inactive ones can stay as examples:

- `template_aggregation` (inactive): a starting point for your own entry.
- `immune_human_2x` (active): the public demo, combining two human GEM wells with optional modules disabled.
- `brain_mouse` and `ENCODE_heart_LV_6x` (inactive): mouse and ENCODE validation examples.
- `mixed_human_31x` (inactive): the aggregation behind the [output gallery](gallery.md), with all optional modules enabled.

## Minimal entry

``` {.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  aggregation_GEM_well_IDs: [my_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
  is_active: true
```

Every other parameter is optional or has a default, listed in the [parameter reference](#parameter-reference) below. Add a parameter to the entry only to change its default.

## Required keys

- [`aggregation_GEM_well_IDs`](parameters.html#aggregation_GEM_well_IDs): the `GEM_well_ID` values to combine. Each must be an active row of the [GEM well table](reference_GEM_wells.md), and all must use the same Cell Ranger ARC reference.
- [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv): the [donor metadata table](reference_donor_metadata.md) for these GEM wells.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes): marker genes per expected cell type, used for cluster annotation and marker plots.

[`is_active`](parameters.html#is_active) (default `true`) controls whether targets are constructed for the aggregation. Set it to `false` for aggregations you are not ready to run before an unqualified `targets::tar_make()`; this does not delete existing results.

## Marker genes and transcription factors

Replace the placeholder genes with symbols appropriate for the tissue. List at least two cell types, and use gene names from the Cell Ranger ARC reference. A gene without a suffix or with a `+` suffix is a positive marker; a `-` suffix marks a gene that should be absent. Each cell type needs at least one positive marker. The optional [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) names transcription factors per cell type for the motif-accessibility plots at checkpoint 7. [Cluster annotation](review_outputs.md#cluster-annotation) describes how the marker lists are used.

## QC filters after peak calling

[`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object) lists filter expressions over the peak-based ATAC metrics. Nuclei for which an expression is `TRUE` are removed, for example:

``` {.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

Omit it until the peak QC plots at [checkpoint 4](main_running.md#checkpoint-4) have shown the distributions; checkpoint 5 then shows which nuclei the filters remove.

## Optional modules

Omit [`modules`](parameters.html#modules) for the first run. After reviewing the primary module's results, list the optional modules to run:

``` {.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [differential_analyses]
```

Each listed module also needs an entry named after the aggregation in its own configuration file, in the same directory. The module pages describe these entries:

| Module | Configuration file | Page |
|---|---|---|
| `differential_analyses` | `cfg_module_differential_analyses.yaml` | [Differential analyses](downstream_differential_analyses.md) |
| `genetic_enrichment` | `cfg_module_genetic_enrichment.yaml` | [Genetic enrichment](downstream_genetic_enrichment.md) |
| `peak_gene_correlation` | `cfg_module_peak_gene_correlation.yaml` | [Peak–gene correlation](downstream_peak_gene_correlation.md) |

## Parameter reference {#parameter-reference}

The [parameter browser](parameters.html) lists every parameter of the primary module and the optional modules with its default, type and an example. Choose a workflow, then search by name or purpose. The [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_aggregations.yaml) shows complete entries.

<details>

<summary>Show the public <code>immune_human_2x</code> example</summary>

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_yaml_entry(aggregations_config_file, "immune_human_2x")
```

</details>

# Aggregation configuration

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`cfg_aggregations.yaml` in the [selected configuration directory](main_overview.md#configuration-directory) has one top-level entry per aggregation: a joint GEX, ATAC and WNN analysis of one or more GEM wells. You add an entry at [step 3 of Configure your data](main_running.md#steps); the running guide explains when to revise each setting at its checkpoint.

The public configuration contains these entries; inactive ones can stay as examples:

- `template_aggregation` (inactive): a starting point for your own entry.
- `immune_human_2x` (active): the public demo, combining two human PBMC GEM wells with optional modules disabled.
- `brain_mouse` (inactive): a mouse example.
- `ENCODE_heart_LV_6x` (inactive): six ENCODE left-ventricle GEM wells with a differential-analysis example.

`configuration_dev/` holds the larger public-data aggregations: `mixed_human_31x`, the aggregation behind the [output gallery](gallery.md) with all optional modules enabled, and `comparison_1x` to `comparison_20x`, nested subsets of its GEM wells for the manuscript's resource comparisons. The seven 10x Genomics and 24 ENCODE GEM wells of `mixed_human_31x` are processed locally with Cell Ranger ARC 2.1.0 and GRCh38-2024-A, so its PBMC well `healthy_PBMC_human_2024A` is separate from the demo's downloaded `healthy_PBMC_human`.

## Minimal entry

Copy `template_aggregation`, rename the copy and fill it in. The result needs at least:

``` {.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  is_active: true
  aggregation_GEM_well_IDs: [my_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
```

`is_active: true` replaces the template's `false`; an entry written from scratch is active by default. Every other parameter is optional or has a default, listed in the [parameter browser](parameters.html#workflow=aggregation). Add a parameter to the entry only to change its default.

## Keys

- [`aggregation_GEM_well_IDs`](parameters.html#aggregation_GEM_well_IDs): the `GEM_well_ID` values to combine, from the [GEM well table](reference_GEM_wells.md).
- [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv): the [donor metadata table](reference_donor_metadata.md) for these GEM wells.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes): marker genes per expected cell type, used for cluster annotation and marker plots.
- [`is_active`](parameters.html#is_active) (default `true`): whether targets are constructed for the aggregation. Set it to `false` for aggregations you are not ready to run before an unqualified `targets::tar_make()`; this does not delete existing results.

## Marker genes and transcription factors

Replace the placeholder genes with symbols appropriate for the tissue. List at least two cell types, and use gene names from the Cell Ranger ARC reference. A gene without a suffix or with a `+` suffix is a positive marker; a `-` suffix marks a gene that should be absent. Each cell type needs at least one positive marker. The optional [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) names transcription factors per cell type for the motif-accessibility plots at checkpoint 7. [Cluster annotation](review_outputs.md#cluster-annotation) describes how the marker lists are used.

## Rules checked {#rules}

When the pipeline starts, it checks every entry against the [parameter manifest](change_the_pipeline.md#parameter-manifest): unknown parameters fail, and the resolved values are checked for missingness, type, cardinality and allowed values. It also checks that every listed GEM well exists and is active, and that each listed module has an entry in its configuration file. When the targets run, all GEM wells of an aggregation must share the same Cell Ranger ARC reference, and every marker gene must be present in the GEX count matrix.

## QC filters after peak calling {#qc-filters-after-peak-calling}

Leave [`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object) unset for the first run. After reviewing the peak QC distributions at [checkpoint 4](main_running.md#checkpoint-4), list filters and run checkpoint 5, which applies them to the GEX-retained nuclei. Each filter is a complete R expression over the per-nucleus metadata:

``` {.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

A nucleus for which any expression is `TRUE` is excluded, and each expression is reported as a separate exclusion reason in the checkpoint 5 UpSet and retention plots. Filters can use the metrics that [`QC_metric_manifest.tsv`](https://github.com/koefoeden/multiomeR/blob/main/QC_metric_manifest.tsv) lists as available from checkpoint 4 or earlier. These cutoffs are examples, not recommendations for your tissue.

## Inherit settings from another aggregation {#inheritance}

An entry can start from the settings of other entries listed in `inherits`, for example to rerun an aggregation with more GEM wells:

``` {.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  aggregation_GEM_well_IDs: [GEM_well_1, GEM_well_2]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]

my_aggregation_all_wells:
  inherits: my_aggregation
  aggregation_GEM_well_IDs: [GEM_well_1, GEM_well_2, GEM_well_3]
```

The entry starts from the parameter defaults, applies each parent in the listed order and then its own values. Every parameter is inherited, including `is_active`, and a value replaces the inherited one as a whole: here the child's `aggregation_GEM_well_IDs` replaces the parent's list rather than extending it. Entries in the module configuration files can inherit in the same way.

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

## Public example {#public-example}

The public demo's entry is shown below; the [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_aggregations.yaml) shows every entry.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_yaml_entry(aggregations_config_file, "immune_human_2x")
```

# Genetic enrichment

```{r setup, include = FALSE}
pipeline_name <- "genetic_enrichment"
source("helpers/_setup.R")
```

## When to use

Use this module to ask which cell types or nuclei have accessible chromatin overlapping genetic evidence for a human trait. It weights ATAC peaks by the fine-mapping probabilities of the GWAS variants they contain, scores enrichment per cell type and per nucleus, and uses the [`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE) method to propagate nucleus scores through the WNN neighbor graph. [Method details](methods_genetic_enrichment.md)

A **credible set** contains candidate causal variants at a GWAS locus, with probabilities from fine-mapping. Enrichment helps prioritize cellular contexts; it does not by itself identify a causal cell type, gene or mechanism.

The public demo leaves this module off. Its configuration already contains a study list for `immune_human_2x`, shown under [Parameter reference](#parameter-reference); add `genetic_enrichment` to that aggregation's `modules` to try it.

## Prerequisites

Before enabling the module, confirm that:

- the aggregation is human and you have reviewed it through [checkpoint 8](main_running.md#checkpoint-8), including its cell-type annotations;
- each configured study represents the intended trait and population;
- the machine can download and keep about 4 GB of Open Targets Parquet data, stored once under `<store>/files/OpenTargets/`; and
- you have decided whether to interpret cell-type summaries, individual nuclei or graph-propagated scores; these answer related but different questions.

## Outputs

| Question | Output |
|---|---|
| Which variants and peaks carry each study's evidence? | Peak-weight summaries and study and variant-to-peak tables |
| Which cell types are enriched, and which matters most? | Cell-type chromVAR deviation heatmaps, which outline each study's top set of cell types |
| Which loci and variants drive a cell type's enrichment? | Locus-contribution heatmaps and bars, and variant detail plots |
| Which nuclei are trait-relevant? | SCAVENGE trait-relevance scores on the WNN graph, summarized by cell type and cluster |

## Configure

Add `genetic_enrichment` to the aggregation's [`modules`](parameters.html#modules) list, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [genetic_enrichment]
```

Then add an entry for `my_aggregation` to `cfg_module_genetic_enrichment.yaml` in the [selected configuration directory](main_overview.md#configuration-directory):

```{.yaml filename="cfg_module_genetic_enrichment.yaml"}
my_aggregation:
  genetic_enrichment_GWAS_studies:
    lymphocyte_count:
      Category: positive_control
      sourceId: GCST90002388
      finemappingMethod: auto
```

Each entry under [`genetic_enrichment_GWAS_studies`](parameters.html#genetic_enrichment_GWAS_studies) is one study. Its name labels the results, and `Category` groups studies in the plots.

A `sourceId` names an Open Targets GWAS study with credible sets, read from the Open Targets release `26.03` pinned by the workflow: a GWAS Catalog accession such as `GCST90002388`, or a FinnGen study such as `FINNGEN_R12_T2D`. The run fails if the release has no credible sets for the study. `finemappingMethod: auto`, the default, selects the first available method in the order `SuSie`, `SuSiE-inf`, `PICS`. Name a method instead when the method is part of the analysis; the run then fails if that method is unavailable for the study.

A `sourceId` ending in `.parquet` is instead the path of a local file, relative to the repository root. Omit `finemappingMethod` for local files: the method, like the study ID, genome build and provenance, is read from the file.

<details>
<summary>Columns of a local fine-mapped GWAS file</summary>

The file holds one GWAS on GRCh38, with one row per credible-set variant. Its columns follow the [Open Targets credible-set dataset](https://platform.opentargets.org/downloads/credible_set/metadata), with the `locus` array flattened into one row per variant, and that dataset documents their meaning. multiomeR uses these columns:

| Column | Content | Used for |
|---|---|---|
| `studyId`, `finemappingMethod`, `confidence`, `sampleSize`, `sourceRelease` | One value per file | Study labels and metadata |
| `credibleSetProbability` | One value, greater than 0 and at most 1 | Study metadata |
| `studyLocusId`, `credibleSetIndex` | Credible-set identifiers | Grouping variants into loci |
| `variantId` | `chromosome_position_REF_ALT` | Matching variants between studies |
| `chromosome`, `position` | GRCh38 coordinates, chromosome without `chr` | Overlap with peaks |
| `posteriorProbability` | 0 to 1 | Peak weights |
| `beta` | Effect size, or missing | Effect-weighted enrichment |
| `locusStart`, `locusEnd` | Bounds that contain each variant's `position` | Locus plots |

The file must also contain `schemaVersion` (`1`), `studyType` (`gwas`), `genomeBuild` (`GRCh38`), `sourceUrl`, `sourceSha256` (lowercase SHA-256) and `sourcePublication`, each with one value, and the columns `variantRepresentation`, `logBF`, `pValueMantissa`, `pValueExponent`, `standardError`, `r2Overall`, `is95CredibleSet` and `is99CredibleSet`. Each `studyId`, `studyLocusId` and `variantId` combination must be unique. A file with `credibleSetProbability` 0.95 or 0.99 may contain only variants in that set.

</details>

## Run

Preview the targets tagged for this module:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Runtime and disk use grow with the number of studies, nuclei and attributed loci.

## Review

Review the plots ([examples](gallery.md#genetic-enrichment)); each plot's subtitle and caption say what to look for.

```text
<store>/plots/my_aggregation/genetic_enrichment/
├── GWAS_peak_weights_barplot.png
├── cell_type_pseudobulk/
│   ├── chromVAR_deviation_heatmaps/
│   ├── chromVAR_locus_contribution_per_GWAS_heatmaps/
│   ├── chromVAR_locus_contribution_per_GWAS_faceted_bars_plots/
│   ├── chromVAR_absolute_effect_locus_contribution_per_GWAS_faceted_bars_plots/
│   └── chromVAR_variant_contribution_detail_plots/
└── single_nucleus/SCAVENGE/WNN_harmony_SNN/
    ├── TRS_heatmap/
    ├── TRS_UMAPs/
    └── TRS_cluster_summary_plot/
```

Read the tables behind the plots in R, for example the cell-type enrichment scores:

```{.r filename="R"}
targets::tar_read(
  chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.my_aggregation
)
```

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=genetic_enrichment). [Genetic-enrichment methods](methods_genetic_enrichment.md) describes the analyses and their key fixed values, links to the source files and shows the target structure.

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_yaml_entry(module_config_file, "immune_human_2x")
```

</details>

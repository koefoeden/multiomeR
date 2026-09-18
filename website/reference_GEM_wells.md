# GEM well table

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

<integrate the text below into the rest to produce a much more succinct guide>
a unique GEM_well_ID and an appropriate GEM_well_dataset label;

GEM_well_cellranger_arc_count_dir;

donors: GEM_well_n_donors and GEM_well_donor_id for a single-donor well, or GEM_well_donors_VCF_file to demultiplex several donors by genotype;

GEM_well_add_cellbender and GEM_well_cellbender_h5_file to replace the Cell Ranger counts with CellBender output;

further GEM_well_ columns for library-level variables you will need as batch covariates or plot variables later;

GEM_well_QC_exclude_list: NA for the first run, then filter expressions such as TSS.enrichment < 4 ;; nCount_RNA < 250;

GEM_well_is_active: TRUE.
`cfg_GEM_wells.tsv` has one row per GEM well: one `cellranger-arc count` output together with its donor assignment, optional inputs, and pre-aggregation QC filters. Aggregations refer to rows by `GEM_well_ID`. Edit the committed file directly; the demo rows can stay as worked examples. This page describes the columns. When to edit them, and how to review the effect, is step 1 of [Run your own analysis](main_running.md#steps).

## Minimal row

The example below describes one non-multiplexed GEM well from one donor. This vertical view is a reading aid; the saved TSV has one GEM well per row.

::: {.scrollable-table}

| Column | Example value |
|---|---|
| `GEM_well_ID` | `your_GEM_well` |
| `GEM_well_dataset` | `your_dataset` |
| `GEM_well_donor_id` | `donor_1` |
| `GEM_well_n_donors` | `1` |
| `GEM_well_cellranger_arc_count_dir` | `/path/to/your_GEM_well` |
| `GEM_well_add_cellbender` | `FALSE` |
| `GEM_well_cellbender_h5_file` | `NA` |
| `GEM_well_donors_VCF_file` | `NA` |
| `GEM_well_QC_exclude_list` | `NA` |
| `GEM_well_is_active` | `TRUE` |
| `GEM_well_multiplex_batch` | `batch_1` |

:::

## Identity and grouping

`GEM_well_ID` must be unique and is used in target names, output folders, and `aggregation_GEM_well_IDs`. `GEM_well_dataset` is a reporting label, not the donor or the biological design. Cross-well QC comparisons follow the wells selected for an aggregation.

## Input files

`GEM_well_cellranger_arc_count_dir` points to the directory containing `outs/`, not to `outs/` itself. The baseline pipeline requires:

``` {.text filename="Text"}
<GEM_well_cellranger_arc_count_dir>/outs/
|-- summary.csv
|-- filtered_feature_bc_matrix.h5
|-- atac_fragments.tsv.gz
|-- atac_fragments.tsv.gz.tbi
`-- per_barcode_metrics.csv
```

multiomeR identifies the reference from the FASTA/GTF hashes in `outs/atac_fragments.tsv.gz` and selects the matching `reference.json` under `reference_metadata/`. For another reference, add its JSON there. The fragment header must be intact and exactly one JSON must match. All GEM wells in an aggregation must use the same reference. Deployments that require an explicit `GEM_well_cellranger_arc_reference_json` column must supply that path for each well; follow the checkout README for its input contract.

`GEM_well_add_cellbender` set to `TRUE` replaces the Cell Ranger GEX counts with a CellBender H5 file produced outside multiomeR, given in `GEM_well_cellbender_h5_file` or found in the standard layout.

## Donors

`GEM_well_n_donors` states how many donors the GEM well contains. For a non-multiplexed well, set it to `1` and put the donor's identifier in `GEM_well_donor_id`; every called nucleus receives that ID, which must match a `donor_id` row in the [donor metadata table](reference_donor_metadata.md).

For a multiplexed well, supply the donors' genotypes in `GEM_well_donors_VCF_file`. `atac_possorted_bam.bam` is then also required in `outs/` for `cellsnp-lite`, and nuclei are assigned to donors by genotype.

## Pre-aggregation QC filters

`GEM_well_QC_exclude_list` contains zero or more complete R filter expressions separated by `;;`, for example:

```{.text filename="cfg_GEM_wells.tsv"}
TSS.enrichment < 4 ;; nucleosome_signal > 4 ;; nCount_RNA < 250
```

Expressions are evaluated individually against per-barcode metadata, preserving their order and their separate exclusion reasons. An empty field (`NA`) applies no pre-aggregation QC filters. Leave it empty for a new GEM well until step 1 of [Run your own analysis](main_running.md#steps) has shown its distributions; do not treat the demo's numerical cutoffs as recommendations for your tissue.

AMULET doublet detection runs as part of the standard QC calculation; see the [implementation and validation](implementation/algorithm_validation.html#bpcells-native-amulet) for details.

## Active flag

`GEM_well_is_active` controls whether processing targets are constructed for the GEM well. Every GEM well selected by an active aggregation must also be active. Before an unqualified `targets::tar_make()`, deactivate every GEM well you are not ready to run, including demo rows you do not need.

## Extra metadata columns

Add library-, run-, or batch-specific variables as further columns with a `GEM_well_` prefix, such as `GEM_well_multiplex_batch` or `GEM_well_cell_sorting`. They are joined to the cell metadata and can be used as Harmony covariates or plotting variables. Donor-level phenotypes belong in the [donor metadata table](reference_donor_metadata.md) instead. Apart from their key columns, the two tables must not reuse column names.

## Column dictionary

Use the [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_GEM_wells.tsv) to check the full set of columns.

The table below is a documentation snapshot of the two public demo wells, showing the core columns and one optional annotation. Bold columns must be present in the TSV; some allow an NA value. Scroll horizontally and focus or hover over a column's **i** button for its meaning. The other inactive rows and metadata columns in the public configuration remain available as examples.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_GEM_well_demo_table(
  GEM_well_config_file = "website/data/demo_GEM_wells.tsv",
  dictionary_file = "website/data/GEM_well_columns.tsv"
)
```

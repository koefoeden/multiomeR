# Donor metadata table

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

The donor metadata table is a TSV with one row per donor, keyed by `donor_id`. You prepare it at [step 2 of Configure your data](main_running.md#steps), and each aggregation points to its file with [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv).

## Minimal table

``` {.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

## Columns

- `donor_id` (required): one row for each donor ID that the nuclei receive; see [Matching donors to nuclei](#matching-donors-to-nuclei).
- Every other column is a donor-level variable, such as condition, age or sex. The variables become cell metadata: select them for plots with [`aggregation_categorical_vars`](parameters.html#aggregation_categorical_vars) or [`aggregation_continuous_vars`](parameters.html#aggregation_continuous_vars), and use them in the model formulas of the [differential analyses](downstream_differential_analyses.md).

Put library-, run- or batch-specific variables in the [GEM well table](reference_GEM_wells.md) instead, with a `GEM_well_` prefix.

## Matching donors to nuclei {#matching-donors-to-nuclei}

Every nucleus receives a `donor_id` from its GEM well: `GEM_well_donor_id` for a single-donor well, or the VCF sample name that Vireo assigns in a [pooled well](reference_GEM_wells.md#pooled-wells). Give each of those IDs one row, with identical spelling; the nuclei-per-donor plot at [checkpoint 1](main_running.md#checkpoint-1) shows the IDs the nuclei received. Nuclei whose donor is missing from the table get `NA` donor variables, and rows for donors without nuclei are ignored.

## Rules checked {#rules}

When the targets that read the table run, they check that:

- `donor_id` values are unique; and
- apart from `donor_id` and `GEM_well_ID`, no column name appears both here and in the [GEM well table](reference_GEM_wells.md).

## Extended table for differential analyses

The [differential analyses](downstream_differential_analyses.md) module can read additional donor-level model variables from a second table given in [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv). It must keep the same unique `donor_id` key. If it is not set, the module uses the aggregation's donor table.

## Public example {#public-example}

The public demo's table also lists the donors of other example aggregations; only `pbmc1` and `pbmc6` have nuclei in the demo:

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_file(
  "website/data/public_defaults/immune_human_dataset_donor_id_metadata.tsv",
  "configuration/immune_human_dataset_donor_id_metadata.tsv"
)
```

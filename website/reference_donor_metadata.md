# Donor metadata table

The donor metadata table is a TSV with one unique row per `donor_id`. Each aggregation points to one such file through [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv) in the [aggregation configuration](reference_aggregations.md). It is created in step 2 of [Run your own analysis](main_running.md#steps).

## Minimal table

``` {.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

## Matching donors to nuclei

Every nucleus receives a `donor_id` from its GEM well: the configured `GEM_well_donor_id` for a non-multiplexed well, or a genotype-based assignment for a well with a configured VCF. Each of those IDs must appear exactly once in this table; see the [GEM well table](reference_GEM_wells.md#fill-in-a-row).

## Which variables belong here

Put donor-specific phenotypes and covariates in this table, for example condition, age, or sex. Put library-, run-, or batch-specific variables in the GEM well table with a `GEM_well_` prefix. Apart from their key columns, the two tables must not reuse column names. Metadata file contents are validated when their targets run, not when the manifest is built.

## Extended table for differential analyses

The [differential analyses](downstream_differential_analyses.md) module can read additional donor-level model variables from a second table given in [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv). It must keep the same unique `donor_id` key. If it is not set, the module inherits the aggregation's donor table.

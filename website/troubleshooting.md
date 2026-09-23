# Troubleshooting

Find the target name and first error message in the run output. Fix that cause, then [rerun the same selection](#rerun-safely).

| What happened? | Start here |
|---|---|
| Setup fails before any data processing | [The manifest does not build](#the-manifest-does-not-build) |
| A target reports an error | [A run reports errored targets](#a-run-reports-errored-targets) |
| You need to check a result that a failed target used | [Inspect one upstream value](#inspect-one-upstream-value) |
| An error names an input file, GEM well, donor, or metadata column | [Input and metadata failures](#input-and-metadata-failures) |
| Workers do not start or are killed | [Controller and scheduler failures](#controller-and-scheduler-failures) |
| A completed result needs rebuilding | [A target is unexpectedly outdated](#a-target-is-unexpectedly-outdated) |

Run the commands below in the repository-root R session. For general techniques, see the [targets debugging guide](https://books.ropensci.org/targets/debugging.html).

## The manifest does not build

Reload the runtime and build the manifest:

```{.r filename="R"}
load_project_runtime()
targets::tar_manifest(callr_function = NULL)
```

Failures at this stage usually mean:

- `configuration.local` names a missing directory, or a `cfg_*` file is missing from the [selected configuration directory](main_overview.md#configuration-directory);
- a YAML parameter is unknown, misspelled, of the wrong type, or missing;
- an aggregation selects a GEM well that is absent from `cfg_GEM_wells.tsv`, or an active aggregation selects an inactive GEM well;
- an aggregation enables an unknown module, or the module's configuration has no entry for it;
- `crew_controllers.R` breaks the [controller rules](performance_distributed_computing.md#controller-rules); or
- an R package fails to load.

Correct the field or file named in the error, then repeat the manifest check.

## A run reports errored targets

List the errors, one target per distinct error message:

```{.r filename="R"}
list_distinct_errored_targets()
```

To add commands and stored tracebacks, filter errored targets by a regular expression that matches a target, aggregation, or module:

```{.r filename="R"}
list_distinct_errored_targets_w_tracebacks(
  target_name_pattern = "my_aggregation"
)
```

Some targets run separately for multiple groups; these runs are called **branches**, and their names end in a hash suffix. Copy the full target or branch name from the listing to inspect its saved workspace:

```{.r filename="R"}
inspect_target_workspace("categorical.UMAPs.8_multimodal_QC.my_aggregation")
```

The inspection reports the target command, stored error and traceback, and short summaries of its dependencies. Empty tables, zero-dimensional matrices, missing model columns, or unexpected labels usually point to an upstream data or configuration problem.

## Inspect one upstream value

Read a dependency by its full name and print only a small part of it:

```{.r filename="R"}
value <- targets::tar_read(metadata_w_cell_types_tibble.WNN.my_aggregation)
dim(value)
head(value)
```

## Input and metadata failures

Metadata files are checked when the targets that read them run. Compare the error with the [GEM well table](reference_GEM_wells.md), [donor metadata table](reference_donor_metadata.md), and [aggregation configuration](reference_aggregations.md) references, and check that:

- each `GEM_well_cellranger_arc_count_dir` contains the required `outs/` files, plus `atac_possorted_bam.bam` for genotype demultiplexing;
- `GEM_well_ID` values in `cfg_GEM_wells.tsv` and `donor_id` values in the donor metadata table are unique;
- every donor ID from the GEM well table or from genotype demultiplexing appears in the donor metadata table with identical spelling; and
- apart from these two keys, no column name appears in both tables.

## Controller and scheduler failures

If workers do not start, are killed, or no controller can run a target:

1. Check `crew_controllers.R` against the [controller rules](performance_distributed_computing.md#controller-rules), then reload it with `load_project_runtime()`.
2. An error starting with `No controller found` means that no tier offers the requested cores, RAM, or GPUs; add or enlarge a tier.
3. For scheduler controllers, read the scheduler's output and error logs, and check the queue, account, wall time, memory, CPU, module, and file-system settings.
4. If local workers are killed for lack of memory, lower their `workers` values as described in [Local execution](performance_distributed_computing.md#local-execution).

## A target is unexpectedly outdated

Ask `targets` which part of a narrow selection needs rebuilding:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation"),
  callr_function = NULL
)
```

Changes to code, configuration, input files, or any upstream target make downstream targets outdated. Look for the most upstream name in the result before assuming that the final target itself is the cause.

## Rerun safely

After fixing the cause, rerun the same selection; up-to-date results are reused:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation")
)
```

Run `targets::tar_make()` without `names` only when you want every active aggregation and enabled module built. Do not delete the store to debug: it holds the error records, and deleting it forces unrelated recomputation.

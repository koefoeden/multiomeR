# Introduction

Use this book to trace a result back to its code or change how multiomeR works. For installation, configuration, execution, and output inspection, start with the [user manual](../). You do not need to read this book to run the demo.

Use this book when you need to trace a configuration value into mapped targets, understand how the simplified graph views relate to the real `{targets}` graph, or decide where an implementation change belongs.

## Where to start

For a first implementation pass:

1. Read [Reading the graph views](graph_methodology.md) and follow its configuration-to-target trace.
2. Open the [primary module](implementation_main.md) graph for the modality or checkpoint you plan to change.
3. Use [Implementation conventions](implementation_conventions.md) to understand the relevant manifest, mapping, symbol, tag, and runtime contracts.
4. Read [Background and design philosophy](background_philosophy.md) when you need the rationale for the editable-workflow design.

The [differential analyses](implementation_differential_analyses.md), [genetic enrichment](implementation_genetic_enrichment.md) and [peak–gene correlation](implementation_peak_gene_correlation.md) chapters cover the optional module graphs.

The **Methods and parameters** chapters, from [Preprocessing and nucleus QC](methods_preprocessing_and_QC.md) to [Genetic enrichment](methods_genetic_enrichment.md), describe each stage and list every fixed value and every configurable setting, using the table layout defined in [Implementation conventions](implementation_conventions.md#methods-and-parameter-tables). Read them when you need the exact behaviour behind a result, or when deciding whether a change is a configuration edit or a code edit.

## Common entry points

| Change | Start with |
|---|---|
| Add or revise a YAML parameter | `cfg_pipeline_parameters.tsv`, then the owning config reader or target. |
| Find whether a threshold is configurable or fixed | The stage's chapter under **Methods and parameters**. |
| Change GEM well preprocessing | `_targets.R` mapping plus `extra_targets/per_GEM_well_targets.R`. |
| Change aggregation GEX, ATAC, or WNN processing | The corresponding graph section and `extra_targets/*_targets.R` file. |
| Inspect existing review selections | `[checkpoint:<name>]` description tags and the steps in [Run your own analysis](../main_running.html#steps). |
| Add a graph-visible target | Existing `[part_of_graph:<graph_id>]` tags and graph-pruning rules. |
| Change resource routing | `crew_controllers.R`, `packages/multiomeRCore/R/resource_helpers.R`, and the runtime bootstrap convention. |

If you are trying to run multiomeR rather than modify it, start with the [main manual](../).

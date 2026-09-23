# Introduction

Use this book to trace a result back to its code or change how multiomeR works. For installation, configuration, execution, and output inspection, start with the [user manual](../). You do not need to read this book to run the demo.

## Design

multiomeR keeps the analysis steps in an editable repository. Configuration covers common choices such as inputs, markers, dimensions, and models; R helpers and target definitions are available when a study needs a change beyond those settings. This flexibility also means that users must review which methods and assumptions fit their data.

- **Reuse completed work.** [`targets`](https://books.ropensci.org/targets/) records dependencies between results, so a change rebuilds only the affected parts of an analysis, and independent tasks such as GEM wells run concurrently when workers are available.
- **Keep large matrices on disk.** [BPCells](https://bnprks.github.io/BPCells/) provides disk-backed matrices and streaming operations; some steps still need substantial RAM. [Choose where the analysis runs](../performance_distributed_computing.html#what-to-expect) gives multiomeR examples.
- **Keep the analysis inspectable.** Separate targets make intermediate tables, matrices, and files available for inspection, and Seurat/Signac exports allow exploration outside the pipeline.

## Where to start

1. Read [Reading the graph views](graph_methodology.md) and follow its configuration-to-target trace.
2. Use [Implementation conventions](implementation_conventions.md) for the manifest, mapping, symbol, tag, and runtime contracts.
3. Open the page for the stage you plan to change.

The **Primary module** and **Optional modules** pages each describe what a stage does, state the key fixed values that are not exposed as parameters, link to the source files, and show the target graph. Configurable settings and their defaults are listed in the [parameter browser](../parameters.html). [Algorithmic implementations, deviations and validation](algorithm_validation.md) records how the reimplemented reference algorithms differ from their references and how they are tested.

## Common entry points

| Change | Start with |
|---|---|
| Add or revise a YAML parameter | `cfg_pipeline_parameters.tsv`, then the owning config reader or target. |
| Find whether a threshold is configurable or fixed | The [parameter browser](../parameters.html), then the stage's page and its source links. |
| Change GEM well preprocessing | `_targets.R` mapping plus `extra_targets/per_GEM_well_targets.R`. |
| Change aggregation GEX, ATAC, or WNN processing | The corresponding graph section and `extra_targets/*_targets.R` file. |
| Inspect existing review selections | `[checkpoint:<name>]` description tags and the steps in [Run your own analysis](../main_running.html#steps). |
| Add a graph-visible target | Existing `[part_of_graph:<graph_id>]` tags and graph-pruning rules. |
| Change resource routing | `crew_controllers.R`, `packages/multiomeRCore/R/resource_helpers.R`, and the runtime bootstrap convention. |

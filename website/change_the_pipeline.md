# Change the pipeline

multiomeR keeps its analysis steps in an editable repository. Configuration covers common choices such as inputs, markers, dimensions and models; when a study needs a change beyond those settings, edit the R helpers and target definitions. This page explains how configuration becomes target definitions, so that a change preserves the contracts other code relies on. What each analysis does is described on the methods pages, starting with [Primary-module methods](methods_primary_module.md).

## Design

- **Reuse completed work.** [`targets`](https://books.ropensci.org/targets/) records dependencies between results, so a change rebuilds only the affected parts of an analysis, and independent tasks such as GEM wells run concurrently when workers are available.
- **Keep large matrices on disk.** [BPCells](https://bnprks.github.io/BPCells/) provides disk-backed matrices and streaming operations; some steps still need substantial RAM. [Choose where the analysis runs](performance_distributed_computing.md#what-to-expect) gives multiomeR examples.
- **Keep the analysis inspectable.** Separate targets make intermediate tables, matrices, and files available for inspection, and Seurat/Signac exports allow exploration outside the pipeline.

This flexibility also means that users must review which methods and assumptions fit their data.

## Where to start a change

| Change | Start with |
|---|---|
| Add or revise a YAML parameter | `manifests/cfg_pipeline_parameters.tsv`, then the owning config reader or target; see [Parameter manifest](#parameter-manifest). |
| Find whether a threshold is configurable or fixed | The [parameter browser](parameters.html), then the stage's methods section and its source links. |
| Change GEM well preprocessing | The `_targets.R` mapping and `extra_targets/per_GEM_well_targets.R`; see [Mapping tibbles](#mapping-tibbles). |
| Change aggregation GEX, ATAC, or WNN processing | The stage's section in [Primary-module methods](methods_primary_module.md) and its `extra_targets/*_targets.R` file. |
| Inspect existing review selections | `[checkpoint:<name>]` [description tags](#target-metadata-tags) and the checkpoints in [Run your own analysis](main_running.md). |
| Add a graph-visible target | `[part_of_graph:<graph_id>]` tags; see [Graph views](#graph-views). |
| Change resource routing | `crew_controllers.R`, `packages/multiomeRCore/R/resource_helpers.R`, and the [runtime bootstrap](#runtime-bootstrap). |

## Trace one configured aggregation

The quickest way to understand the implementation is to follow one value across the graph:

1. `immune_human_2x` is a key in `cfg_aggregations.yaml`.
2. `read_aggregation_config_tibble()` resolves manifest defaults and inheritance into one aggregation row.
3. `build_aggregation_tibble()` filters active rows and adds symbols for GEM-well upstream targets, including those used for aggregation-level QC summaries.
4. The root `_targets.R` passes that row through `tar_map(names = aggregation, delimiter = ".")`.
5. A base target such as `multimodal_Seurat_object.8_multimodal_QC` becomes `multimodal_Seurat_object.8_multimodal_QC.immune_human_2x`.
6. Description tags make selected targets discoverable as checkpoints or graph nodes, while structured file helpers derive output paths from the active target name.

Inspect the exact target command and description without running it:

```r
targets::tar_manifest(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]8_multimodal_QC[.]immune_human_2x$"
  ),
  fields = c(name, command, description),
  callr_function = NULL
)
```

The following sections describe each step. Preserve these contracts unless a change is meant to replace one of them.

## Parameter manifest

`manifests/cfg_pipeline_parameters.tsv` is the schema for YAML-backed pipeline configuration. Each row defines one parameter for one scope: `aggregation`, or the name of an optional module. Its columns record the type, cardinality, default, missing-value rule, allowed values and description of the parameter. The YAML files then only need to specify values that differ from the manifest defaults, plus values that are required because their resolved value may not be missing.

Configuration readers resolve file names with `configuration_path()` in the [selected configuration directory](main_overview.md#configuration-directory). Aggregation and enabled-module settings are resolved during graph construction; disabled modules do not read their configuration.

At read time, the pipeline loads the manifest for a scope and parses each `default_value` as YAML. This allows defaults to be literal scalars, `NULL`, YAML lists, or evaluated YAML expressions such as `!expr 1:30`. Manifest defaults seed every config row before inheritance and row-specific overrides are applied.

The resolution order is:

1.  Start with manifest defaults for the requested scope.
2.  Resolve each parent listed in `inherits`.
3.  Overlay parent values onto the defaults.
4.  Overlay the child row onto the inherited values.
5.  Validate the fully resolved row.

[Inherit settings from another aggregation](reference_aggregations.md#inheritance) shows an example.

Validation is manifest-driven and happens before target construction. Unknown YAML parameters fail early. Resolved values are then checked for missingness, cardinality, type, and allowed values.

``` text
scalar      one non-list value
vector      atomic vector
list        list
named_list  list with non-empty names
```

The `data_type` column checks the R type after YAML parsing. `path` and `regex` are currently character-like schema labels; the validator does not check file existence or compile regular expressions. `allowed_values` is a comma-separated allow-list checked after coercing resolved values to character.

The [parameter browser](parameters.html) is generated from a snapshot of the public manifest; [`website/data/README.md`](https://github.com/koefoeden/multiomeR/blob/main/website/data/README.md) describes how to refresh it.

Module-specific YAML uses the same mechanism. Aggregations opt into modules through the aggregation config, and each enabled aggregation must have a matching module config row. Some cross-scope fallbacks are still implemented by target code rather than by manifest inheritance; for example, a module parameter may intentionally allow `NULL` and then fall back to an aggregation-level path during module setup.

## Mapping tibbles

The root `_targets.R` builds the target graph from mapping tibbles. Each mapping tibble is a row-wise contract: one row becomes one set of mapped target instances, and columns in that row become local symbols inside the corresponding `tarchetypes::tar_map()` block.

The core mapping flow is:

1.  `GEM_well_tibble_all` reads only the pre-aggregation processing columns from every row in the canonical `cfg_GEM_wells.tsv`.
2.  `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
3.  `aggregation_tibble` keeps active aggregations, validates their GEM well references against the complete view, and adds upstream target-symbol columns.
4.  `GEM_well_tibble` keeps GEM wells whose `GEM_well_is_active` value is true.
5.  `_targets.R` expands active GEM wells and aggregations with `tar_map()`, then appends module target files. Cross-GEM-well QC summaries use the aggregation's selected wells.

Aggregation targets read GEM-well annotations through keyed projection targets that expose only the columns a consumer needs. These projections are cache boundaries: editing an unrelated column of `cfg_GEM_wells.tsv` does not invalidate expensive consumers.

``` r
tarchetypes::tar_map(
  values = GEM_well_tibble,
  names = GEM_well_ID,
  delimiter = ".",
  source("extra_targets/per_GEM_well_targets.R")$value
)
```

With `GEM_well_ID = "healthy_PBMC_human"`, a target named `cellranger_summary_file` becomes `cellranger_summary_file.healthy_PBMC_human`. The same dot-delimited suffix convention is used for aggregations, module targets, and nested module maps.

``` text
active cfg_GEM_wells.tsv row -> GEM_well_tibble row    -> per-GEM-well targets
cfg_aggregations.yaml key   -> aggregation_tibble row -> per-aggregation targets
```

Aggregation rows may opt into optional modules through `modules`. `_targets.R` validates module names against the known module list, and module target files then filter `aggregation_tibble` to the active aggregations that requested that module. Each opted-in aggregation must have a matching module config row.

The naming convention is therefore compositional:

``` text
<target>.<GEM_well_ID>
<target>.<aggregation_name>
<module_target>.<module_name>.<aggregation_name>
<nested_module_target>.<nested_suffix>.<module_name>.<aggregation_name>
```

Because these suffixes become target names and cache identity, config keys should be stable, human-readable, and free of unnecessary punctuation. In particular, avoid dots in GEM well, aggregation, and module IDs unless there is a compelling reason.

## Target-symbol columns

Mapped target tables sometimes need to carry references to other mapped targets. multiomeR represents those references as columns of `rlang` symbols. Each row stores the upstream target symbols that should be spliced into downstream target commands generated for that row.

`add_GEM_well_target_syms()` adds one such list-column per base target name, named `aggregation_<target>_syms`, from each aggregation's `aggregation_GEM_well_IDs`:

``` r
aggregation_tibble |>
  add_GEM_well_target_syms("GEX_counts_BPCells_matrix")
```

For an aggregation whose `aggregation_GEM_well_IDs` are `c("rx1", "rx2")`, this creates a row value equivalent to:

``` r
rlang::syms(c(
  "GEX_counts_BPCells_matrix.rx1",
  "GEX_counts_BPCells_matrix.rx2"
))
```

The aggregation target can then consume the row-local symbol list directly:

``` r
combined_counts_matrix <- purrr::reduce(
  aggregation_GEX_counts_BPCells_matrix_syms,
  cbind
)
```

Column names should describe the downstream scope, the upstream target, and the fact that the value is a symbol list. The `aggregation_*_syms` columns, such as `aggregation_GEX_counts_BPCells_matrix_syms`, splice per-GEM-well targets into aggregation-level targets.

Module target files also need aggregation-specific references to primary-module targets. For this, `add_aggregation_target_syms()` creates one symbol per row, suffixed by the aggregation name. These columns are named like the target they replace rather than with `*_syms`, because each cell is a single symbol rather than a list.

``` r
differential_analyses_tibble |>
  add_aggregation_target_syms(c(
    "metadata_w_cell_types_tibble.WNN",
    "pseudobulk_counts_matrix.GEX",
    "organism_chr"
  ))
```

For aggregation `PBMC`, the column `metadata_w_cell_types_tibble.WNN` contains the symbol `metadata_w_cell_types_tibble.WNN.PBMC`. Inside a module target, the command can be written against the unsuffixed local name; `tar_map()` resolves it to the aggregation-specific upstream target for that row.

## Target metadata tags

multiomeR stores lightweight target metadata in the `description` argument of `targets::tar_target()` and `tarchetypes::tar_file()` calls. The descriptions should remain readable prose, with bracketed tags appended when a target needs to be discoverable from the manifest.

``` r
targets::tar_target(
  name = harmony_embeddings_matrix.GEX,
  description = "Harmony-corrected SCTransform GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:WNN]",
  command = ...
)
```

The tag families are:

``` text
[checkpoint:<name>]             review or execution checkpoint
[part_of_graph:<graph_id>]      curated membership in an implementation graph
[resource_observation:<note>]   compact empirical resource note
```

`[checkpoint:<name>]` marks targets selectable with `targets::tar_described_as()`. The numbered primary-module groups are listed in `manifests/QC_checkpoint_manifest.tsv`; optional module groups remain unnumbered. Selection matches description substrings; include the closing `]` to match a complete checkpoint tag. Dependencies still come from the target commands. [Run your own analysis](main_running.md) explains each checkpoint.

Numbered checkpoint plot targets end in the checkpoint name, with hyphens replaced by underscores, before the mapped dataset or aggregation suffix; [Output folders](review_outputs.md#output-folders) explains how the name sets the plot folder. Only plot targets use this naming convention; computational and metadata targets retain their modality suffixes.

`[part_of_graph:<graph_id>]` marks targets that should stay visible in a named implementation graph after graph-pruning helpers remove less informative intermediate nodes. This is the strictest tag family: `graph_id` must contain only letters, numbers, and underscores, and helper code parses these tags directly from target descriptions. A target may belong to several graph views.

``` r
description = paste(
  "Build the lightweight BPCells-backed chromVAR RSE",
  "[part_of_graph:ATAC]",
  "[part_of_graph:seurat_export]",
  "[part_of_graph:genetic_enrichment_single_nucleus]"
)
```

`[resource_observation:<note>]` keeps a short, dated runtime or memory observation next to the target that produced it. It is documentation, not a resource-estimation system.

Use tags only when they create a durable handle for readers, graph helpers, or checkpoint commands. Ordinary internal dependencies can stay untagged.

## Graph views {#graph-views}

The target graphs on the methods pages are simplified views of the real `targets` dependency graph, meant to make the workflow easier to reason about before reading the target code. [`graphs_v2.R`](https://github.com/koefoeden/multiomeR/blob/main/website/figures/human_curated/graphs_v2.R) generates one view per `[part_of_graph:<graph_id>]` tag: it keeps the tagged targets, bypasses untagged intermediate targets while preserving the dependencies between tagged ones, replaces configured suffixes with placeholders such as `<aggregation_name>`, and merges duplicate labels.

The views are orientation aids, not alternate target definitions. A target can be absent because it lacks the tag for that view or was bypassed as a lower-level implementation detail. Use `targets::tar_manifest()`, `targets::tar_network()` or the source target files when exact completeness matters. The [figure README](https://github.com/koefoeden/multiomeR/blob/main/website/figures/human_curated/README.md) describes how to regenerate the views.

## Runtime bootstrap

multiomeR assumes that the repository runtime is bootstrapped before the target graph is inspected or run. The root `.Rprofile` is intentionally minimal: it sources `R/bootstrap_helpers.R` and calls `load_project_runtime()`.

`load_project_runtime()` is the single entry point for:

1.  loading core workflow packages and conflict preferences,
2.  sourcing generally reusable helpers from `packages/multiomeRCore/R`,
3.  sourcing pipeline-specific helpers from the root `R/` directory,
4.  applying global plotting and `{targets}` options,
5.  sourcing `crew_controllers.R` from the selected configuration directory and installing controller resources.

The nested `multiomeRCore` directory is both ordinary editable pipeline source and an installable package boundary for standalone repositories. multiomeR does not install or attach that package itself: `targets::tar_source()` loads the same implementation files before the root helpers. Keep domain-specific code under `R/`, but do not duplicate the generally reusable implementations there.

For commands that intentionally bypass startup side effects, source the bootstrap helper directly and then load the runtime:

``` r
source("R/bootstrap_helpers.R")
load_project_runtime()
targets::tar_manifest(callr_function = NULL)
```

`load_project_runtime()` keeps no state: each call reloads the packages, helpers, options, and controllers. Call it again when the current R session may be stale, such as after changing helper files, switching checkout roots, editing `crew_controllers.R`, or reusing a long-lived interactive session.

Project-root detection walks upward from the current working directory until it finds `pixi.toml`. Bootstrap commands should therefore be run from inside the multiomeR checkout.

Controller loading is part of the runtime contract, not a later execution detail. `crew_controllers.R` must return a named list with `controller_resources_tibble` and `controller_list`. The bootstrap validates that shape, installs a grouped `crew` controller into `{targets}`, and stores the resource table for `get_tar_resources()`.

``` r
targets::tar_target(
  example_target,
  example_function(),
  resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
)
```

If `get_tar_resources()` is called before controller resources are loaded, it fails deliberately with an instruction to call `load_project_runtime()` first. Scheduler-specific examples are in [Choose where the analysis runs](performance_distributed_computing.md); the implementation contract is that target code can request resources declaratively once the runtime has been loaded.

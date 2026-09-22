# Implementation conventions

This chapter explains how configuration becomes target definitions: the parameter manifest supplies defaults and validation rules, mapping tables define repeated analyses, and target symbols connect their dependencies. Description tags support result selection and graph views. The final section covers project startup and resource configuration.

## Target metadata tags

multiomeR stores lightweight target metadata in the `description` argument of `targets::tar_target()` and `tarchetypes::tar_file()` calls. The descriptions should remain readable prose, with bracketed tags appended when a target needs to be discoverable from the manifest.

``` r
targets::tar_target(
  name = harmony_embeddings_matrix.GEX,
  description = "Harmony-corrected SCTransform GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:WNN]",
  command = ...
)
```

The currently meaningful tag families are:

``` text
[checkpoint:<name>]             review or execution checkpoint
[part_of_graph:<graph_id>]      curated membership in an implementation graph
[resource_observation:<note>]   compact empirical resource note
```

`[checkpoint:<name>]` marks targets selectable with `targets::tar_described_as()`. The eight numbered main-pipeline groups are listed in `QC_checkpoint_manifest.tsv`; optional module groups remain unnumbered. UMAP parameter sweeps belong to their modality's numbered checkpoint, and compatibility objects belong to GEX checkpoint 3 or multimodal checkpoint 8. Selection matches description substrings; include the closing `]` to match a complete checkpoint tag. Dependencies still come from the target commands. [Run your own analysis](../main_running.html#steps) explains each boundary; acceptance criteria depend on the study.

Numbered checkpoint plot targets end in the checkpoint name with hyphens replaced by underscores, before the mapped dataset or aggregation suffix. For example, `VizDimLoadings_plots.2_GEX_PCA_QC.my_aggregation` writes beneath `<store>/plots/my_aggregation/2_GEX_PCA_QC/`. Only plot targets use this naming convention; computational and metadata targets retain their modality suffixes.

`[part_of_graph:<graph_id>]` marks targets that should stay visible in a named implementation graph after graph-pruning helpers remove less informative intermediate nodes. This is the strictest tag family: `graph_id` must contain only letters, numbers, and underscores, and helper code parses these tags directly from target descriptions. A target may belong to several graph views.

``` r
description = paste(
  "Build the lightweight BPCells-backed chromVAR RSE",
  "[part_of_graph:ATAC]",
  "[part_of_graph:seurat_export]",
  "[part_of_graph:genetic_enrichment_single_nucleus]"
)
```

`[resource_observation:<note>]` is currently best treated as provisional documentation. It is useful when a target has a compact empirical runtime or memory observation worth keeping near the target definition, but it is not yet a structured resource-estimation system. Keep these notes short, dated when relevant, and self-explanatory.

Use tags only when they create a durable handle for readers, graph helpers, or checkpoint commands. Ordinary internal dependencies can stay untagged.

## Parameter manifest

`cfg_pipeline_parameters.tsv` is the schema for YAML-backed pipeline configuration. Each row defines one parameter for one scope:

``` text
aggregation
differential_analyses
genetic_enrichment
peak_gene_correlation
```

For each parameter, the manifest records its name, type, cardinality, default value, missing-value rule, allowed values, example values, topic, graph/module ownership, and human description. The YAML files then only need to specify values that differ from the manifest defaults, plus values that are required because their resolved value may not be missing.

At read time, the pipeline loads the manifest for a scope and parses each `default_value` as YAML. This allows defaults to be literal scalars, `NULL`, YAML lists, or evaluated YAML expressions such as `!expr 1:30`. Manifest defaults seed every config row before inheritance and row-specific overrides are applied.

The resolution order is:

1.  Start with manifest defaults for the requested scope.
2.  Resolve each parent listed in `inherits`.
3.  Overlay parent values onto the defaults.
4.  Overlay the child row onto the inherited values.
5.  Validate the fully resolved row.

``` yaml
immune_human_2x:
  aggregation_GEM_well_IDs: [healthy_PBMC_human, lymphoma_lymph_human]
  aggregation_GEX_marker_genes:
    B: [MS4A1, CD79A]
    T: [TRAC, CD3D]

PBMC_human_6x:
  inherits: immune_human_2x
  aggregation_GEM_well_IDs:
    - healthy_PBMC_human
    - pbmc_10k_chromium_controller
  modules: [genetic_enrichment]
```

Validation is manifest-driven and happens before target construction. Unknown YAML parameters fail early. Resolved values are then checked for missingness, cardinality, type, and allowed values.

``` text
scalar      one non-list value
vector      atomic vector
list        list
named_list  list with non-empty names
```

The `data_type` column checks the R type after YAML parsing. `path` and `regex` are currently character-like schema labels; the validator does not check file existence or compile regular expressions. `allowed_values` is a comma-separated allow-list checked after coercing resolved values to character.

The website renders parameter tables from a generated snapshot of the public runtime manifest. Refresh that snapshot with the shared documentation inputs when public defaults change; private configuration is never used for these tables.

Module-specific YAML uses the same mechanism. Aggregations opt into modules through the aggregation config, and each enabled aggregation must have a matching module config row. Some cross-scope fallbacks are still implemented by target code rather than by manifest inheritance; for example, a module parameter may intentionally allow `NULL` and then fall back to an aggregation-level path during module setup.

## Mapping tibbles

The root `_targets.R` builds the target graph from mapping tibbles. Each mapping tibble is a row-wise contract: one row becomes one set of mapped target instances, and columns in that row become local symbols inside the corresponding `tarchetypes::tar_map()` block.

The core mapping flow is:

Configuration readers use `configuration_path()` to resolve a basename within `configuration/` or the directory named by the ignored `configuration.local`. Relative selections are anchored at the repository root. This selection does not change data-path interpretation or `_targets.yaml`. The selected GEM-well path is a graph global consumed by the file target, so changing directories also changes its dependency. Aggregation and enabled-module settings are resolved during graph construction. Disabled modules do not read their configuration.

1.  `GEM_well_tibble_all` reads only the pre-aggregation processing columns from every row in the canonical `cfg_GEM_wells.tsv`.
2.  `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
3.  `aggregation_tibble` keeps active aggregations, validates their GEM well references against the complete view, and adds upstream target-symbol columns.
4.  `GEM_well_tibble` keeps GEM wells whose `GEM_well_is_active` value is true.
5.  `_targets.R` expands active GEM wells and aggregations with `tar_map()`, then appends module target files. Cross-GEM-well QC summaries use the aggregation's selected wells.

Within each aggregation, `GEM_well_metadata_tibble` reads the same canonical file, subsets it to `aggregation_GEM_well_IDs`, and preserves that order. Cheap keyed projection targets then expose only the columns requested for SCT, Harmony or configured analyses. Complete non-processing annotations are joined only for explicit export objects. These projection targets are cache boundaries: a newly added or edited online column can update the canonical table without changing expensive consumers whose selected view is identical.

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
active cfg_GEM_wells.tsv row -> GEM_well_tibble row    -> per GEM well targets
cfg_aggregations.yaml key   -> aggregation_tibble row -> per-aggregation targets
```

Aggregation rows may opt into optional modules through `modules`. `_targets.R` validates module names against the known module list, and module target files then filter `aggregation_tibble` to the active aggregations that requested that module. Each opted-in aggregation must have a matching module config row.

The naming convention is therefore compositional:

``` text
<target>.<GEM_well_ID>
<target>.<dataset_name>
<target>.<aggregation_name>
<module_target>.<module_name>.<aggregation_name>
<nested_module_target>.<nested_suffix>.<module_name>.<aggregation_name>
```

Because these suffixes become target names and cache identity, config keys should be stable, human-readable, and free of unnecessary punctuation. In particular, avoid dots in GEM well, dataset, aggregation, and module IDs unless there is a compelling reason.

## Target-symbol columns

Mapped target tables sometimes need to carry references to other mapped targets. multiomeR represents those references as columns of `rlang` symbols. Each row stores the upstream target symbols that should be spliced into downstream target commands generated for that row.

The compact constructor is `target_sym_col()`. It records a base target name, the source column containing suffixes, the separator, and an optional transform. `add_target_sym_cols()` then turns those specifications into list-columns of `rlang::syms()`.

``` r
aggregation_tibble |>
  add_target_sym_cols(
    aggregation_GEX_counts_BPCells_matrix_syms =
      target_sym_col("GEX_counts_BPCells_matrix", "aggregation_GEM_well_IDs")
  )
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

Column names should describe the downstream scope, the upstream target, and the fact that the value is a symbol list. The `aggregation_*_syms` columns, such as `aggregation_GEX_counts_BPCells_matrix_syms`, splice per GEM well targets into aggregation-level targets.

Module target files also need aggregation-specific references to main-pipeline targets. For this, `add_aggregation_target_syms()` creates one symbol per row, suffixed by the aggregation name. These columns are named like the target they replace rather than with `*_syms`, because each cell is a single symbol rather than a list.

``` r
differential_analyses_tibble |>
  add_aggregation_target_syms(c(
    "metadata_w_cell_types_tibble.WNN",
    "pseudobulk_counts_matrix.GEX",
    "organism_chr"
  ))
```

For aggregation `PBMC`, the column `metadata_w_cell_types_tibble.WNN` contains the symbol `metadata_w_cell_types_tibble.WNN.PBMC`. Inside a module target, the command can be written against the unsuffixed local name; `tar_map()` resolves it to the aggregation-specific upstream target for that row.

## Runtime bootstrap

multiomeR assumes that the repository runtime is bootstrapped before the target graph is inspected or run. The root `.Rprofile` is intentionally minimal: it sources `R/bootstrap_helpers.R` and calls `load_project_runtime()`.

`load_project_runtime()` is the single entry point for:

1.  loading core workflow packages and conflict preferences,
2.  sourcing generally reusable helpers from `packages/multiomeRCore/R`,
3.  sourcing pipeline-specific helpers from the root `R/` directory,
4.  applying global plotting and `{targets}` options,
5.  sourcing `crew_controllers.R` and installing controller resources.

The nested `multiomeRCore` directory is both ordinary editable pipeline source and an installable package boundary for standalone repositories. multiomeR does not install or attach that package itself: `targets::tar_source()` loads the same implementation files before the root helpers. Keep domain-specific code under `R/`, but do not duplicate the generally reusable implementations there.

For commands that intentionally bypass startup side effects, source the bootstrap helper directly and then load the runtime:

``` r
source("R/bootstrap_helpers.R")
load_project_runtime(force = TRUE)
targets::tar_manifest(callr_function = NULL)
```

Bootstrap state is cached in `bootstrap_state_env`. This avoids reloading packages, re-sourcing helpers, reapplying target options, reassigning patches, and reloading controllers on every call. Use `force = TRUE` when the current R session may be stale, such as after changing helper files, switching checkout roots, editing `crew_controllers.R`, or reusing a long-lived interactive session.

Project-root detection walks upward from the current working directory until it finds `pixi.toml`. Bootstrap commands should therefore be run from inside the multiomeR checkout.

Controller loading is part of the runtime contract, not a later execution detail. `crew_controllers.R` must return a named list with `controller_resources_tibble` and `controller_list`. The bootstrap validates that shape, installs a grouped `crew` controller into `{targets}`, and stores the resource table for `get_tar_resources()`.

``` r
targets::tar_target(
  example_target,
  example_function(),
  resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
)
```

If `get_tar_resources()` is called before controller resources are loaded, it fails deliberately with an instruction to call `load_project_runtime()` first. Scheduler-specific examples belong in the main manual's [Choose where the analysis runs](../performance_distributed_computing.html) page; the implementation contract is that target code can request resources declaratively once the runtime has been loaded.

## Methods and parameter tables {#methods-and-parameter-tables}

The chapters in the **Methods and parameters** part describe each analysis stage and list every setting that determines its result. They are the single place where exact values are recorded; the manuscript supplement describes the same algorithms without values. Each chapter uses one table layout:

| Column | Content |
|---|---|
| Step | The analysis step, in the order the targets run. |
| Setting | The quantity or method choice. |
| Status | `Configurable` or `Hardcoded`, as defined below. |
| Value | Shown only for hardcoded settings. |
| Source | Where the setting lives. |

A setting is **configurable** when a row in `cfg_pipeline_parameters.tsv` controls it, directly or as a field inside a nested manifest parameter such as a model specification, or when a column of `cfg_GEM_wells.tsv` controls it. Configurable rows link to the [parameter browser](../parameters.html), which renders the current default from the public manifest snapshot. They never repeat the default in prose, so a default change needs no edit outside the manifest.

A setting is **hardcoded** when no manifest row or GEM-well column controls it. Hardcoded rows show the value and the source, which is one of:

``` text
helper default   default argument of an R helper that the calling target does not override
target literal   a literal passed explicitly in a target command
inline literal   a constant inside a function body
```

Helper defaults are the easiest to expose as parameters later; inline literals require a code change. Changing any hardcoded value invalidates the affected targets on the next run, like any other code change.

The demonstration settings in the manuscript supplement are resolved values for one aggregation and are not repeated here.

## How to read the rest of the implementation book

These conventions are the connective tissue behind the graph chapters. The parameter manifest explains why config rows can be compact. Mapping tibbles explain why target names have stable suffixes. Target-symbol columns explain how mapped targets pass sets of upstream targets across graph levels. Target metadata tags explain why some nodes remain visible in curated graph views. The bootstrap contract explains why helper functions, controller resources, and target options are available before `_targets.R` is evaluated.

When modifying the implementation, preserve these contracts unless the change is explicitly meant to replace one of them.

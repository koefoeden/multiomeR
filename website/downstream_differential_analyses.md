# Differential analyses

## When to use

Use this module to test how cell-type proportions, gene expression or chromatin accessibility differ with a donor condition or phenotype. Donors are the biological replicates: proportions are modelled per donor, and molecular measurements are summed into **pseudobulks**, one per cell type and donor.

The module does not create replication. The donors, covariates, formula and contrasts must support the intended comparison. For this reason the public demo leaves the module off: one healthy PBMC donor and one lymphoma lymph-node donor cannot separate condition, donor and tissue effects.

## Prerequisites

Before enabling the module, confirm that:

- you have reviewed the aggregation through [checkpoint 8](main_running.md#checkpoint-8), including its cell-type annotations;
- the donor metadata has one row per `donor_id` and every variable used in a model;
- model variables describe donors, not individual nuclei; and
- each compared group has enough donors for the design and contrasts.

If the models need variables beyond the aggregation's donor metadata, supply a second table in [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv) with the same unique `donor_id` key.

## Outputs

| Question | Output family |
|------------------------------------|------------------------------------|
| Do cell-type proportions differ? | `cell_type_composition` |
| Which genes change expression? | `gene_expression` |
| Which peaks change accessibility? | `chromatin_accessibility` |
| Which motif families change accessibility? | `motif_family_accessibility` (JASPAR) |
| Which regulators show altered expression-based activity? | `transcription_factor_activity` (CollecTRI) |

Motif-family accessibility summarizes ATAC data; transcription-factor activity is inferred from GEX data. The module also produces pseudobulk-depth and model diagnostics, comparisons across the feature families, and Hallmark and Reactome gene-set tests for gene expression.

## Configure

Add `differential_analyses` to the aggregation's [`modules`](parameters.html#modules) list, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [differential_analyses]
```

Then add an entry for `my_aggregation` to `cfg_module_differential_analyses.yaml` in the [selected configuration directory](main_overview.md#configuration-directory):

```{.yaml filename="cfg_module_differential_analyses.yaml"}
my_aggregation:
  differential_analyses_cell_type_composition_models:
    condition_abundance:
      formula: ~ condition
      contrast_specs_vec:
        treated_vs_control: conditiontreated
      plot_phenotype_vars: condition
      color_by: condition
  differential_analyses_pseudobulk_models:
    condition_model:
      formula: ~ 0 + cluster + condition
      contrast_specs_vec:
        treated_vs_control: conditiontreated
```

Each branch holds named models; the model and contrast names label the outputs.

- **Abundance models** ([`differential_analyses_cell_type_composition_models`](parameters.html#differential_analyses_cell_type_composition_models)) fit each cell type's share of a donor's nuclei. Use a one-sided, fixed-effects formula; the response is added for you. `plot_phenotype_vars` and `color_by` choose the plotted variables. Omit the field to skip this branch.
- **Feature models** ([`differential_analyses_pseudobulk_models`](parameters.html#differential_analyses_pseudobulk_models)) test the four feature families. The formula can use `cluster`, the cell type of each pseudobulk, and donor metadata variables. The example estimates one condition effect across cell types; to test within one cell type, add `cell_type_subset` and use `~ condition`.
- **Contrasts** in `contrast_specs_vec` are coefficient names of the model matrix, or linear combinations of them. A text variable's alphabetically first value is the reference, so `conditiontreated` is treated minus control. Replace the example variable, formula and contrast with your own.

Optional fields narrow a model:

- `donor_ids` restricts the donors in either branch;
- `cell_type_subset` restricts the cell types whose pseudobulks a feature model tests;
- `GEM_well_IDs` restricts the GEM wells that define an abundance model's population; and
- `cell_types_to_test` restricts the cell types an abundance model tests, while every retained nucleus still counts towards its donor's total.

Donors missing a model variable are excluded. The branches use different models, so their effect sizes are not directly comparable. The methods describe [abundance models](implementation/methods_differential_analyses.html#cell-type-composition) and the [feature-model routes](implementation/methods_differential_analyses.html#molecular-pseudobulk-analyses), including random effects, custom design and contrast functions and paired cell-type designs.

## Run

Preview the targets tagged for this module:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

## Review

Review the plots ([examples](gallery.md#differential-analyses)); each plot's subtitle and caption say what to look for.

```text
<store>/plots/my_aggregation/differential_analyses/
├── pseudobulk_depth_distribution_plot.png
├── cell_type_composition/model_plots_condition_abundance/
├── gene_expression/
│   ├── volcano_plots/condition_model/
│   ├── gene_set_enrichment/{Hallmark,Reactome}/enrichment_plots/condition_model/
│   ├── significant_elements_plot.png
│   └── p_value_distribution_plot.png
├── chromatin_accessibility/
├── motif_family_accessibility/
├── transcription_factor_activity/
├── significant_elements_modality_distribution_plots/
└── CollecTRI_JASPAR/activity_accessibility_concordance_plots/
```

The other feature families follow the `gene_expression/` layout without gene sets.

Read feature-model results in R, for example for gene expression:

```{.r filename="R"}
targets::tar_read(
  results_tibble.gene_expression.differential_analyses.my_aggregation
)
```

The run also writes cohort tables recording each model's included and excluded donors, and the counts and results of each abundance model, under `<store>/files/my_aggregation/differential_analyses/`.

Use [Troubleshooting](troubleshooting.md) if a formula, contrast or metadata join fails.

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=differential_analyses). The [differential analyses implementation page](implementation/methods_differential_analyses.html) describes the methods and their key fixed values, links to the source files and shows the target structure.

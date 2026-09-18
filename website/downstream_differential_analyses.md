# Differential analyses

```{r setup, include = FALSE}
pipeline_name <- "differential_analyses"
source("helpers/_setup.R")
```

## When to use this module

Use this module to ask how cell-type proportions, gene expression, or chromatin accessibility differ with a condition or donor phenotype. It requires biological replication across donors or samples. Molecular measurements are combined into **pseudobulks**: counts summarized for each cell type within a donor or sample.

The module does not create biological replication. The donor structure, covariates, design formula, and contrasts must be defensible for the intended analysis before the workflow is run.

See the [Differential analyses gallery](gallery_differential_analyses.md) for representative diagnostics and the [implementation graph](implementation/implementation_differential_analyses.html) for target structure.

## Prerequisites

Before enabling the module, confirm that:

- you have reviewed the main aggregation and its cell-type annotations;
- WNN cell-type metadata and GEX and ATAC pseudobulk matrices are available;
- the donor metadata contains one unique row per `donor_id` and every variable used in a model;
- model variables are donor- or pseudobulk-sample-level variables, not duplicated cell-level measurements; and
- the number and distribution of donors support the specified design and contrasts.

Use `differential_analyses_extended_donor_id_metadata_tsv` when the modelling table needs variables beyond the aggregation's normal donor metadata. It must retain the same unique `donor_id` key.

## Outputs

Choose the output that matches your question:

| Question | Output family |
|------------------------------------|------------------------------------|
| Do cell-type proportions differ? | `cell_type_composition` |
| Which genes change expression? | `gene_expression` |
| Which peaks change accessibility? | `chromatin_accessibility` |
| Which motif families change accessibility? | `motif_family_accessibility` (JASPAR) |
| Which regulators show altered expression-based activity? | `transcription_factor_activity` (CollecTRI) |

The module also produces model diagnostics, comparisons across modalities, and gene-set tests for Hallmark and Reactome pathways. Motif-family accessibility summarizes ATAC evidence; transcription-factor activity is inferred from gene expression using CollecTRI. Interpret each in the context of its measurement.

Plot directories use these descriptive family names below `plots/<aggregation>/differential_analyses/`. Gene-set plots appear under `gene_expression/gene_set_enrichment/Hallmark/enrichment_plots/<model>/` or the corresponding `Reactome` directory. Volcano outputs use `<family>/volcano_plots/<model>/`; saved plot targets omit redundant `_file` and `_files` suffixes. Renaming targets creates new cache entries and output paths on the next run; existing output directories are not migrated.

See the [method details](implementation/implementation_differential_analyses.html#method-details) for activity inference, motif-family definitions, and gene-set testing.

## Configure

Add `modules` to the existing aggregation entry, keeping its input and marker settings:

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

Then create a matching row directly in `configuration/cfg_module_differential_analyses.yaml`.

``` {.yaml filename="configuration/cfg_module_differential_analyses.yaml"}
your_aggregation:
  differential_analyses_cell_type_composition_models:
    condition_abundance:
      formula: ~ condition
      contrast_specs_vec:
        treated_vs_control: conditiontreated
      plot_phenotype_vars: condition
      color_by: condition
  differential_analyses_pseudobulk_models:
    condition_model:
      cell_type_subset: NULL
      design_matrix_func_name: NULL
      formula: ~ 0 + cluster + condition
      random_effect: NULL
      contrast_specs_vec:
        treated_vs_control: conditiontreated
```

Both branches use named models, donor eligibility checks and named contrasts. Abundance models use `differential_analyses_cell_type_composition_models`; feature models use `differential_analyses_pseudobulk_models`. Omitting abundance models disables that branch. The former aggregation-wide cell-type composition formula, phenotype and colour settings have been replaced by fields inside each named model.

For mixed tissues, set `GEM_well_IDs` inside an abundance model to define its population, for example the six left-ventricle wells. Optional `donor_ids` can further restrict donors in either branch. Donors with missing model metadata or no selected samples are excluded and recorded in model-specific cohort TSVs. Feature cohorts also report retained pseudobulk sample counts and depth-filter exclusions.

By default, abundance models test all cell-type labels observed in their eligible population, including unassigned labels. Optional `cell_types_to_test` restricts the response cell types **without changing the denominator**: every retained nucleus in the selected wells contributes to its donor's total. Zero donor–cell-type counts remain in the analysis. In contrast, feature-model `cell_type_subset` selects the cells represented by the pseudobulks. Feature matrices already pool wells within donors and cell types, so they cannot support a late `GEM_well_IDs` filter; the module rejects that field for feature models.

Abundance models fit a separate fixed-effects beta-binomial logit model per cell type. Use a one-sided predictor formula, `formula: ~ ...`, and named `contrast_specs_vec`. The implementation supplies the fixed response `cbind(n_nuclei, n_other_nuclei)`; two-sided formulas are rejected. Custom design/contrast functions and random-effects formulas are not supported in this branch. Contrast tables report log-odds effects, Wald uncertainty, donor counts and BH FDR across tested cell types within each model/contrast. Failed fits are explicitly marked non-estimable. Counts, cohort tables and plots use the same eligible donors and denominators. Plot-only variables do not exclude donors from the fit.

The module selection below requests both configured abundance and pseudobulk outputs. Formula terms and contrast coefficients must match columns produced by the model matrix. The two branches retain their distinct response construction and fitting methods; sharing configuration does not make their effect estimates interchangeable.

The model example assumes `condition` distinguishes treated and control donors. Check which group is the reference and what each model coefficient represents before using `conditiontreated` as a contrast. Replace the example formula and contrast to match your study.

## Run

Preview the selected module outputs before running them:

``` {.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

``` {.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".your_aggregation")
)
```

## Review

Open the configured model outputs listed above and the [differential gallery](gallery_differential_analyses.md). Interpretation and method details are included in the plot subtitles and captions.

Runtime depends on donors, cell types, models, contrasts, and gene-set analyses. Use [Troubleshooting](troubleshooting.md) if a formula, contrast, or metadata join fails.

## Parameter reference

The OLINK and bulk-RNA path fields are reserved optional integration inputs and are not consumed by the current public differential-analysis selection. Leave them `NULL` unless the corresponding integration is implemented in your downstream workflow.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_parameter_overview("differential_analyses")
```

The public demos leave this module disabled. Comparing one healthy PBMC donor with one lymphoma lymph-node donor cannot separate condition, donor, and tissue effects. Configure differential analyses for a design with biological replication.
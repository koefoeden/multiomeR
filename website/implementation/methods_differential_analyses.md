# Differential analyses

This chapter covers the optional `differential_analyses` module. The target structure is shown in the [differential analyses graph](implementation_differential_analyses.md), and the prerequisites and configuration in [Differential analyses](../downstream_differential_analyses.html). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

## Cell-type composition

{{< include _shared_methods/cell_type_composition.md >}}

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Input | Annotation column counted per donor | Fixed | `WNN_harmony_SNN_cluster_cell_type` | `model_data.cell_type_composition` |
| Population | GEM wells defining the population | Configurable |  | `GEM_well_IDs` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Population | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Response | Classes tested; denominators always use all retained nuclei | Configurable |  | `cell_types_to_test` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Model | Fixed-effects formula | Configurable |  | `formula` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Model | Family, link and response | Fixed | beta-binomial with logit link on `cbind(n_nuclei, n_other_nuclei)`; no dispersion or zero-inflation formula, so package defaults apply | `fit_cell_type_composition_model()`, `glmmTMB::glmmTMB()` |
| Contrasts | Named linear contrasts | Configurable |  | `contrast_specs_vec` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Plot | Phenotype panels and colour variable | Configurable |  | `plot_phenotype_vars`, `color_by` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |

## Molecular pseudobulk analyses

{{< include _shared_methods/pseudobulk_differential_analyses.md >}}

### Pseudobulk construction

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Summation | Grouping and inputs | Fixed | counts summed (`method = "sum"`) per `WNN_harmony_SNN_cluster_cell_type` and donor; the ATAC input is the peak-QC-filtered consensus-peak matrix | `get_BPCells_pseudobulk_matrix()`, `pseudobulk_counts_BPCells_matrix_dir.GEX`, `pseudobulk_counts_BPCells_matrix_dir.ATAC` |
| DTFA | Motif families | Fixed | 233 JASPAR 2026 CORE vertebrate familial root motifs; membership map of 1,019 motifs | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf`, `resources/JASPAR2026_vertebrate_motif_families.tsv` |
| DTFA | Scores | Fixed | betterChromVAR analytic deviations with `compute = c("deviations", "z")` on peaks with non-zero pseudobulk counts; z-scores column-centred and quantile-normalised with `limma::normalizeBetweenArrays(method = "quantile")`, switched on by the project helper argument `normalize = TRUE` | `get_pseudobulk_chromVAR_background_record()`, `compute_pseudobulk_chromVAR_deviation_SE()`, `get_pseudobulk_chromVAR_accessibility_matrix()`, `get_pseudobulk_motif_family_accessibility_matrix()` |
| DCTA | Regulon network | Fixed | CollecTRI from `https://rescued.omnipathdb.org/CollecTRI.csv`, SHA-256 `86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0`; 43,536 signed interactions from 1,189 regulators | `CollecTRI_human_network_csv`, `read_CollecTRI_human_network()` |
| DCTA | Expression preprocessing | Fixed | `edgeR::filterByExpr(group = cluster)`; TMM normalisation with `normLibSizes()`; `cpm(log = TRUE, prior.count = 2)` | `get_pseudobulk_CollecTRI_TF_activity_matrix()` |
| DCTA | Inference | Fixed | `decoupleR::run_ulm()` with at least 5 measured targets per regulator (`min_targets = 5`) | `get_pseudobulk_CollecTRI_TF_activity_matrix()` |

### Sample and feature filtering

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Donors | Variables that must be non-missing | Configurable |  | `formula` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Extended donor metadata table | Configurable |  | [`differential_analyses_extended_donor_id_metadata_tsv`](../parameters.html#differential_analyses_extended_donor_id_metadata_tsv) |
| Samples | Annotation-class subset | Configurable |  | `cell_type_subset` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Samples | Minimum ATAC depth, DTFA only | Configurable |  | [`differential_analyses_motif_family_accessibility_min_ATAC_counts`](../parameters.html#differential_analyses_motif_family_accessibility_min_ATAC_counts) |
| Count routes | Sample and feature filtering | Fixed | samples with zero counts removed; `filterByExpr(design = design_matrix)` and TMM `normLibSizes()`, otherwise package defaults | `fit_pseudobulk_feature_matrix_model()`, `fit_pseudobulk_cell_type_matrix()`, `edgeR::filterByExpr()`, `edgeR::normLibSizes()` |

### Model routes

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Design | Formula or custom design function | Configurable |  | `formula`, `design_matrix_func_name` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Design | Contrasts | Configurable |  | `contrast_specs_vec` and custom contrast functions inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Route | Count detection | Fixed | integer check on up to 10 × 10 randomly sampled entries | `is_count_matrix()` |
| Route | Correlation route trigger | Configurable |  | `random_effect` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| edgeR route | Dispersion and test | Fixed | `estimateDisp()` with package defaults; `glmQLFit(robust = TRUE)`; `glmQLFTest()` | `fit_pseudobulk_feature_matrix_model()`, `get_pseudobulk_feature_model_results()`, `edgeR::estimateDisp.DGEList()` |
| limma routes | Moderation | Fixed | `eBayes()` with package defaults, without `robust` or `trend` | `get_pseudobulk_feature_model_results()`, `get_pseudobulk_cell_type_contrast_statistics()`, `limma::eBayes()` |
| Paired route | Trigger and pairing | Configurable |  | `cell_type_formula`, `pairing_variable`, `correlation_block` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Paired route | Per-class count fit | Fixed | `edgeR::voomLmFit()` with `normalize.method = "none"` after the filtering and TMM normalisation above | `fit_pseudobulk_cell_type_matrix()` |
| Paired route | Residual correlation between classes | Fixed | shared donors, more than coefficients + 2 required; up to 2,000 evenly spaced common features; Fisher-z mean trimmed at 0.15 | `estimate_pseudobulk_cell_type_residual_correlations()` |
| Paired route | Cross-class test | Fixed | t from the two class estimates and their covariance; df is the smaller per-class total df | `get_pseudobulk_paired_cell_type_contrast_statistics()` |

### Multiplicity, significance and gene sets

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Significance | FDR and significant features | Fixed | BH across the tested features of each model and contrast; significant when FDR \< 0.05 and `logFC != 0` | `get_pseudobulk_feature_model_results()`, `get_pseudobulk_cell_type_model_results()`, `get_pseudobulk_differential_significant_elements_tibble()` |
| Top features | Features labelled and queried per contrast | Fixed | 40 with the smallest nominal P | `top_features_tibble` |
| Gene sets | Collections | Fixed | MSigDB Hallmark (`H`) and Reactome (`C2`, `CP:REACTOME`); human gene sets for `Homo_sapiens`, mouse for `Mus_musculus` | `gene_sets`, `get_msigdb_gene_sets()` |
| Gene sets | Test | Fixed | `limma::cameraPR(inter.gene.cor = 0.01)`; at least 10 tested genes per set; BH within contrast and collection | `get_gene_set_enrichment_results()` |
| Open Targets | Trait identifier; empty skips the query | Configurable |  | [`differential_analyses_pseudobulk_OT_GWAS_efo_id`](../parameters.html#differential_analyses_pseudobulk_OT_GWAS_efo_id) |
| Open Targets | Query | Fixed | Open Targets Platform GraphQL API (`https://api.platform.opentargets.org/api/v4/graphql`), queried at run time for the top gene-expression features only | `top_feature_open_targets_evidence_tibble`, `get_OT_GWAS_gene_evidence_tibble()` |

### Cross-modality comparison

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Crosswalk | Complex members | Fixed | AP1: FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND; NFKB: NFKB1, NFKB2, REL, RELA, RELB | `get_CollecTRI_JASPAR_family_map()` |
| Summary | Family statistics and concordance | Fixed | median CollecTRI and TF-expression t per family; Spearman correlation with the motif-family t when at least 3 families are mapped in a contrast | `get_CollecTRI_JASPAR_family_comparison_tibble()`, `get_CollecTRI_JASPAR_concordance_tibble()` |

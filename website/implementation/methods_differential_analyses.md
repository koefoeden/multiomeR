# Differential analyses

This chapter describes the optional `differential_analyses` module and lists every setting that determines its results, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [differential analyses graph](implementation_differential_analyses.md); the user-facing prerequisites and configuration are in [Differential analyses](../downstream_differential_analyses.html).

Donors, not nuclei, are the biological replicates. Both branches take their annotation classes from the WNN cell-type label carried in the final WNN metadata. The module retains its own metadata file and full-tibble targets, then projects a canonical analysis view containing only donors in the aggregation and the columns required by the configured models and composition plots. Rows are ordered by `donor_id` and non-key columns by name, so changes to unused columns or out-of-aggregation donors stop at this projection boundary.

## Cell-type composition

Nuclei are counted by donor and annotation class. Every observed class is tested unless `cell_types_to_test` restricts the response classes; missing donor–class combinations are completed with zero counts, so no class is dropped for being observed in few donors. Each class is fitted separately as a beta-binomial model with logit link on the two-column response of nuclei in the class versus nuclei in all other classes, using the configured one-sided formula. Random effects, custom design functions and two-sided formulas are rejected in this branch. Contrasts are named linear combinations of the fixed effects and are tested with normal Wald statistics. Fits that error, fail to converge or lack a positive-definite Hessian are marked non-estimable. Benjamini–Hochberg (BH) correction is applied across the tested classes within each model and contrast; there is no adjustment across contrasts or models.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Input | Annotation column counted per donor | Hardcoded: target literal | `WNN_harmony_SNN_cluster_cell_type` | `model_data.cell_type_composition` in `module_differential_analyses/setup_and_cell_type_composition_targets.R` |
| Population | GEM wells defining the population | Configurable |  | `GEM_well_IDs` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Population | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Response | Classes tested; denominators always use all retained nuclei | Configurable |  | `cell_types_to_test` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Response | Two-column response prepended to the formula | Hardcoded: inline literal | `cbind(n_nuclei, n_other_nuclei)` | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Model | Fixed-effects formula | Configurable |  | `formula` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Model | Family and link; no dispersion or zero-inflation formula | Hardcoded: inline literal | `glmmTMB::betabinomial(link = "logit")`, package defaults | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Model | Random effects, design functions, contrast functions, paired fields | Hardcoded: inline literal | rejected with an error | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Model | Design validity | Hardcoded: inline literal | full rank, finite, more donors than coefficients | `validate_differential_design()` in `R/differential_analysis_helpers.R` |
| Contrasts | Named linear contrasts | Configurable |  | `contrast_specs_vec` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Test | Wald test from the conditional covariance | Hardcoded: inline literal | two-sided normal | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Test | Non-estimable rule | Hardcoded: inline literal | fit error, non-zero convergence code, or non-positive-definite Hessian | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Multiplicity | BH across tested classes within model and contrast | Hardcoded: inline literal | `p.adjust(method = "BH")` | `R/differential_analysis_helpers.R` |
| Plot | Interval and significance colour | Hardcoded: inline literal | ±1.96 SE; FDR \< 0.05 | `R/differential_analysis_helpers.R` |
| Plot | Phenotype panels and colour variable | Configurable |  | `plot_phenotype_vars`, `color_by` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |

## Pseudobulk construction

GEX and ATAC counts are summed within each donor and annotation class in the main pipeline, so the module reuses the same pseudobulk targets as the compatibility export. Sample identifiers combine the class and the donor. Four feature matrices are tested:

- **Gene expression (DGE)**: the GEX pseudobulk count matrix.
- **Chromatin accessibility (DCA)**: the consensus-peak pseudobulk count matrix after peak-level QC.
- **Motif-family accessibility (DTFA)**: betterChromVAR deviations of the JASPAR 2026 familial root motifs computed on the pseudobulk ATAC counts. Peaks with zero pseudobulk counts are dropped before the background model. The z-scores are column-centred and quantile-normalised across samples, so this matrix is continuous.
- **Transcription-factor activity (DCTA)**: signed CollecTRI regulator activities inferred from the GEX pseudobulks with the `decoupleR` univariate linear model. Genes are filtered by expression across classes, library sizes are normalised and log-CPM values are computed before inference. Regulators need a minimum number of measured targets. The CollecTRI network is downloaded from the OmniPath rescue archive and accepted only when it matches the pinned checksum.

Motif families come from the official JASPAR 2026 CORE vertebrate clustering. Each of the 233 families is represented by its published root motif, which is scanned directly against the consensus peaks; individual member motifs are used only as family metadata. The same family-level matrix supports marker plots and the Seurat export.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Summation | Aggregation method | Hardcoded: inline literal | `BPCells::pseudobulk_matrix(method = "sum")` | `get_BPCells_pseudobulk_matrix()` in `R/pseudobulk_helpers.R` |
| Summation | Grouping column for GEX and ATAC | Hardcoded: target literal | `WNN_harmony_SNN_cluster_cell_type` | `pseudobulk_counts_matrix.GEX` and `.ATAC` in `extra_targets/general_aggregation_targets.R` |
| Summation | ATAC input matrix | Hardcoded: target literal | `peak_QC_filtered_BPCells_matrix.ATAC` | `extra_targets/general_aggregation_targets.R` |
| DTFA | Deviation computation | Hardcoded: inline literal | betterChromVAR `compute = c("deviations", "z")`, `normalize = TRUE` | `R/pseudobulk_helpers.R` |
| DTFA | Post-processing of z-scores | Hardcoded: inline literal | column centring, `limma::normalizeBetweenArrays(method = "quantile")` | `R/pseudobulk_helpers.R` |
| DTFA | Peaks dropped before background | Hardcoded: inline literal | zero pseudobulk row sum | `R/pseudobulk_helpers.R` |
| DTFA | Motif universe | Hardcoded: target literal | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf`, 233 families checked | `extra_targets/setup_targets.R`, `extra_targets/ATAC_targets.R` |
| DTFA | Family membership table | Hardcoded: target literal | `resources/JASPAR2026_vertebrate_motif_families.tsv`, 1,019 motifs | `extra_targets/setup_targets.R` |
| DCTA | Network source and checksum | Hardcoded: target literal | `https://rescued.omnipathdb.org/CollecTRI.csv`, SHA-256 `86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0` | `module_differential_analyses/shared_targets.R` |
| DCTA | Network validation | Hardcoded: inline literal | 43,536 interactions, 1,189 sources, `mor` in {−1, 1} | `R/TF_activity_helpers.R` |
| DCTA | Gene filter before inference | Hardcoded: inline literal | `edgeR::filterByExpr(group = cluster)` | `R/TF_activity_helpers.R` |
| DCTA | Normalisation and log-CPM | Hardcoded: inline literal | `normLibSizes()` (TMM), `cpm(log = TRUE, prior.count = 2)` | `R/TF_activity_helpers.R` |
| DCTA | Minimum measured targets per regulator | Hardcoded: helper default | `min_targets = 5` | `R/TF_activity_helpers.R`, not overridden in `module_differential_analyses/targets.R` |

## Sample and feature filtering

For each configured pseudobulk model, donors missing any variable in the formula, donors outside an optional `donor_ids` list, and samples outside an optional annotation-class subset are removed. The DTFA branch additionally removes samples below the configured minimum ATAC depth. The remaining steps depend on the model route, described next. For the routes that model counts, zero-depth samples are removed, features are filtered with the design-aware `edgeR::filterByExpr()`, and library sizes are normalised with `edgeR::normLibSizes()`.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Donors | Variables that must be non-missing | Configurable |  | `formula` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Extended donor metadata table | Configurable |  | [`differential_analyses_extended_donor_id_metadata_tsv`](../parameters.html#differential_analyses_extended_donor_id_metadata_tsv) |
| Samples | Annotation-class subset | Configurable |  | `cell_type_subset` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Samples | Minimum ATAC depth, DTFA only | Configurable |  | [`differential_analyses_motif_family_accessibility_min_ATAC_counts`](../parameters.html#differential_analyses_motif_family_accessibility_min_ATAC_counts) |
| Samples | Zero-depth removal, count routes | Hardcoded: inline literal | column sum \> 0 | `R/pseudobulk_helpers.R` |
| Features | Expression filter, count routes | Hardcoded: inline literal | `edgeR::filterByExpr(design = design_matrix)`, package defaults | `R/pseudobulk_helpers.R` |
| Features | Library-size normalisation, count routes | Hardcoded: inline literal | `edgeR::normLibSizes()`, TMM | `R/pseudobulk_helpers.R` |
| Residual df | Minimum residual degrees of freedom, checked before and after filtering | Hardcoded: inline literal | samples − coefficients \> 0 | `R/pseudobulk_helpers.R` |

## Model routes

Model matrices come from the configured formula or a custom design function, and contrasts from named linear expressions or custom contrast functions. Whether a matrix holds counts is decided by checking a random sample of entries for integers. The route is then selected as follows:

1.  A configured `random_effect` selects the correlation route: `limma::voom()` for counts, then one pass of `limma::duplicateCorrelation()` with the random effect as block, then `lmFit()` with the consensus correlation.
2.  Otherwise, a count matrix with more than one coefficient selects the edgeR route: `estimateDisp()` with package defaults, `glmQLFit(robust = TRUE)` and `glmQLFTest()`.
3.  Otherwise, the basic limma route: `lmFit()`, `contrasts.fit()` and `eBayes()` with package defaults. Note that a count matrix with a single coefficient takes this route on the raw counts, without the filtering and normalisation described above.

An optional paired-cell-type route is selected by a `cell_type_formula`. It fits each annotation class separately, permits at most one pseudobulk per pairing unit and class, and requires at least two classes. Count matrices are fitted per class with `edgeR::voomLmFit()` blocked on the correlation block; continuous matrices use `duplicateCorrelation()` and `lmFit()`. Residual correlations between classes are estimated on shared donors from the common features and combined into cross-class contrast tests. A class pair needs more shared donors than coefficients plus two.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Design | Formula or custom design function | Configurable |  | `formula`, `design_matrix_func_name` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Design | Contrasts | Configurable |  | `contrast_specs_vec` and custom contrast functions inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Design | Coefficient-name sanitising | Hardcoded: inline literal | `:` to `.`, `-` to `_` | `R/pseudobulk_helpers.R` |
| Route | Count detection | Hardcoded: inline literal | integer check on up to 10 × 10 sampled entries | `is_count_matrix()` in `R/data_transformations.R` |
| Route | Correlation route trigger | Configurable |  | `random_effect` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Correlation route | Fit | Hardcoded: inline literal | `limma::voom()` for counts, one `duplicateCorrelation()` pass, `lmFit(correlation, block)` | `R/pseudobulk_helpers.R` |
| edgeR route | Dispersion and fit | Hardcoded: inline literal | `estimateDisp()` defaults, `glmQLFit(robust = TRUE)`, `glmQLFTest()` | `R/pseudobulk_helpers.R` |
| limma routes | Moderation | Hardcoded: inline literal | `eBayes()` defaults, no `robust` or `trend` | `R/pseudobulk_helpers.R` |
| Paired route | Trigger and pairing | Configurable |  | `cell_type_formula`, `pairing_variable`, `correlation_block` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Paired route | One pseudobulk per pairing unit and class; at least two classes | Hardcoded: inline literal | error otherwise | `R/pseudobulk_helpers.R` |
| Paired route | Per-class count fit | Hardcoded: inline literal | `edgeR::voomLmFit(block, normalize.method = "none")` | `R/pseudobulk_helpers.R` |
| Paired route | Minimum shared donors per class pair | Hardcoded: inline literal | more than coefficients + 2 | `R/pseudobulk_helpers.R` |
| Paired route | Residual-correlation estimate | Hardcoded: helper default | up to 2,000 evenly spaced common features; Fisher-z trimmed mean, `trim = 0.15` | `R/pseudobulk_helpers.R` |
| Paired route | Cross-class statistic | Hardcoded: inline literal | t from the two class estimates and their covariance; df is the minimum per-class total df | `R/pseudobulk_helpers.R` |
| Paired route | Parallel class fits | Hardcoded: inline literal | up to 6 forks | `R/pseudobulk_helpers.R`; `cores_req = 6` in `module_differential_analyses/pseudobulk_differential_targets.R` |
| Failure | Conditions that fail the branch | Hardcoded: inline literal | no samples or features, residual df below one, invalid block or contrast, too few shared donors | `R/pseudobulk_helpers.R` |

## Multiplicity, significance and gene sets

For each model and contrast, BH FDR is calculated across the tested features. A feature is counted as significant when its FDR is below the threshold and its effect is non-zero. Differential-expression statistics are tested against the MSigDB Hallmark and Reactome collections with the competitive `cameraPR` test. The statistic is the moderated t where the route provides one and otherwise a signed normal quantile derived from the nominal P-value, which is the case for the edgeR route. A set must contain the minimum number of tested genes, and FDR is calculated within each contrast and collection. The optional Open Targets annotation queries the platform for the top features of the gene-expression branch when an EFO identifier is configured.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Multiplicity | FDR across features per contrast | Hardcoded: inline literal | `p.adjust(method = "BH")` | `R/pseudobulk_helpers.R` |
| Significance | FDR threshold for significant counts | Hardcoded: helper default | `FDR_threshold = 0.05`, plus `logFC != 0` | `R/pseudobulk_helpers.R`, not overridden in `module_differential_analyses/pseudobulk_differential_targets.R` |
| Significance | Threshold used in plots and comparisons | Hardcoded: inline literal | 0.05 | `R/pseudobulk_helpers.R`, `R/TF_activity_helpers.R` |
| Top features | Features labelled and queried per contrast | Hardcoded: target literal | `n = 40` by nominal P | `module_differential_analyses/pseudobulk_differential_targets.R` |
| Gene sets | Collections | Hardcoded: target literal | MSigDB `H`; `C2` with `CP:REACTOME` | `module_differential_analyses/targets.R` |
| Gene sets | Species mapping | Hardcoded: inline literal | `Homo_sapiens` to human, `Mus_musculus` to mouse | `R/pseudobulk_helpers.R` |
| Gene sets | msigdbr version | Hardcoded: environment pin | `r-msigdbr >=26.1.0,<27` | `pixi.toml` |
| Gene sets | Test | Hardcoded: inline literal | `limma::cameraPR(inter.gene.cor = 0.01)` | `R/pseudobulk_helpers.R` |
| Gene sets | Minimum tested genes per set | Hardcoded: helper default | `min_genes_per_set = 10` | `R/pseudobulk_helpers.R`, not overridden in `module_differential_analyses/gene_set_enrichment_targets.R` |
| Gene sets | Multiplicity | Hardcoded: inline literal | BH within contrast and collection | `R/pseudobulk_helpers.R` |
| Gene sets | Terms shown per plot | Hardcoded: inline literal | 25 smallest P | `R/pseudobulk_helpers.R` |
| Open Targets | Trait identifier; empty skips the query | Configurable |  | [`differential_analyses_pseudobulk_OT_GWAS_efo_id`](../parameters.html#differential_analyses_pseudobulk_OT_GWAS_efo_id) |
| Open Targets | Scope and endpoint | Hardcoded: inline literal | gene-expression branch only; `https://api.platform.opentargets.org/api/v4/graphql` | `module_differential_analyses/pseudobulk_differential_targets.R`, `R/pseudobulk_helpers.R` |

## Cross-modality comparison

The cross-modality fragment maps CollecTRI regulators to JASPAR families and compares model t-statistics, not raw activity scales. AP1 and NFKB stay intact as complex regulons during activity inference; their canonical members are used only to associate the complexes with motif families. Source-level tables retain TF expression as a third reference. Family-level summaries use the median CollecTRI regulator t-statistic and report whether any mapped source is FDR-significant.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Crosswalk | Complex members | Hardcoded: inline literal | AP1: FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND; NFKB: NFKB1, NFKB2, REL, RELA, RELB | `R/TF_activity_helpers.R` |
| Summary | Family-level statistics | Hardcoded: inline literal | median CollecTRI t, median TF-expression t, first motif-family t and FDR | `R/TF_activity_helpers.R` |
| Summary | Concordance | Hardcoded: inline literal | Spearman correlation, at least 3 mapped families per contrast | `R/TF_activity_helpers.R` |

## Diagnostics and resources

Diagnostic outputs report pseudobulk depth with the ATAC threshold marked, cohort tables per model, contrast support (samples, donors, paired donors, smallest group), P-value distributions, signed significant-feature counts and volcano plots.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Resources | Filtering, fitting and DCTA allocations | Hardcoded: target literal | filter 60 GB; fit 6 cores, 60 GB; DCTA 32 GB | `module_differential_analyses/pseudobulk_differential_targets.R`, `module_differential_analyses/targets.R` |
| Environment | Package pins | Hardcoded: environment pin | `bioconductor-edger >=4.8.2`, `bioconductor-limma >=3.66.0`, `r-glmmtmb >=1.1.14`, `bioconductor-decoupler >=2.16.0` | `pixi.toml` |
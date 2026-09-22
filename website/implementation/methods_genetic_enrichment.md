# Genetic enrichment

This chapter describes the optional `genetic_enrichment` module and lists every setting that determines its results, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [genetic enrichment graph](implementation_genetic_enrichment.md); the user-facing configuration is in [Genetic enrichment](../downstream_genetic_enrichment.html). The sparse SCAVENGE reimplementation and its validation are described in [Algorithmic implementations](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

## GWAS inputs and peak weights

Study identifiers of the Open Targets form are resolved against the pinned platform release, whose study, credible-set, credible-set evidence and target datasets are downloaded once. Any other identifier is treated as a local Parquet file that must satisfy the local schema. With automatic fine-mapping selection, the first available method in a fixed priority order is used. Variants with posterior probability above the configured cutoff are retained and mapped to the consensus peaks that carry chromVAR state. Variant weights mapping to the same peak are summed and capped at one.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Studies | Study label, category, source identifier, fine-mapping method | Configurable | | fields inside [`genetic_enrichment_GWAS_studies`](../parameters.html#genetic_enrichment_GWAS_studies) |
| Open Targets | Platform release | Hardcoded: target literal | `26.03` | `module_genetic_enrichment/shared_targets.R` |
| Open Targets | Datasets | Hardcoded: target literal | `credible_set`, `study`, `evidence_gwas_credible_sets`, `target` | `module_genetic_enrichment/shared_targets.R` |
| Open Targets | Source classification | Hardcoded: inline literal | `^GCST[0-9]+$` is Open Targets; anything else is a local path | `R/GWAS_chromVAR_input_helpers.R` |
| Open Targets | Study type | Hardcoded: inline literal | `studyType == "gwas"` | `R/GWAS_chromVAR_input_helpers.R` |
| Fine-mapping | Automatic priority | Hardcoded: inline literal | SuSie, SuSiE-inf, PICS | `resolve_open_targets_GWAS_input_tibble()` in `R/GWAS_chromVAR_input_helpers.R` |
| Credible sets | Probability | Hardcoded: inline literal | 0.95 | `R/GWAS_chromVAR_input_helpers.R` |
| Credible sets | Chromosomes | Hardcoded: inline literal | autosomes, X, Y, MT | `R/GWAS_chromVAR_input_helpers.R` |
| Local files | Validator | Hardcoded: inline literal | 29 required columns, schema version 1, GRCh38, SHA-256 field | `validate_local_finemapped_GWAS_tibble()` in `R/GWAS_chromVAR_input_helpers.R` |
| Variants | Posterior probability cutoff | Configurable | | [`genetic_enrichment_posterior_probability_cutoff`](../parameters.html#genetic_enrichment_posterior_probability_cutoff) |
| Variants | Comparison | Hardcoded: inline literal | strictly greater than the cutoff | `R/GWAS_chromVAR_helpers.R` |
| Peak weights | Peak set | Hardcoded: target literal | rows of the ATAC chromVAR object | `module_genetic_enrichment/gchromVAR_targets.R` |
| Peak weights | Combination | Hardcoded: helper default | `weight_transform = "cap_1"`: sum, capped at 1 | `R/GWAS_chromVAR_helpers.R` |

## Nucleus-level deviations

The trait peak weights form a chromVAR annotation. Deviations and z-scores per nucleus are computed analytically with betterChromVAR on the ATAC chromVAR object reused from the main pipeline, with the GC-bias background described in [Cell-type annotation and motif accessibility](methods_annotation_and_motifs.md#motif-families-and-motif-accessibility).

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Object | chromVAR object and chunk context | Hardcoded: target literal | `chromVAR_obj.ATAC`, `chromVAR_chunk_context_records.ATAC` | `module_genetic_enrichment/targets.R` |
| Deviations | Method | Hardcoded: inline literal | `computeDeviationsAnalytic(compute = "z")` | `R/GWAS_chromVAR_helpers.R` |
| Background | Bias and bins | Hardcoded: environment pin | `addGCBias()`; `getBackgroundBins()` defaults | `R/celltype_labeling_helpers.R` |

## Annotation-class pseudobulk deviations

Annotation-class deviations are calculated separately rather than by averaging nucleus-level results. ATAC counts are summed by the GEX-derived cell-type label carried into the final WNN metadata, a betterChromVAR background is fitted to the peak-by-class matrix, and the raw deviation is the observed-minus-background accessibility of the weighted peaks relative to their expected accessibility. A relative deviation standardizes the raw deviations across classes within each trait. The analytic z-score uses the background variance; one-sided P-values and BH-adjusted values over the whole table are retained. Plot labels use unadjusted z thresholds. An absolute-effect branch weights variants by posterior probability times effect size when effect sizes are available.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Pseudobulks | Grouping column | Hardcoded: target literal | `PCA_harmony_SNN_cluster_cell_type` from the WNN metadata | `module_genetic_enrichment/GWAS_chromVAR_cell_type_targets.R` |
| Pseudobulks | Aggregation | Hardcoded: inline literal | `BPCells::pseudobulk_matrix(method = "sum")`, 6 threads | `R/pseudobulk_helpers.R` |
| Pseudobulks | Peak filter | Hardcoded: inline literal | zero-count peaks removed before the background | `R/pseudobulk_helpers.R` |
| Deviations | Raw, relative and z | Hardcoded: inline literal | weighted observed minus background over expected; standardized within trait; z from background variance | `R/GWAS_chromVAR_contribution_helpers.R` |
| Deviations | P-values and adjustment | Hardcoded: inline literal | one-sided normal; BH over all trait × class rows | `R/GWAS_chromVAR_contribution_helpers.R` |
| Plots | Support labels | Hardcoded: inline literal | `**` for z ≥ 2.326, `*` for z ≥ 1.645, unadjusted | `R/GWAS_chromVAR_contribution_helpers.R` |
| Plots | Compartment grouping | Configurable | | [`genetic_enrichment_compartment_patterns`](../parameters.html#genetic_enrichment_compartment_patterns) |
| Plots | Unmatched classes and ordering | Hardcoded: inline literal | `Other`; hierarchical clustering within compartments of more than two classes | `R/GWAS_plot_helpers.R` |
| Absolute effect | Eligibility | Hardcoded: inline literal | variant-level effects when all variants have them, else locus-level, else skipped | `R/GWAS_chromVAR_absolute_effect_helpers.R` |

## Locus attribution and detail plots

Per-class contributions are attributed to loci and variants and reconciled against the class totals. Loci are labelled with the highest-scoring locus-to-gene (L2G) genes from the Open Targets evidence. Detail plots are drawn for the top loci of classes passing a configurable z screen.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| L2G | Score filter and label | Hardcoded: inline literal | score ≥ 0.05; top 3 genes | `R/GWAS_chromVAR_contribution_helpers.R` |
| Reconciliation | Tolerances | Hardcoded: inline literal | 1e-10 peak–variant; 1e-8 level sums | `R/GWAS_chromVAR_contribution_helpers.R` |
| Detail plots | Minimum z | Configurable | | [`genetic_enrichment_variant_detail_min_z`](../parameters.html#genetic_enrichment_variant_detail_min_z) |
| Detail plots | Loci, flank, coverage | Hardcoded: helper default | 3 loci per class; 25 kb flank; 500 bins; 0.999 coverage cap | `R/GWAS_chromVAR_contribution_helpers.R` |
| Attribution plots | Loci shown | Hardcoded: target literal | heatmaps 15; bar plots 5 | `module_genetic_enrichment/GWAS_chromVAR_contribution_targets.R` |

## SCAVENGE trait-relevance propagation

Trait-relevance scores are computed from the nucleus-level z-scores on the WNN SNN graph only. Nuclei whose one-sided normal-tail probability is at or below the seed cutoff are seeds, capped at the configured fraction; when no nucleus qualifies, all scores are zero. The graph's nonzero support becomes binary adjacency, degree-zero nuclei are removed, and seed mass is propagated by a random walk with restart until the L1 change is below the tolerance. Degree-matched seed permutations are sampled sequentially in R and evaluated by a native worker that streams per-nucleus exceedance counts. Cell-level empirical P-values are the exceedance fraction; cluster-level tests compare medians with add-one P-values and BH within each grouping. Scores are capped at an upper quantile, min–max scaled and multiplied by the mean z of the top nuclei.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Graph | Representation | Hardcoded: target literal | `WNN_harmony_SNN` only | `module_genetic_enrichment/targets.R` |
| Graph | Binarization | Hardcoded: inline literal | nonzero support set to 1; column-normalized transition matrix | `R/SCAVENGE_helpers.R` |
| Seeds | P cutoff | Hardcoded: helper default | `p_value_cutoff = 0.05` | `R/SCAVENGE_helpers.R` |
| Seeds | Maximum seed fraction | Configurable | | [`genetic_enrichment_SCAVENGE_seed_percent`](../parameters.html#genetic_enrichment_SCAVENGE_seed_percent) |
| Seeds | z pre-filter | Hardcoded: helper default | non-finite dropped; `max_z_score = 1000` | `R/SCAVENGE_helpers.R` |
| Walk | Restart probability | Configurable | | [`genetic_enrichment_SCAVENGE_restart_prob`](../parameters.html#genetic_enrichment_SCAVENGE_restart_prob) |
| Walk | Convergence | Hardcoded: helper default | L1 change ≤ 1e-5; at most 10,000 iterations | `R/SCAVENGE_helpers.R`, `src/scavenge_random_walk.cpp` |
| Permutations | Count | Configurable | | [`genetic_enrichment_SCAVENGE_permutation_times`](../parameters.html#genetic_enrichment_SCAVENGE_permutation_times) |
| Permutations | Cores and chunks | Hardcoded: target literal | 15 cores; 4 chunks per core | `module_genetic_enrichment/SCAVENGE_graph_targets.R` |
| Cell P | Definition and call | Hardcoded: inline literal | exceedances (strictly greater) / permutations; significant at P ≤ 0.05 | `R/SCAVENGE_helpers.R` |
| Cluster P | Definition | Hardcoded: inline literal | median score; (exceedances + 1)/(permutations + 1) with ≥; BH within grouping | `R/SCAVENGE_helpers.R` |
| Score | Cap, scaling and scale factor | Hardcoded: helper default | cap at the 0.95 quantile; min–max; × mean z of the top 1 % | `R/SCAVENGE_helpers.R` |
| Summaries | Groupings | Hardcoded: helper default | `WNN_harmony_SNN_cluster_named`, `WNN_harmony_SNN_cluster_cell_type` | `R/SCAVENGE_helpers.R` |
| Summaries | Statistics | Hardcoded: inline literal | count, significant count and proportion, median, mean, quartiles, range, null median and 0.95 quantile | `R/SCAVENGE_helpers.R` |
| Plots | Heatmap labels | Hardcoded: inline literal | on BH-adjusted cluster P: `***` ≤ 0.001, `**` ≤ 0.01, `*` ≤ 0.05 | `R/GWAS_plot_helpers.R` |

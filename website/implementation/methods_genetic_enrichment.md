# Genetic enrichment

This chapter covers the optional `genetic_enrichment` module. The target structure is shown in the [genetic enrichment graph](implementation_genetic_enrichment.md), and the configuration in [Genetic enrichment](../downstream_genetic_enrichment.html). The sparse SCAVENGE implementation is compared with its reference in [Algorithmic implementations](algorithm_validation.md#sparse-scavenge-propagation-and-significance). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

{{< include _shared_methods/genetic_enrichment.md >}}

## GWAS inputs and peak weights

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Studies | Study label, category, source identifier, fine-mapping method | Configurable |  | fields inside [`genetic_enrichment_GWAS_studies`](../parameters.html#genetic_enrichment_GWAS_studies) |
| Open Targets | Platform release and datasets | Fixed | release `26.03`: `credible_set`, `study`, `evidence_gwas_credible_sets`, `target` | `open_targets_credible_set_dataset_path`, `open_targets_study_dataset_path`, `open_targets_gwas_credible_sets_evidence_dataset_path`, `open_targets_target_dataset_path` |
| Open Targets | Study selection | Fixed | identifiers matching `^GCST[0-9]+$`, anything else is a local file; `studyType == "gwas"` | `classify_GWAS_source()`, `resolve_open_targets_GWAS_input_tibble()` |
| Fine-mapping | Automatic priority | Fixed | SuSie, SuSiE-inf, PICS | `resolve_open_targets_GWAS_input_tibble()` |
| Credible sets | Open Targets variants | Fixed | every variant in each credible-set `locus`, without an `is95CredibleSet` filter, on autosomes, X, Y or MT | `get_open_targets_credible_set_variants_tibble()` |
| Variants | Posterior probability cutoff | Configurable |  | [`genetic_enrichment_posterior_probability_cutoff`](../parameters.html#genetic_enrichment_posterior_probability_cutoff) |
| Variants | Cutoff comparison | Fixed | posterior probability strictly greater than the cutoff | `filter_credible_set_variants()` |
| Peak weights | Peaks and combination | Fixed | peaks of the ATAC chromVAR object; weights of variants in the same peak summed and capped at 1 (`weight_transform = "cap_1"`) | `genetic_enrichment_peak_ranges`, `get_GWAS_chromVAR_peak_weight_record()` |

## Nucleus-level deviations

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Deviations | Nucleus-level statistic | Fixed | analytic z-scores only (`compute = "z"`) | `get_GWAS_chromVAR_z_score_chunk_record()` |

## Annotation-class pseudobulk deviations

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Pseudobulks | Grouping | Fixed | ATAC counts summed (`method = "sum"`) per `PCA_harmony_SNN_cluster_cell_type` from the WNN metadata | `cell_type_pseudobulk_counts_BPCells_matrix_dir.ATAC`, `get_BPCells_group_pseudobulk_matrix()` |
| Pseudobulks | Peak filter | Fixed | peaks with zero pseudobulk counts removed before the background | `get_pseudobulk_chromVAR_background_record()` |
| Plots | Support labels | Fixed | `**` for z ≥ 2.326, `*` for z ≥ 1.645, unadjusted | `chromVAR_Z_support_labels()` |
| Plots | Compartment grouping | Configurable |  | [`genetic_enrichment_compartment_patterns`](../parameters.html#genetic_enrichment_compartment_patterns) |
| Absolute effect | Eligibility | Fixed | variant-level effects when every variant has one, otherwise locus-level effects when every locus has one; otherwise skipped | `infer_GWAS_absolute_effect_weighting()` |

## Locus attribution and detail plots

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| L2G | Gene labels | Fixed | L2G score ≥ 0.05; top 3 genes per locus | `get_open_targets_GWAS_locus_to_gene_tibble()`, `get_GWAS_chromVAR_locus_contribution_tibble()` |
| Detail plots | Minimum z | Configurable |  | [`genetic_enrichment_variant_detail_min_z`](../parameters.html#genetic_enrichment_variant_detail_min_z) |
| Detail plots | Loci | Fixed | classes with a positive deviation passing the z screen; per class the top 3 loci by absolute contribution plus the top 3 by combined contribution and effect-size percentile; 25 kb flank | `prepare_GWAS_variant_contribution_detail_records()`, `select_GWAS_detail_loci()` |
| Attribution plots | Loci shown individually | Fixed | 15 in heatmaps; 5 in bar plots | `chromVAR_locus_contribution_per_GWAS_heatmaps.cell_type_pseudobulk`, `chromVAR_locus_contribution_per_GWAS_faceted_bars_plots.cell_type_pseudobulk` |

## SCAVENGE trait-relevance propagation

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Graph | Representation | Fixed | WNN SNN graph (`WNN_harmony_SNN`) only | `graph_matrix` |
| Seeds | Selection | Fixed | one-sided normal P ≤ 0.05; nuclei with non-finite z or z \> 1000 excluded first | `get_SCAVENGE_seed_index()`, `get_SCAVENGE_result_from_chromVAR_z_score_record()` |
| Seeds | Maximum seed fraction | Configurable |  | [`genetic_enrichment_SCAVENGE_seed_percent`](../parameters.html#genetic_enrichment_SCAVENGE_seed_percent) |
| Walk | Restart probability | Configurable |  | [`genetic_enrichment_SCAVENGE_restart_prob`](../parameters.html#genetic_enrichment_SCAVENGE_restart_prob) |
| Walk | Convergence | Fixed | L1 change ≤ 1e-5; at most 10,000 iterations | `run_sparse_random_walk_with_restart()`, `run_SCAVENGE_permutation_statistics()` |
| Permutations | Count | Configurable |  | [`genetic_enrichment_SCAVENGE_permutation_times`](../parameters.html#genetic_enrichment_SCAVENGE_permutation_times) |
| Significance | Cell and cluster P-values | Fixed | cell: strictly greater permuted scores / permutations, significant at P ≤ 0.05; cluster: median score, (exceedances + 1) / (permutations + 1) counting ties, BH within each grouping | `get_SCAVENGE_result_from_chromVAR_z_score_record()`, `summarize_SCAVENGE_cluster_permutations()` |
| Score | Cap and scale | Fixed | capped at the 0.95 quantile, min–max scaled, multiplied by the mean z of the top 1 % of nuclei | `get_SCAVENGE_result_from_chromVAR_z_score_record()`, `get_SCAVENGE_scale_factor()` |
| Summaries | Groupings | Fixed | `WNN_harmony_SNN_cluster_named`, `WNN_harmony_SNN_cluster_cell_type` | `summarize_SCAVENGE_TRS_by_groups()`, `get_SCAVENGE_cluster_index_record()` |
| Plots | Heatmap stars | Fixed | on BH-adjusted cluster P: `***` ≤ 0.001, `**` ≤ 0.01, `*` ≤ 0.05 | `add_SCAVENGE_heatmap_significance()` |

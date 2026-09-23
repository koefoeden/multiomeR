# Peak–gene correlation

This chapter covers the optional `peak_gene_correlation` module. The target structure is shown in the [peak–gene correlation graph](implementation_peak_gene_correlation.md), and the configuration in [Peak–gene correlation](../downstream_peak_gene_correlation.html). The compiled kernels are compared with their references in [Algorithmic implementations](algorithm_validation.md). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

{{< include _shared_methods/peak_gene_correlation.md >}}

## Cell groups, candidate pairs and donor–state pseudobulks

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cell groups | Annotation column and minimum class size | Fixed | `WNN_harmony_SNN_cluster_cell_type`; at least 200 nuclei | `make_peak_gene_correlation_cell_groups()` |
| Pairs | TSS and search window | Fixed | TSS at the gene start, or the gene end on the minus strand; peak centre within 250 kb of the TSS | `make_peak_gene_correlation_gene_TSS_tibble()`, `make_peak_gene_correlation_candidate_pairs()` |
| Pairs | Self-promoter window and link classes | Fixed | strand-aware −1,500 to +500 bp around the TSS; self-promoter, gene body, proximal (≤ 10 kb), distal, in that precedence | `make_peak_gene_correlation_candidate_pairs()` |
| Donors | Donor and depth columns | Fixed | `donor_id`; first available of `nCount_RNA`, `gex_umis_count` and of `nCount_ATAC`, `atac_fragments` | `make_peak_gene_correlation_donor_state_record()` |
| States | Embedding | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs), [`aggregation_harmony_correction_metadata_col_names`](../parameters.html#aggregation_harmony_correction_metadata_col_names), [`aggregation_extra_harmony_covars_ATAC`](../parameters.html#aggregation_extra_harmony_covars_ATAC) |
| States | Dimensions and scaling | Fixed | LSI dimensions 2–20 present in the embedding; zero-variance dimensions dropped and the rest standardized within the class | `make_peak_gene_correlation_donor_state_record()` |
| States | Number of bins | Fixed | `min(20, floor(median nuclei per eligible donor / 20))`; classes with fewer than 2 bins are skipped | `make_peak_gene_correlation_donor_state_record()` |
| States | k-means | Fixed | Lloyd algorithm, `iter.max = 1000`, `nstart = 1`, seed 1 | `make_peak_gene_correlation_donor_state_record()` |
| Pseudobulks | Minimum sizes | Fixed | donors need at least 40 nuclei before binning; at least 20 nuclei per donor–state, 2 states per donor, and 1 donor per class and per state | `make_peak_gene_correlation_donor_state_record()` |

## Normalization, eligibility and measurement-support filtering

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Normalization | Scaling | Fixed | counts per million (`scale_factor = 1e6`) of the metadata-derived pseudobulk depth, then `log1p` | `normalize_peak_gene_correlation_aggregate_matrices()` |
| Eligibility | Detection | Fixed | genes detected and peaks accessible in at least 5 % of pseudobulks | `prepare_peak_gene_correlation_branch()` |
| Eligibility | Branch size | Fixed | at least 10 pseudobulks; pseudobulks − design rank − 1 ≥ 5 | `prepare_peak_gene_correlation_branch()` |
| Nuisance design | Covariates | Fixed | donor fixed effects; standardized `log1p` GEX and ATAC depth when non-constant; QR-pruned to full rank | `make_peak_gene_correlation_design_matrix()` |
| Support filter | Preset | Configurable |  | [`peak_gene_correlation_filter`](../parameters.html#peak_gene_correlation_filter) |
| Support filter | Preset thresholds | Fixed | RNA count, ATAC count, supporting aggregates, shared donors, supporting aggregates per donor. Lenient: 5, 3, `max(6, 10%)`, 2, 2; moderate: 10, 5, `max(6, 10%)`, 2, 2; strict: 10, 5, `max(10, 20%)`, 3, 3 | `peak_gene_filter_settings()` |
| Support filter | Depth scaling of count thresholds | Fixed | `max(2, count × aggregate depth / median depth)` | `filter_peak_gene_candidate_pairs()` |

## Conditional and hierarchical tests

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| HC3 | Test distribution | Fixed | two-sided t with pseudobulks − design rank − 1 degrees of freedom | `score_peak_gene_correlations_for_cell_group()` |
| Multiplicity | BH family | Fixed | within each class; HC3 over the non-missing P-values; hierarchical over all eligible pairs, including unreliable fits | `finalize_peak_gene_correlation_results()`, `finalize_peak_gene_hierarchical_results()` |
| Links | Conditional link | Fixed | correlation ≥ 0.15, FDR \< 0.05, not self-promoter | `make_peak_gene_correlation_links()` |
| Hierarchical | Donor requirement | Fixed | at least 2 donors with within-donor peak variation | `score_peak_gene_hierarchical_associations()` |
| Hierarchical | Unreliable fits | Fixed | kernel diagnostic raised or Kenward–Roger df \< 1: estimate kept, P-value missing | `score_peak_gene_hierarchical_associations()` |

## Prioritization, top links and plots

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| SuSiE | Genes and candidate peaks | Fixed | genes with a conditional link; up to 50 genes per branch by best FDR and 500 peaks per gene by absolute correlation; at least 2 variable peaks | `finemap_peak_gene_correlations_for_branch()` |
| SuSiE | Model settings | Fixed | `L = min(10, n_peaks)`, `intercept = FALSE`, `standardize = TRUE`, `estimate_residual_variance = TRUE`, `max_iter = 100`; credible-set coverage 0.95 | `finemap_peak_gene_correlations_for_branch()` |
| SuSiE | Records retained | Fixed | PIP ≥ 0.01 or credible-set member, otherwise the top peak | `finemap_peak_gene_correlations_for_branch()` |
| Top links | Links per cell group | Configurable |  | [`peak_gene_correlation_top_links_per_cell_group`](../parameters.html#peak_gene_correlation_top_links_per_cell_group) |

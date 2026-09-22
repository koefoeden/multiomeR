# Peak–gene correlation

This chapter describes the optional `peak_gene_correlation` module and lists every setting that determines its results, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [peak–gene correlation graph](implementation_peak_gene_correlation.md); the user-facing configuration is in [Peak–gene correlation](../downstream_peak_gene_correlation.html). Only two manifest parameters belong to this module; everything else is fixed. The module also depends indirectly on the ATAC embedding settings of the aggregation.

## Cell groups, candidate pairs and donor–state pseudobulks

The analysis runs separately within WNN annotation classes with enough nuclei. Consensus peaks are paired with gene transcription start sites (TSSs) on the same chromosome within a fixed window, measured from the peak centre; each pair is classified as self-promoter, gene-body, proximal or distal. Within a class, nuclei from eligible donors are partitioned into mutually exclusive ATAC-state bins by k-means on the scaled ATAC LSI or Harmony dimensions; the number of bins adapts to the median number of nuclei per eligible donor. Counts are summed within each donor–state combination. Pseudobulks below the minimum size are dropped, and donors and states are pruned iteratively until every donor contributes the minimum number of states. Per-nucleus GEX and ATAC depth columns are required; groups lacking donor or depth metadata are skipped with a diagnostic.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cell groups | Annotation column | Hardcoded: helper default | `WNN_harmony_SNN_cluster_cell_type` | `make_peak_gene_correlation_cell_groups()` in `R/peak_gene_correlation_helpers.R` |
| Cell groups | Minimum nuclei per class | Hardcoded: helper default | 200 | `R/peak_gene_correlation_helpers.R` |
| Pairs | TSS definition | Hardcoded: inline literal | gene start, or end on the minus strand | `R/peak_gene_correlation_helpers.R` |
| Pairs | Maximum peak-centre to TSS distance | Hardcoded: helper default | `max_distance = 250000` | `R/peak_gene_correlation_helpers.R` |
| Pairs | Self-promoter window | Hardcoded: inline literal | strand-aware −1,500 to +500 bp around the TSS | `R/peak_gene_correlation_helpers.R` |
| Pairs | Link classes | Hardcoded: inline literal | self-promoter, gene body, proximal (≤ 10 kb), distal, in that precedence | `R/peak_gene_correlation_helpers.R` |
| Donors | Donor column | Hardcoded: helper default | `donor_id` | `make_peak_gene_correlation_donor_state_record()` in `R/peak_gene_correlation_helpers.R` |
| Donors | Depth columns | Hardcoded: inline literal | first of `nCount_RNA`, `gex_umis_count`; first of `nCount_ATAC`, `atac_fragments` | `R/peak_gene_correlation_helpers.R` |
| States | Embedding source | Configurable |  | ATAC LSI or Harmony dimensions from [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) and the Harmony covariates |
| States | Dimensions used | Hardcoded: helper default | `2:20`, intersected with available dimensions | `R/peak_gene_correlation_helpers.R` |
| States | Preprocessing | Hardcoded: inline literal | zero-variance dimensions dropped; columns standardized within the group | `R/peak_gene_correlation_helpers.R` |
| States | Maximum bins and bin rule | Hardcoded: helper default | 20; `min(20, floor(median nuclei per eligible donor / 20))`; fewer than 2 skips the group | `R/peak_gene_correlation_helpers.R` |
| States | k-means | Hardcoded: inline literal | Lloyd, `iter.max = 1000`, `nstart = 1`, seed 1; non-convergence is an error | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Minimum nuclei per donor–state | Hardcoded: helper default | 20 | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Minimum states per donor | Hardcoded: helper default | 2 | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Donor eligibility before binning | Hardcoded: inline literal | at least 40 nuclei | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Minimum donors per group and per state | Hardcoded: helper default | 1 and 1 | `R/peak_gene_correlation_helpers.R` |

## Normalization, eligibility and measurement-support filtering

GEX and ATAC pseudobulk counts are separately scaled to counts per million using the metadata-derived pseudobulk depth and log1p-transformed. A chromosome branch requires a minimum number of pseudobulks and residual degrees of freedom after the nuisance design, and genes and peaks must be detected in a minimum fraction of pseudobulks. The configurable measurement-support filter then removes hypotheses whose gene and peak are not both supported in enough shared donors; count thresholds scale with each aggregate's depth relative to the class median. Excluded hypotheses never enter the BH family.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Normalization | Scale factor and transform | Hardcoded: helper default | `scale_factor = 1e6`, `log1p` | `R/peak_gene_correlation_helpers.R` |
| Eligibility | Gene detection fraction | Hardcoded: helper default | 0.05 | `prepare_peak_gene_correlation_branch()` in `R/peak_gene_correlation_helpers.R` |
| Eligibility | Peak accessibility fraction | Hardcoded: helper default | 0.05 | `R/peak_gene_correlation_helpers.R` |
| Eligibility | Minimum pseudobulks per branch | Hardcoded: helper default | 10 | `R/peak_gene_correlation_helpers.R` |
| Eligibility | Minimum residual df | Hardcoded: inline literal | pseudobulks − design rank − 1 ≥ 5 | `R/peak_gene_correlation_helpers.R` |
| Nuisance design | Covariates | Hardcoded: inline literal | donor fixed effects; standardized `log1p` GEX and ATAC depth when non-constant; QR-pruned to full rank | `make_peak_gene_correlation_design_matrix()` in `R/peak_gene_correlation_helpers.R` |
| Support filter | Preset | Configurable |  | [`peak_gene_correlation_filter`](../parameters.html#peak_gene_correlation_filter) |
| Support filter | Preset thresholds | Hardcoded: inline literal | lenient: RNA 5, ATAC 3 counts, `max(6, 10%)` aggregates, 2 shared donors, 2 aggregates per donor; moderate: 10, 5, `max(6, 10%)`, 2, 2; strict: 10, 5, `max(10, 20%)`, 3, 3 | `peak_gene_filter_settings()` in `R/peak_gene_filter_helpers.R` |
| Support filter | Depth scaling of count thresholds | Hardcoded: inline literal | `max(2, count × depth / median depth)` | `R/peak_gene_filter_helpers.R` |

## Conditional and hierarchical tests

The conditional analysis residualizes GEX and ATAC values against the nuisance design and reports the Pearson correlation of the residuals with an HC3 heteroskedasticity-robust regression test. The hierarchical analysis fits donor fixed intercepts, the depth covariates and a donor-varying slope for the within-donor-centred peak value, using project-owned compiled kernels for profiled restricted maximum likelihood and Kenward–Roger inference; the lme4 and pbkrtest route exists only as a test reference. It requires within-donor peak variation in at least two donors. Fits that fail the kernel diagnostics keep their estimates but no inferential P-value. BH correction runs within each annotation class over the complete eligible pair family, counting unreliable tests. Because every support preset requires at least two shared donors, single-donor classes yield diagnostics but no tests.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| HC3 | Correlation and test | Hardcoded: inline literal | Pearson on residuals; OLS slope with HC3 sandwich SE; two-sided t on residual df | `R/peak_gene_correlation_helpers.R` |
| HC3 | Numerical guards | Hardcoded: inline literal | leverage ≥ 1 − 1e-8 gives NA; residual-variation tolerance 1e-12 | `R/peak_gene_correlation_helpers.R` |
| HC3 | Multiplicity | Hardcoded: inline literal | BH by class over non-missing P-values | `R/peak_gene_correlation_helpers.R` |
| HC3 | Link definition for summary outputs | Hardcoded: inline literal | correlation ≥ 0.15, FDR \< 0.05, not self-promoter | `make_peak_gene_correlation_links()` in `R/peak_gene_correlation_helpers.R` |
| Hierarchical | Model | Hardcoded: inline literal | `y ~ 0 + design + x + (0 + x | donor)`, REML, x within-donor centred | `score_peak_gene_hierarchical_associations()` in `R/peak_gene_hierarchical_helpers.R` |
| Hierarchical | Engine | Hardcoded: inline literal | `src/peak_gene_REML.cpp`, `src/peak_gene_KR.cpp` (pbkrtest 0.5.5 equations) | `R/peak_gene_hierarchical_helpers.R` |
| Hierarchical | Minimum donors and within-donor variation | Hardcoded: inline literal | ≥ 2 donors; ≥ 2 donors with within-donor peak variation | `R/peak_gene_hierarchical_helpers.R` |
| Hierarchical | Variance-ratio search | Hardcoded: inline literal | grid over `expm1(0..32)` then golden section, tolerance 1e-9; zero allowed | `src/peak_gene_REML.cpp` |
| Hierarchical | Unreliable-fit rules | Hardcoded: inline literal | kernel diagnostics (no residual variation, ratio out of range, conditioning \< 1e-12, invalid KR covariance) or KR df \< 1 | `R/peak_gene_hierarchical_helpers.R`, `src/peak_gene_KR.cpp` |
| Hierarchical | P-value | Hardcoded: inline literal | F test with 1 and KR degrees of freedom | `src/peak_gene_KR.cpp` |
| Hierarchical | Multiplicity | Hardcoded: inline literal | BH by class, family size including unreliable tests | `R/peak_gene_hierarchical_helpers.R` |

## Prioritization, top links and plots

SuSiE fine-mapping prioritizes peaks for genes that have at least one conditional link, using the donor- and depth-residualized values. Top-link figures rank estimable positive hierarchical slopes by nominal P-value, excluding self-promoter peaks; gene-body peaks remain eligible, and no significance cutoff is applied. Each figure shows the gene context, the focal cell type's insertion coverage and the donor-residual scatter.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| SuSiE | Input and gene screen | Hardcoded: inline literal | HC3 results; genes with a link at correlation ≥ 0.15, FDR \< 0.05, not self-promoter | `R/peak_gene_finemapping_helpers.R` |
| SuSiE | Genes per branch and peaks per gene | Hardcoded: helper default | 50 genes by best FDR; 500 peaks by absolute correlation; at least 2 variable peaks | `R/peak_gene_finemapping_helpers.R` |
| SuSiE | Model settings | Hardcoded: inline literal | `L = min(10, n_peaks)`, `intercept = FALSE`, `standardize = TRUE`, `estimate_residual_variance = TRUE`, `max_iter = 100`, credible-set coverage 0.95 | `R/peak_gene_finemapping_helpers.R` |
| SuSiE | Records retained | Hardcoded: inline literal | PIP ≥ 0.01 or credible-set members, else the top peak | `R/peak_gene_finemapping_helpers.R` |
| Top links | Links per cell group | Configurable |  | [`peak_gene_correlation_top_links_per_cell_group`](../parameters.html#peak_gene_correlation_top_links_per_cell_group) |
| Top links | Selection | Hardcoded: inline literal | estimable, coefficient \> 0, not self-promoter; ordered by hierarchical P, gene, peak | `R/peak_gene_hierarchical_helpers.R` |
| Plots | Window padding and coverage | Hardcoded: inline literal | 2 kb padding; 500 bins; focal WNN cell type; fragments as read counts | `R/peak_gene_plot_helpers.R` |
| Plots | Coverage clip | Hardcoded: helper default | `clip_quantile = 0.999` | `R/ATAC_tracks_helpers.R` |
| Plots | Histogram bin width and distance bins | Hardcoded: target literal | 0.025; 5 kb bins capped at 245 kb | `module_peak_gene_correlation/correlation_targets.R`, `R/peak_gene_correlation_helpers.R` |
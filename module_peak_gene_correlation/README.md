# Peak–gene correlation

The hierarchical branch scans all pairs passing the existing distance,
detection, aggregate-eligibility and measurement-support filters, using WNN-derived cell types.
There is no HC3 significance, correlation, promoter or top-N screen before
fitting. The former hierarchical-candidate-limit parameter is removed.

The model includes donor fixed intercepts, log-depth covariates and a Gaussian
random peak slope by donor. `src/peak_gene_REML.cpp` profiles residual variance
and independently searches the slope/residual variance ratio, including zero.
Genes sharing a peak reuse an immutable spectral decomposition.
`src/peak_gene_KR.cpp` supplies batched Kenward–Roger inference. Workers compile
both kernels into private temporary caches; file targets track source changes.
The fresh lme4/pbkrtest helper remains a numerical reference.

Results retain slope estimates, donor-slope SD, degrees of freedom, raw p-values
and numerical diagnostics. Inferential p-values are unavailable for insufficient
within-donor peak information, failed fits or effective degrees of freedom below
1. The latter remains a provisional reliability rule. Zero random-slope variance
alone does not invalidate a fit. Finite-search and conditioning failures are
recorded rather than silently accepted.

BH correction is performed across the complete eligible pair family within
each cell type, including promoters and counting unreliable tests in the family
size; unreliable tests retain missing inferential p/FDR values. These Gaussian
models and their approximate tests have not yet undergone broader calibration,
power or influence validation and do not guarantee robustness to outliers.

Top-link figures rank positive, estimable nonpromoter slopes by hierarchical
p-value across the full scan, without a significance cutoff, and report the
hierarchical BH FDR. Existing HC3 results and summary figures remain available
as a distinct analysis. Contextual arcs show HC3-filtered links and HC3 FDR
widths; the focal hierarchical link is included even when absent from that
filter. Missing HC3 FDR uses minimum width. Plot captions distinguish the two
analyses. Hierarchical fitting is cached separately from plot generation.

## Measurement-support filtering

`peak_gene_correlation_filter` is a module setting with allowed values `lenient`
(default in `cfg_pipeline_parameters.tsv`), `moderate`, and `strict`. Override
it for an aggregation in `module_peak_gene_correlation/cfg.yaml`.

| Preset | RNA / ATAC counts at median library depth | Minimum supported aggregates per feature | Shared supported donors | Supported aggregates per feature/donor |
|---|---:|---|---:|---:|
| lenient | 5 / 3 | max(6, ceiling(10%)) | 2 | 2 |
| moderate | 10 / 5 | max(6, ceiling(10%)) | 2 | 2 |
| strict | 10 / 5 | max(10, ceiling(20%)) | 3 | 3 |

Count thresholds scale with each aggregate's library depth relative to the
cell-type median, with a two-count raw floor. Both features must qualify in
the required shared donors, but their qualifying aggregates need not coincide.
The filter removes hypotheses before both HC3 and hierarchical testing; all
observations remain in retained regressions. Excluded hypotheses do not enter
BH correction. Unreliable tests among retained hypotheses still count in the
hierarchical family size. Filtering is based on measurement support, not on
obtaining discoveries; null calibration remains a separate validation task.

`peak_gene_correlation_filter_records.WNN` caches selected candidates and
retention diagnostics for all three presets. `filter_retention_plot` compares
retained fractions by cell type and identifies the active setting. Diagnostic
exclusion counts for RNA, ATAC and donor support overlap and should not be added.
Strict filtering necessarily excludes cell types with fewer than three donors.

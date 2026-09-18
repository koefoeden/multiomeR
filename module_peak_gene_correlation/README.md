# Peak–gene correlation

The hierarchical branch scans all pairs passing the existing distance,
detection and aggregate-eligibility filters, using WNN-derived cell types.
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

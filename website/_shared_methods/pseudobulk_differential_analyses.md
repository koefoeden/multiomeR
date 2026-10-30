For molecular analyses, GEX and ATAC counts of the final WNN nuclei are summed within each donor–annotation-class combination to form pseudobulk samples.
Three feature matrices are tested:

- **Gene expression** (differential gene expression, DGE): gene pseudobulk counts.
- **Chromatin accessibility** (differential chromatin accessibility, DCA): consensus-peak pseudobulk counts.
- **Motif-family accessibility** (differential transcription-factor motif accessibility, DTFA): betterChromVAR bias-corrected deviations of the 233 JASPAR motif families, calculated from the pseudobulk ATAC counts [@germain2026_betterchromvar].
  This matrix is therefore continuous.
  Unlike chromVAR *z*-scores, deviations do not grow with sample depth, and no per-sample centering or cross-sample normalization is applied.

For each model, donors that miss a model variable or fall outside an optional donor list are removed.
Samples outside an optional annotation-class subset are also removed, and the DTFA analysis by default removes samples with fewer than 5,000 ATAC counts.
For count matrices, zero-depth samples are removed, features are filtered with the design-aware `edgeR::filterByExpr()` [@chen2025_edger_v4], and library sizes are TMM-normalized with `edgeR::normLibSizes()` [@robinson2009_edger; @robinson2010_tmm].
TMM assumes that most features do not change in one direction.
A broad one-sided change, such as a global opening of chromatin, would shrink the estimated effects and make unchanged features appear to change the other way.

Model matrices are built from a configured formula or custom design function, and contrasts from named linear expressions or custom contrast functions.
A configured random effect selects `limma::duplicateCorrelation()` with the random effect as block [@smyth2005_duplicatecorrelation; @ritchie2015_limma].
For counts, voom precision weights [@law2014_voom] are re-estimated for the block and correlation (`edgeR::voomLmFit()`).
Weighted fits with a block use a closed-form solution of limma's generalized least squares, computed for all features at once, which agrees with `limma::lmFit()` to floating-point error.
A contrast of several coefficients takes each feature's exact standard error from its own weighted design, rather than `limma::contrasts.fit()`'s approximation, which uses one coefficient correlation for all features.
A model with several samples of one donor in an annotation class must block on that donor, as its samples are not independent replicates.
Wherever an intra-block correlation is estimated, the consensus is taken from at most 50,000 evenly spaced features, enough to estimate this one value pooled across all features [@smyth2005_duplicatecorrelation] at a fraction of the runtime.
The models are still fitted to all features.
Without a random effect, count designs with more than one coefficient are fitted with edgeR quasi-likelihood models and tested with quasi-likelihood F-tests [@lund2012_quasi_likelihood; @chen2025_edger_v4].
These models use robust empirical Bayes moderation of the dispersions [@phipson2016_robust_eb].
Continuous matrices use limma.
edgeR cannot test a single-coefficient design, so such a count design is fitted with voom–limma after the same filtering and normalization.
An optional paired-cell-type route fits each annotation class separately with voom–limma, or limma for continuous matrices.
It allows at most one pseudobulk per pairing unit and class.
A cross-class contrast is tested on a combined response, restricted to the donors present and the features tested in every class involved.
This response is the contrast-weighted sum of the classes' log-expression.
It is weighted by the inverse of its variance, which is computed from the classes' precision weights.
The combined response is fitted with the same donor-level design and correlation block, so each feature has its own covariance between classes.

For each model and contrast, Benjamini–Hochberg FDR is calculated across the tested features, and features with an FDR below 0.05 are significant.
DGE contrasts are also tested against the Hallmark and Reactome collections of MSigDB release 2026.1 [@liberzon2015_hallmark; @milacic2024_reactome] with competitive `camera` tests [@wu2012_camera; @dolgalev2026_msigdbr].
These tests use the contrast's expression, precision weights and design, or the combined response for cross-class contrasts.
Each set's inter-gene correlation is estimated from the residuals, and each set needs at least ten tested genes.
camera cannot account for a random effect or correlation block, so such models instead test the blocked fit's moderated *t*-statistics with `limma::cameraPR()`, passing camera's estimated correlations, floored at zero.
The 13 mitochondrially encoded genes are left out of the gene-set tests, but kept in the per-gene results: in nuclei, their reads mostly come from ambient and cytoplasmic RNA.
FDR is calculated within each contrast and collection.
A set's direction compares its genes with all other tested genes, including genes that belong to no set.
The sets of one contrast therefore often lean the same way, even under the null.

Transcription-factor regulons of the CollecTRI network [@muller_dott2023_collectri] are tested as signed gene sets in the same competitive way.
The network is used only when its checksum matches the recorded one, and complexes such as AP1 and NF-κB are kept as single regulons.
For every contrast, `limma::cameraPR()` compares the moderated *t*-statistics of a regulon's targets, with those of repressed targets negated, against all other tested genes.
Models fitted with edgeR quasi-likelihood tests, which report no *t*-statistic, use the signed normal quantile of each gene's *P* value instead.
Each regulon's inter-gene correlation is camera's estimate for the signed set, floored at zero.
Each class fit sets the gene universe, and a regulon needs at least ten tested targets.
A significant regulon shows that its targets moved together in the expected directions, not that the transcription factor itself changed.
The results therefore name each regulon's leading targets and the factor's expression, pooled over the pseudobulks of each class the contrast involves, and flag factors below 0.5 counts per million.
The network is human, so mouse aggregations stop with an error.
A cross-modality comparison maps CollecTRI regulators to JASPAR families.
It compares each regulon's signed enrichment, the normal quantile of its *P* value, with the motif-family accessibility *t*-statistics, and summarizes families by the median over their regulators.
Each contrast's Spearman correlation is shown against ±1.96/√(*n* − 1), the 95% range for unrelated rankings of its *n* families.
Where gene expression is fitted with limma, the comparison also reports TF expression (edgeR tests report no *t*-statistic).
An analysis stops with an error rather than returning a model when no samples or testable features remain or residual degrees of freedom are below one.
It also fails when a correlation block or contrast is invalid, or when too few shared donors remain for a requested paired-cell-type comparison.

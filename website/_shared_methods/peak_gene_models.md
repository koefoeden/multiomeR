GEX and ATAC values are residualized against donor, scaled log RNA and ATAC depth and the class's accessibility components described below.
The Pearson correlation of the residuals is reported as a descriptive summary.
The hierarchical analysis fits donor fixed intercepts, depth covariates, the accessibility components and a donor-varying slope for the within-donor-centered peak value.
It uses profiled restricted maximum likelihood, with Kenward–Roger inference for the average slope [@kenward1997_kenward_roger].
When the slope varies between donors, donors rather than pseudobulks are the units of inference.
Kenward–Roger adjusts its standard error and degrees of freedom for the small number of donors.
When the estimated slope variance is small, its degrees of freedom move towards the number of pseudobulks.
A peak is tested only when its accessibility varies within at least eight donors.
With fewer, the random-slope test has too few degrees of freedom to be informative, so its pairs remain untested.
Fits that fail numerical checks keep their estimates but get no P-value.
Peaks and genes that both vary along shared chromatin states correlate whether or not they interact.
The accessibility components absorb the leading states, but not necessarily a gradient within the class that shapes only part of its accessibility, such as a mix of subtypes.
They are the top ten principal components, or one per ten residual degrees of freedom when fewer.
They are computed from all candidate peaks, across chromosomes, that are accessible in at least 5% of the class's pseudobulks.
Each peak enters as its donor- and depth-residualized, standardized log1p CPM.
Genes that vary along the remaining states still correlate with many peaks.
Each candidate is therefore compared with its own gene's background, like Signac's comparison of each peak–gene correlation with matched background peaks [@stuart2021_signac].
This background consists of up to 200 candidate peaks that lie on the same chromosome but more than 10 Mb from the gene's TSS, scored with the same model.
The candidate's signed normal score is the normal quantile of its P-value, signed by the slope.
This score is standardized by the mean and SD of the gene's background scores, given at least ten of them.
It is then rescaled by the class's background SD within its peak-detection quintile.
A candidate's empirical P-value is the share of the class's rescaled background scores that are at least as extreme.
These P-values are Benjamini–Hochberg-corrected within each annotation class over the complete eligible pair family, counting untested and unreliable pairs.
A link requires an estimable positive average slope above the gene's background mean and an FDR below 0.05, excluding self-promoter peaks.
Links whose peak lies in the target gene's own body are reported separately as intragenic associations, and the others as candidate enhancer links.
Much of a gene's body can vary with its expression, so the strongest such peak need not mark an enhancer.
These model-based associations do not establish causal enhancer–gene regulation.
Background peaks more than 10 Mb away also cannot separate a specific enhancer from covariation within a local chromatin domain.
An empirical P-value cannot fall below one over the number of background scores plus one, so in classes with many tested pairs a few strong pairs cannot pass the first Benjamini–Hochberg step alone, and link counts can change in jumps.

Top-link plots rank estimable positive hierarchical associations by their uncalibrated Kenward–Roger P-value, without a significance cutoff.
They exclude self-promoter peaks but keep gene-body peaks.
Each plot shows the gene context, the plotted class's insertion coverage and the donor-residual scatter.

Each candidate pair also gets an activity-by-contact (ABC) score [@fulco2019_activity_by_contact].
It needs no variation between donors, so it also covers cell types that the hierarchical model cannot test.
Scores are computed for each labelled WNN cell type with at least 100 nuclei, the genetic-enrichment floor for a cell type's summed profile.
A gene is scored in a cell type when its summed GEX counts there reach 10 counts and 1 count per million, the rule of the S-LDSC gene sets.
A peak's activity is its ATAC counts per million in the cell type, as in the ABC pipeline without H3K27ac data.
Its contact with a gene falls with the distance from the peak centre to the TSS as a power law with exponent 1.024, the ABC pipeline's default [@nasser2021_enhancer_maps].
Peaks closer than 5 kb get the contact at 5 kb.
The score is activity times contact, divided by the sum of this product over all consensus peaks within 5 Mb of the TSS, including the gene's own promoter.
The score is reported without a threshold.
The published threshold of 0.02 was calibrated with DNase-seq, H3K27ac and Hi-C data in K562 cells, so it does not carry over to this data.

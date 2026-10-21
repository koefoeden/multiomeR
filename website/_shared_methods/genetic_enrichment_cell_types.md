Annotation-class pseudobulk deviations are calculated separately rather than by averaging the nucleus-level results.
ATAC counts are summed by the final WNN cell-type labels, which the other analyses also use.
Classes with fewer than 100 nuclei, such as isolated unassigned nuclei that a small minimum cluster size keeps, are left out of these and the S-LDSC analyses, because their profiles are too sparse and every class weighs equally in the background model.
A betterChromVAR background model is fitted to the resulting peak-by-class matrix.
The expected accessibility of each peak averages the depth-normalized class profiles with equal weight.
Each class is thus compared with the average class rather than with a reference dominated by the most abundant classes.
Background peaks are matched on this expectation and on GC content.
The deviation is the accessibility of the weighted peaks minus that of the matched background peaks, relative to their expected accessibility.
It thus corresponds to chromVAR's bias-corrected deviation.
Each trait's deviations average about zero across classes.
A class's mean deviation over many traits therefore reflects the panel of traits rather than a property of the class.
A relative deviation also standardizes the deviations across annotation classes within each trait.
Significance is tested against a null that keeps the clustering of real GWAS loci, which often place several correlated peaks at one locus.
Each trait's credible sets are moved along their chromosomes 10,000 times and scored like the real ones.
Credible sets whose lead variants lie within 500 kb of each other move together by 0.5 to 2 Mb, keeping their variant spacing and posterior probabilities.
No moved credible set may land on an original one, in the extended MHC (chr6:25–34 Mb), in a peak-free stretch of at least 1 Mb such as a centromere, or off its chromosome.
A block without such a move within 2 Mb moves elsewhere on its chromosome.
Moved loci fall in fewer and less accessible peaks than real ones, so their raw deviations vary more.
The real and moved loci are therefore compared by their background z-scores.
This is the deviation's numerator divided by its standard deviation under independently drawn matched background peaks, which scales with the squared peak weights.
betterChromVAR 1.1.10 overestimates this background variance for deep counts such as pseudobulks, so it is recomputed from the counts.
A class's one-sided P-value is the share of moved copies, plus the real trait, whose background z-score is at least as large.
The reported z-score standardizes the real background z-score by the mean and standard deviation of the moved copies.
This null is skewed to the right, so normal P-values would be too small.
Benjamini–Hochberg-adjusted P-values across all traits and annotation classes are reported.
Heatmaps and locus-contribution plots mark classes with one star at P < 0.05 and with two at q < 0.05.
When effect sizes are available, an absolute-effect version weights variants by posterior probability times absolute effect size.
When only loci carry effects, each locus's lead-variant effect is used.
These weights are rescaled to sum to the kept variants' posterior probabilities.
The absolute-effect version is tested against the same moved copies, weighted in the same way.
Each class deviation is split into additive contributions from its loci and variants, which sum to the class totals.
The loci and variants are labeled with the highest-scoring Open Targets locus-to-gene genes [@mountjoy2021_l2g].
Detail plots are drawn for classes with a positive deviation and a z-score of at least 1 by default, in either weighting.
They show the top loci by absolute contribution and by combined contribution and effect-size rank.

Each trait with a class at q < 0.05 gets a top set of classes.
Its anchor is the class with the largest deviation among those at q < 0.05.
The top set holds the anchor and every class not significantly below it, so classes with larger but nonsignificant deviations join it.
Credible sets whose lead variants in peaks lie within 500 kb of each other form one locus.
In 1,000 bootstrap samples of these loci, the sampled loci's contributions are summed for each class.
Each class keeps the background expectation and denominator of all loci, so the samples approximate rather than repeat the analysis.
They give the standard error of each class's shortfall from the anchor.
A class is excluded when a one-sided normal test of its shortfall reaches P < 0.05 divided by the number of other classes.
Whichever class is the anchor, the Bonferroni inequality ensures that the top set contains the class with the largest true deviation with at least 95% probability.
Traits with fewer than 10 loci get no top set, because so few loci give no reliable bootstrap.
Heatmaps outline the top sets.

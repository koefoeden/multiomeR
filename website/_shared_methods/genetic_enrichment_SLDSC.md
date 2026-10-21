Stratified LD score regression (S-LDSC) tests whether a GWAS's heritability is enriched in variants within the peaks accessible in each annotation class, using genome-wide summary statistics rather than fine-mapped loci [@finucane2015_stratified_ldsc; @finucane2018_ldsc_seg].
It uses the peaks and the ATAC counts summed by the final WNN cell-type labels of the pseudobulk chromVAR analysis.
Each class's counts are divided by its total, and each peak's depth-normalized values are divided by their sum across classes.
A peak's annotation value for a class is thus its share of the summed class profiles: the chromVAR background model's observed-over-expected ratio divided by the number of classes.
A peak equally accessible in all classes therefore gets the same small value in each.
An all-peaks annotation marks every peak with reads.
Reference variants inside a peak take its values, and all other variants take zero.
LD scores of these annotations are computed from the 1000 Genomes European samples on GRCh38 [@auton2015_1000_genomes] for HapMap3 regression variants outside the MHC, with ldsc's bias-corrected squared correlation summed over a 1 cM window [@bulik-sullivan2015_ldsc].
Summary statistics are taken from the GWAS Catalog, from its harmonized file where one exists and otherwise from the submitted file [@cerezo2025_gwas_catalog].
A GWAS can instead name another Catalog study or its own summary statistics in GWAS-SSF format on GRCh38 with rsIDs [@hayhurst2022_gwas_ssf], such as a consortium release that the Catalog lacks.
The statistics are cleaned largely as ldsc does, without its imputation-quality filter.
Variants need an rsID, a P-value within the double-precision range, alleles that match the reference without strand ambiguity, and a minor-allele frequency above 0.01 when one is reported.
Variants whose sample size is below two-thirds of the 90th percentile are removed.
Studies without per-variant sample sizes use the total of the samples in the summary statistics' GWAS-SSF metadata, which also gives the sample size and ancestry shown for each trait.
Each class is fitted separately, jointly with the all-peaks annotation and the 97 annotations of the baseline-LD v2.2 model [@gazal2017_baseline_ld].
The class coefficient therefore measures heritability per variant beyond accessible chromatin in general and the baseline functional categories.
Regression weights follow ldsc but are computed once per GWAS from the shared annotations.
Variants with a chi-square statistic above 80, or above 0.001 times the largest sample size when that is larger, are excluded.
Standard errors come from a jackknife over 200 blocks of adjacent variants.
A one-sided P-value tests for a positive coefficient, and Benjamini–Hochberg-adjusted values across all traits and classes are reported.
The standardized effect size tau* scales the coefficient by the annotation's standard deviation and the number of common variants relative to total heritability.
It is left missing when the total heritability estimate is not positive.
GWAS without genome-wide summary statistics in the GWAS Catalog or a named file are not analyzed, and non-European GWAS are analyzed with the European reference.

Cell-type gene expression is tested in the same way, following LDSC-SEG [@finucane2018_ldsc_seg].
The GEX counts of the donor pseudobulks are summed by class and divided by each class's total.
A gene's specificity for a class is its share of these depth-normalized values summed across classes [@skene2018_brain_cell_types].
A gene counts as expressed in a class with at least 10 counts and one count per million there.
Each class's gene annotation marks the reference variants within 100 kb of the gene bodies of its most specific expressed genes.
These number 10% of the autosomal protein-coding genes expressed in any class, and an all-genes annotation marks the windows of all of them.
Where windows overlap, a variant takes the largest value, which for these binary annotations is their union.
Each class is fitted in three models with the baseline-LD annotations: its peak annotation with the all-peaks annotation, its gene annotation with the all-genes annotation, and both together with both controls.
In the joint model, each coefficient measures the heritability that one modality explains beyond the other.
All models share the regression weights of the baseline and all-peaks annotations and the same jackknife blocks.
Benjamini–Hochberg adjustment is applied separately to each model and annotation.

For each model and tested annotation, each trait with a class at q < 0.05 gets a top set of classes.
Its anchor is the class with the largest tau* among those at q < 0.05.
The top set holds the anchor and every class not significantly below it, so classes with larger but nonsignificant tau* join it.
All classes share the jackknife blocks.
The standard error of each class's shortfall from the anchor therefore comes from the paired delete-one-block estimates of both tau* values, each scaled by its full-data total heritability.
A class is excluded when a one-sided normal test of its shortfall reaches P < 0.05 divided by the number of other classes.
Whichever class is the anchor, the Bonferroni inequality ensures that the top set contains the class with the largest true tau* with at least 95% probability.
Heatmaps of the peak model colour each class by tau*.
They mark classes with one star at P < 0.05 and with two at q < 0.05, and outline the top sets.

The robustness of each coefficient at q < 0.05 is read from the jackknife.
The blocks whose single deletion lowers the coefficient most are deleted together, 1, 5 or 10 at a time, and the model is refitted.
The share of the coefficient they carry is reported, with the z-score without them from a jackknife over the remaining blocks.
A coefficient that spread evenly over the genome would lose a share close to the deleted fraction of blocks, plus noise.
Blocks span about 15 Mb, so this measures concentration only coarsely.
Variants above the chi-square cap never enter the fits, so the loci with the largest effects never contribute.

Annotation strategies are compared at matched sizes, because enrichment depends on how much of the genome an annotation covers.
For each class, peaks are ranked by their accessibility share and gene windows by expression specificity.
A variant in overlapping windows takes the highest specificity among them.
Peaks or window pieces are added in this order until they hold 0.5%, 1%, 2% or 5% of the common (MAF > 5%) reference variants.
Pieces are kept whole and ties keep genomic order, so an annotation can exceed its share by one piece.
All consensus peaks together hold about 10% of common variants, so larger peak annotations would lose their class specificity.
A reference annotation of each strategy ranks the same peaks or genes at random.
Each binary annotation is fitted with the 97 baseline-LD annotations, which also set the regression weights.
Its heritability share sums the fitted per-variant heritability of the whole model over the annotation's common variants, as ldsc's Prop._h2 does [@finucane2015_stratified_ldsc].
Heritability enrichment divides this share by the annotation's share of common variants, with a ratio jackknife standard error over the same 200 blocks.
The chi-square cap removes the strongest loci from this statistic.
A complementary statistic uses the fine-mapped credible sets of each GWAS instead.
Fine-mapped enrichment divides the annotation's share of the summed posterior inclusion probabilities of autosomal credible-set variants by its share of common variants.
Variants are assigned by position, so indels and variants missing from the reference count.
Its uncertainty comes from 1,000 bootstrap samples of loci, which group the credible sets whose centres lie within 500 kb of each other.
Fine-mapping methods differ between sources and some use functional priors, so each result keeps its source's method.
Pairs of GWAS and class are called enriched when a single-modality model reaches q < 0.05.
Summaries across pairs or GWAS report medians with 95% intervals from 1,000 bootstrap samples of GWAS.
A strategy's specificity is each GWAS's median enrichment in its enriched classes divided by that in its other classes, when both are positive.
Small peak annotations hold only the credible-set variants inside them, while gene windows hold whole credible sets, so imprecise fine-mapping lowers the fine-mapped enrichment of peaks more.

Agreement with pseudobulk chromVAR is summarized for each GWAS by the Spearman correlation of the two methods' z-scores across classes.
For each class, it is the correlation across GWAS after each method's z-scores are standardized within each GWAS.

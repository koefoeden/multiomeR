Each multiomeR version uses one fixed Open Targets Platform release, currently 26.03, and downloads its study, credible-set, credible-set evidence and target datasets [@buniello2025_open_targets_platform].
Users can supply their own fine-mapping results as Parquet files in a fixed format.
For each configured GWAS, multiomeR retrieves the requested study and fine-mapping result.
When selection is automatic, it takes the first available of SuSiE, SuSiE-inf and PICS fine-mapping [@wang2020_susie; @cui2024_susie_inf; @farh2015_pics].
Open Targets credible sets flagged "Variant not found in LD reference" are left out, because their lead variant was given all the posterior probability without fine-mapping.
Open Targets' other flags on these studies describe provenance, such as fine-mapping with out-of-sample LD or from curated top hits, and are kept.
Variants of the study's credible sets (95% sets for Open Targets) with a posterior probability above 0.001 by default are mapped to consensus ATAC peaks.
A credible set's posterior probabilities sum to its expected number of causal signals, which exceeds one in multi-signal regions such as those of FINEMAP.
The rounded sum counts such a set's signals in the loci tracks, and labels of its lead variant name them.
Variant weights that map to the same peak are summed and capped at one.
These peak weights form the trait annotation from which analytic nucleus-level chromVAR z-scores are calculated, reusing the GC-bias background of the motif analysis [@schep2017_chromvar; @ulirsch2019_gchromvar; @germain2026_betterchromvar].
The z-scores use betterChromVAR's finite-population correction, which accounts for each nucleus's expectation being estimated from its own counts and sums the background variance over the squared peak weights.
To show which traits are not independent, credible sets are also compared pairwise.
Posterior probabilities are summed per variant and normalized to one within each GWAS.
The similarity of two GWAS is the sum, over identical variants, of the smaller of their two values.
It equals one minus the total variation distance and is not a colocalization test.

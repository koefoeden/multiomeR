Transcription-factor motifs are represented by the 233 sequence-similarity families of the JASPAR 2026 CORE vertebrate collection [@ovekbaydar2026_jaspar].
Families are used because related transcription factors often share near-identical motifs, whose matches and accessibility cannot be told apart.
The pipeline includes a copy of the familial root motifs and the family membership table.
Each root motif is trimmed of flank columns whose total count is less than half of its largest column total.
These columns come from a minority of the family's aligned member motifs, so their counts are sparse and can look specific by chance.
Scanning with them would restrict matches to that minority's sub-family.
A family that merges different motif types is thus scored mainly by its most common type, and its other members are captured poorly.
For example, the root motif of the family that holds both forkhead and MEF2 factors matches forkhead sites.
The trimmed root motifs are matched to the consensus peaks with `motifmatchr` at its default threshold of P < 5 × 10^−5^ [@schep2025_motifmatchr].
A motif needs at least eight columns to reach this threshold, so the trim keeps the best-supported eight columns around a shorter core.
Root motifs that are shorter to begin with may match only their best-matching sequence.
That sequence can occur by chance up to six times as often as the threshold allows.
betterChromVAR applies the chromVAR workflow with GC-bias correction to compute analytic motif-family deviations and z-scores for individual nuclei, with its finite-population correction of the z-scores [@schep2017_chromvar; @germain2026_betterchromvar].
All families are tested per ATAC annotation class with one-versus-rest Wilcoxon tests, with Benjamini–Hochberg correction [@benjamini1995_fdr] across all family–class tests.
They are also summarized as cell-weighted means per ATAC cluster, GEX cluster and GEX cell type.
Configured marker transcription factors are scored in the same way by their own JASPAR motifs, and these scores are shown on embeddings.
A transcription factor's name selects all of its motifs, and a name with a motif identifier selects one.
These values measure accessibility associated with a motif or motif family and are not direct measurements of transcription-factor activity.

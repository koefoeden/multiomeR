Locus attribution is an experimental route, off by default, that asks which cell types each credible set points to through its own accessibility.
It uses the peak weights of the cell-type enrichment: each pair of a credible set and a peak carries the posterior probability (PIP) of the set's credible variants in the peak.
A cell type takes a peak in proportion to the peak's share of its ATAC reads, summed by the final WNN cell-type labels.
Each cell type's profile is shrunk toward the pooled profile of all cell types by the reads of 50 average nuclei.
A peak without reads in a small cell type therefore counts as average rather than closed.
A uniform background takes peaks that no cell type explains better than reads spread evenly over all peaks.
Without it, the cell types with the most diffuse reads, such as neurons in a peripheral tissue, would collect every weak peak.
Each pair's PIP is split between the cell types and the background in proportion to these likelihoods, and summed over each credible set's peaks.
A cell type's fold is its share among the cell types divided by an even share.
A fold of 1 thus means that the locus's accessibility does not favour the cell type.
A compartment's fold is the mean of its cell types' folds, so that large compartments do not lead by their size alone.
These shares and folds use no information from other loci.
They are attributions, not tests.
A random peak of the same GC content and accessibility is usually about as specific to some cell type, so a single locus's attribution is weak evidence on its own.

For context, a mixture over the cell types and the background is fitted to all of a trait's pairs, each weighted by its PIP, by expectation–maximization [@dempster1977_em].
Its weights estimate which cell types the trait's loci act through as a whole.
A Dirichlet prior of 0.1, 1, 10 or 100 loci, spread evenly over the cell types, is chosen by the likelihood of the pairs on held-out chromosomes in five folds.
The weights are reported only when they predict the held-out pairs at least one nat better than even weights, which split each pair by its accessibility alone.
This is positive evidence on the scale of Kass and Raftery [@kass1995_bayes_factors].
Traits with fewer loci keep their locus attributions, and their weights are marked as not estimable.
Each cell type's share by chance is the mean of the weights fitted to ten sets of random peaks.
Each random peak comes from the same bin of GC content and expected accessibility as the real peak it replaces, as in the background model of the cell-type enrichment.
The random peaks are drawn with a fixed seed.
Loci are labelled with the highest-scoring Open Targets locus-to-gene genes [@mountjoy2021_l2g].

The differential-analysis module treats donors, rather than nuclei, as the biological replicates.
It takes annotation classes from the WNN cell-type label of the final WNN metadata.
For cell-type composition, nuclei are counted by donor and annotation class, adding missing combinations as zero counts.
Unless the model lists the response classes, the assigned classes observed in at least half of the donors are tested.
Rarer classes and unassigned clusters are not tested, as they would only add to the multiple-testing burden, but they still count in every donor's total.
Each class is modeled separately in a fixed-effects beta-binomial model with logit link, fitted by `glmmTMB` with the configured formula [@brooks2017_glmmtmb].
Its response has two columns: the numbers of nuclei in that class and in all other classes.
Named linear contrasts of the fixed effects are tested with Wald statistics.
Fits that fail, do not converge or lack a positive-definite Hessian are reported as non-estimable.
Benjamini–Hochberg correction is applied across the tested annotation classes within each model and contrast, but not across contrasts or models.

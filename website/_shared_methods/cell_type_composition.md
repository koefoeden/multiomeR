The differential-analysis module treats donors, rather than nuclei, as the biological replicates.
It takes annotation classes from the WNN cell-type label of the final WNN metadata.
For cell-type composition, nuclei are counted by donor and annotation class, adding missing combinations as zero counts.
Unless the model lists the response classes, the assigned classes observed in at least half of the donors are tested.
Rarer classes and unassigned clusters are not tested, as they would only add to the multiple-testing burden, but they still count in every donor's total.
Each class is modeled separately in a fixed-effects beta-binomial model with logit link, fitted by `glmmTMB` with the configured formula [@brooks2017_glmmtmb].
Its response has two columns: the numbers of nuclei in that class and in all other classes.
Named linear contrasts of the fixed effects are tested with Wald statistics against a *t* distribution with the fit's residual degrees of freedom.
The beta-binomial dispersion is estimated from the same donors, so a normal reference would be liberal with few donors.
Each level of a tested factor needs at least three donors; with fewer, no reference distribution is reliable.
Fits that fail, do not converge or lack a positive-definite Hessian are reported as non-estimable.
Benjamini–Hochberg correction is applied across the tested annotation classes within each model and contrast, but not across contrasts or models.

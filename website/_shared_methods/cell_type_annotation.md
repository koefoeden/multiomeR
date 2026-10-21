Cell-type annotation uses GEX marker signatures supplied by the user.
Signatures list positive markers and, optionally, genes expected to be absent.
Genes missing from the count matrix stop the run with an error before scoring.
Scores are computed with a BPCells-based implementation of UCell scoring [@andreatta2021_ucell; @parks2025_bpcells].
It processes the raw GEX counts in blocks of limited size and reproduces UCell's per-cell descending ranks, rank truncation at 1,500 genes, positive and negative signatures and lower-bound clipping.
Working in blocks avoids loading the whole matrix or building a Seurat object only for scoring.
The pipeline does not use the UCell package itself.
Validation instead requires agreement with the UCell reference to within 10^−12^.
Signed signatures are clipped at zero per cell before averaging.
The same scoring annotates the GEX, ATAC and WNN clusters.
WNN labels can use their own marker sets, whose matched controls come from the same GEX reference nuclei.

Raw UCell scores are not comparable between signatures, because some gene sets score highly by chance.
Each label's cluster-level mean score is therefore compared with 999 random control signatures matched on gene abundance and detection.
The control reference samples up to 50 nuclei per GEM well from the GEX metadata before doublet filtering.
The GEX and ATAC annotations thus share the same controls, and the WNN annotation draws its own from the same reference, whatever nuclei each stage keeps.
For each control replicate, all marker genes of the panel are taken in random order.
Each is replaced by a gene drawn from its 50 nearest eligible candidates not yet used in that replicate.
Candidates exclude all marker genes and undetected genes.
A label's adjusted score in a cluster is its observed mean minus the 95th percentile of its matched controls.
The label with the highest adjusted score is the candidate.
Its advantage is the smaller of its lead over zero and its lead over the runner-up.
The candidate is assigned when its advantage is positive and at least a minimum configured separately for the GEX, ATAC and WNN clusters (0.05 by default).
Otherwise, the cluster remains unassigned, with the candidate and reason kept.
Requiring a lead over both the background and the runner-up leaves ambiguous clusters unassigned rather than forcing a label.
Raising the minimum can therefore only make assigned clusters unassigned, never relabel them.
Each unassigned cluster keeps its own label, `Unassigned_cluster_<id>`, so doublet detection, pseudobulks and per-cell-type tests never pool unrelated unassigned clusters.
Because annotation depends on the supplied signatures and their level of detail, labels can represent either cell types or broader groups of cells, such as a lineage or a tissue of origin.

Annotation outputs also report, for each cluster, how many positive markers are detected in at least 10% of its nuclei.

Each GEM well supplies a Cell Ranger ARC count directory.
Gene-expression counts are imported into BPCells from the filtered feature-barcode matrix, or from a CellBender output when configured [@parks2025_bpcells; @fleming2023_cellbender].
ATAC fragments are imported from the Cell Ranger fragment file.
Both carry a GEM-well-prefixed barcode, so nuclei remain distinct across wells.
Per-well quality-control metrics are calculated for the barcodes that Cell Ranger called as cells.
BPCells calculates ATAC metrics such as TSS enrichment and nucleosome signal.
The blacklist fraction is the share of a nucleus's fragments with at least one cut site in a blacklisted region of the genome [@amemiya2019_encode_blacklist; @ogata2023_excluderanges].
Gene-expression metrics, such as the fraction of counts from mitochondrial genes (names beginning with `MT-`, ignoring case), are calculated from the imported matrix.
Optional genotype demultiplexing runs cellsnp-lite on the ATAC alignments, with a minimum minor-allele frequency of 0.1 and 20 reads per SNP [@huang2021_cellsnp_lite].
Vireo then demultiplexes the barcodes with the supplied donor genotypes [@huang2019_vireo].
The most likely singlet donor replaces the GEM-well donor identifier of every barcode.
Vireo's singlet, doublet and unassigned calls and probabilities are recorded as metadata for the exclusions.
A BPCells-based implementation of AMULET reads the stored fragments of up to 1,000 bp one by one, without loading them all into memory [@thibodeau2021_amulet].
This avoids the chromosome-scale `GRanges` objects that scDblFinder builds from the text fragment file [@germain2022_scdblfinder].
It skips fragments on chrM, chrX and chrY and, for GRCh38, fragments in the AMULET authors' list of repeats, segmental duplications and blacklisted regions, where multi-mapping reads cause spurious overlaps.
It reproduces the overlap metrics and q-values of `scDblFinder::amulet()` but applies no automatic filter.
Because AMULET performs well only on deeply sequenced nuclei, filtering on its q-values is left to the user, and its authors recommend q < 0.01 [@thibodeau2021_amulet].

Configured GEM-well exclusions are filter expressions applied to the called cells of each well.
Demultiplexing and AMULET results can therefore be used as exclusion criteria before aggregation-level analysis.
Called cells missing from the imported gene-expression matrix are always excluded.
The remaining nuclei of the selected wells are combined without further filtering.

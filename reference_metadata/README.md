# Cell Ranger ARC reference metadata

These small `reference.json` files are copied from the corresponding official
10x Genomics Cell Ranger ARC reference archives. multiomeR uses their genome,
organism, FASTA/GTF hashes, and Gencode release; it does not need the reference
FASTA, indexes, or annotation files themselves.

| Reference | Upstream archive |
|---|---|
| GRCh38 2020-A | <https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz> |
| GRCh38 2024-A | <https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2024-A.tar.gz> |
| mm10 2020-A | <https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz> |

multiomeR matches the FASTA/GTF hashes in each `atac_fragments.tsv.gz` header
against these JSON files. For another standard or custom reference, copy its
`reference.json` into a subdirectory here. Exactly one JSON must match; no
per-GEM-well reference configuration is needed.

# JASPAR2026 vertebrate familial root motifs

`JASPAR2026_vertebrate_familial_root_motifs.tf` contains the 233 root position
frequency matrices from the official JASPAR2026 CORE vertebrate clustering.
The file is copied without modification from
`vertebrates_CORE_motifs/root_motifs/Root_motifs.tf` in the official clustering
archive:

```text
https://jaspar2026.elixir.no/static/clustering/2026/vertebrates/CORE/interactive_trees/JASPAR_2026_matrix_clustering_vertebrates_archive.zip
```

Its SHA-256 digest is
`b0bb83653ba21a8afc693cf4d941d2a3b699f678ec51bfb2b1f16f0083b29ca0`.

`JASPAR2026_vertebrate_motif_families.tsv` maps the 1,019 non-redundant CORE
vertebrate motifs to those 233 clusters. The rows were collected on 2026-08-10
from the official cluster summaries at:

```text
https://jaspar.elixir.no/static/archetypes/2026/vertebrates/clusters_info_tabs/cluster{1..233}.txt
```

The pipeline scans only the 233 root motifs. Individual family members are
retained as metadata for resolving configured TF markers and mapping external
TF resources such as CollecTRI to the corresponding JASPAR family.

The source is licensed under CC BY 4.0. JASPAR2026 is described in Baydar Ovek
et al., *Nucleic Acids Research* 2026, https://doi.org/10.1093/nar/gkaf1209.

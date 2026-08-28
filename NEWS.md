# multiomeR 0.4.1.9000

- Made fresh clones directly runnable after downloading the public demo data:
  configuration and local controller files are now regular files, and only the
  two `immune_human_2x` GEM wells and that aggregation are active by default.

- Added a restart-safe `setup-demo` Pixi task that installs the environment,
  downloads the two configured public inputs, and installs GitHub-only R
  dependencies behind one quickstart command.

- Replaced the public reaction-based API and configuration vocabulary with 10x
  Genomics GEM well terminology. Existing configurations must migrate to
  `cfg_GEM_wells.tsv`, `GEM_well_ID`, and the corresponding `GEM_well_*`
  parameter names documented in the migration table.

- Made manuscript benchmark wall-time results portable across targets stores
  and removed machine-specific figure paths ([#4](https://github.com/koefoeden/multiomeR/pull/4)).

- Replaced Seurat-backed cell-cycle module scoring with BPCells-native control-binned scoring, restored BPCells-native UCell marker scoring for signed marker sets, and added synthetic parity validation against Seurat and UCell reference implementations.

Initial public beta snapshot.

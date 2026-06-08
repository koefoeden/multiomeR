[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

# multiomeR

multiomeR is a targets-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome data. It is designed as a lean, readable framework that users can adapt to their own studies rather than as a black-box command-line pipeline.

The active workflow is a single root `targets` project driven by `_targets.R`, `cfg_reactions.tsv`, `cfg_datasets.yaml`, and `cfg_aggregations.yaml`.

## Status

multiomeR is in beta. The workflow is actively changing, and public interfaces may still change.

## Documentation

The user manual is built from the Quarto book in `website/`:

- Quickstart: <https://koefoeden.github.io/multiomeR/quickstart.html>
- Selected outputs gallery: <https://koefoeden.github.io/multiomeR/output_gallery.html>
- Full documentation site: <https://koefoeden.github.io/multiomeR/>

## Requirements

- Linux x86_64 system with at least 64 GBs of RAM with plenty of disk-space.
- 10x Genomics Multiome outputs from `cellranger-arc count`.
- Cell Ranger ARC reference directories matching the genomes in your config.

Most R packages and command-line tools are managed through Pixi via `pixi.toml` and `pixi.lock`.

## Contributions

Bug reports and broadly useful feature requests are welcome, especially when they affect users analyzing 10x Multiome data. The project prioritizes lean, inspectable workflow changes over broad abstractions or site-specific convenience layers. See `.github/CONTRIBUTING.md`.

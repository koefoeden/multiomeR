<p align="center">
  <img src="website/figures/multiomeR-logo.svg" alt="multiomeR logo" width="900">
</p>

[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Docs](https://github.com/koefoeden/multiomeR/actions/workflows/docs.yaml/badge.svg)](https://koefoeden.github.io/multiomeR/)
[![R](https://img.shields.io/badge/R-%3E%3D%204.5-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Pixi](https://img.shields.io/badge/env-pixi-f3c638)](https://pixi.sh/)
[![targets](https://img.shields.io/badge/workflow-targets-276DC3)](https://docs.ropensci.org/targets/)
[![BPCells](https://img.shields.io/badge/BPCells-native-2A9D8F)](https://bnprks.github.io/BPCells/)

# multiomeR

multiomeR is a targets-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome data. It is designed as a lean, readable framework that users can adapt to their own studies rather than as a black-box command-line pipeline.

The active workflow is a single root `targets` project driven by `_targets.R`, the settings in [`configuration/`](configuration/README.md). To select a separate settings directory for your project, use the ignored root `configuration.local` file.

## Status

multiomeR 1.0 is the first stable release. The manuscript describing it is under peer review (link to come).

See [release notes and migration steps](NEWS.md) and the [release convention](RELEASES.md).

## User manual

The user manual is built from the Quarto book in `website/`. It includes a quickstart guide, an output gallery and full implementation details: <https://koefoeden.github.io/multiomeR/>

## System requirements

- Linux system with at least 60 GB of RAM, preferably equipped with a job-scheduler supported by the crew.cluster package: SLURM, PBS, SGE or LSf.

## Citation

If you use multiomeR, please cite the release you used: GitHub's **Cite this repository** gives its authors from [`CITATION.cff`](CITATION.cff), and the release's Zenodo record gives its DOI. Once the manuscript describing multiomeR is published, this section will cite it.

multiomeR is built on targets and BPCells, so please also cite them:

- Landau WM. The targets R package: a dynamic Make-like function-oriented pipeline toolkit for reproducibility and high-performance computing. *Journal of Open Source Software* 2021;6:2959. <https://doi.org/10.21105/joss.02959>
- Parks B, Greenleaf W. Scalable high-performance single cell data analysis with BPCells. *bioRxiv* 2025. <https://doi.org/10.1101/2025.03.27.645853>

Please also cite the methods behind the results you report, such as weighted nearest-neighbor integration or SCAVENGE. The In depth pages of the [user manual](https://koefoeden.github.io/multiomeR/) cite each method where it is used.

## Contributions

Bug reports and broadly useful feature requests are welcome, especially when they affect users analyzing 10x Multiome data. The project prioritizes lean, inspectable workflow changes over broad abstractions or site-specific convenience layers. See `.github/CONTRIBUTING.md`.

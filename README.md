[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

# multiomeR

multiomeR is a targets-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome data. It is designed as a lean, readable framework that users can adapt to their own studies rather than as a black-box command-line pipeline.

The active workflow is a single root `targets` project driven by `_targets.R`, `cfg_reactions.tsv`, `cfg_datasets.yaml`, and `cfg_aggregations.yaml`.

## Status

multiomeR is in beta. The workflow is actively changing, and public interfaces may still change.

## Documentation

The user manual is built from the Quarto book in `website/`:

- Documentation site: <https://koefoeden.github.io/multiomeR/>
- Installation: <https://koefoeden.github.io/multiomeR/installation.html>
- Quickstart: <https://koefoeden.github.io/multiomeR/quickstart.html>
- Configuration: <https://koefoeden.github.io/multiomeR/configuration.html>

## Quickstart

```bash
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR
curl -fsSL https://pixi.sh/install.sh | sh # if pixi is not installed yet
export PATH="$HOME/.pixi/bin:$PATH"
pixi install
cp cfg_datasets_template.yaml cfg_datasets.yaml
cp cfg_aggregations_template.yaml cfg_aggregations.yaml
cp cfg_reactions_template.tsv cfg_reactions.tsv
cp crew_controllers_template.R crew_controllers.R
cp module_differential_analyses/cfg_template.yaml module_differential_analyses/cfg.yaml
cp module_genetic_enrichment/cfg_template.yaml module_genetic_enrichment/cfg.yaml
pixi shell
R
```

Then install the remaining GitHub R packages inside that R session:

```r
pak::pkg_install(c("bnprks/BPCells/r", "stuart-lab/signac@v2", "plger/betterChromVAR"))
```

Then follow the [quickstart](https://koefoeden.github.io/multiomeR/quickstart.html) to download public 10x example data and run the example workflow.

## Requirements

- Linux x86_64.
- `git`.
- `pixi`.
- 10x Genomics Multiome outputs from `cellranger-arc count`.
- Cell Ranger ARC reference directories matching the genomes in your config.
- Sufficient CPU, RAM, and disk for file-backed single-cell analysis.

Most R packages and command-line tools are managed through `pixi.toml` and `pixi.lock`.

## Contributions

Bug reports and broadly useful feature requests are welcome, especially when they affect users analyzing 10x Multiome data. The project prioritizes lean, inspectable workflow changes over broad abstractions or site-specific convenience layers. See `.github/CONTRIBUTING.md`.

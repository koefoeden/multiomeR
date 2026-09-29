# Plan your analysis

Before the first run, prepare the inputs below and choose where your settings will live.

## Inputs to prepare

- **Cell Ranger ARC outputs:** one `cellranger-arc count` output directory per GEM well, with the [files the pipeline reads](reference_GEM_wells.md#cellranger-inputs). GEM wells analyzed together must use the same Cell Ranger ARC reference.
- **Donor metadata table:** a TSV with one row per donor and the phenotypes or covariates you want to analyze; see [Donor metadata table](reference_donor_metadata.md).
- **Donor genotypes (pooled GEM wells only):** a VCF with the genotypes of the pooled donors, whose sample names become the donor IDs; see [Pooled GEM wells](reference_GEM_wells.md#pooled-wells).
- **CellBender output (optional):** to use gene-expression counts corrected for ambient RNA, run [CellBender remove-background](https://cellbender.readthedocs.io/en/latest/usage/) first and give its H5 file in the [GEM well table](reference_GEM_wells.md#cellranger-inputs).
- **Compute resources:** for large datasets, use a compute cluster with a job scheduler; see [Choose where the analysis runs](performance_distributed_computing.md).

## Where the configuration lives {#configuration-directory}

The pipeline reads its settings from one configuration directory:

- `cfg_GEM_wells.tsv`: one row per GEM well; see [GEM well table](reference_GEM_wells.md).
- `cfg_aggregations.yaml`: one entry per aggregation, including the path to its donor metadata table; see [Aggregation configuration](reference_aggregations.md).
- `cfg_module_<module>.yaml`: the settings of one optional module, needed only when that module is enabled.

By default, this is `configuration/`, which holds the public demo and example settings. You can edit it directly, or keep your project's settings in a copy and select that copy in `configuration.local`:

```{.bash filename="Bash"}
cp -r configuration configuration_my_project
echo configuration_my_project > configuration.local
```

`configuration.local` lives in the repository root, is ignored by Git and contains one directory path; a relative path is read from the repository root. While it exists, every configuration file comes from the selected directory, with no fallback to `configuration/`. Delete it to return to the defaults. Do not change the selection while a pipeline is running. Relative paths inside the files are still read from the repository root, and the targets store set in `_targets.yaml` stays the same.

Continue to [Run your own analysis](main_running.md#steps).

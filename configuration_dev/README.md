# Development configuration

This directory holds the larger public-data aggregations behind the output
gallery, the release validation of the optional modules and the manuscript
figures. Select it in `configuration.local`:

```bash
printf '%s\n' 'configuration_dev' > configuration.local
```

- `mixed_human_31x`: 31 public 10x Genomics and ENCODE GEM wells with every
  optional module enabled.
- `comparison_1x` to `comparison_20x`: its first 1, 2, 5, 10 and 20 GEM wells,
  as `comparison_` copies without per-well QC filters or CellBender, for the
  manuscript's Seurat/Signac resource comparison and algorithm comparisons.

The Cell Ranger ARC inputs were reprocessed from the raw reads with Cell Ranger
ARC 2.1.0 and the GRCh38 2024-A reference and are not supplied. The ignored
links `example_data/encode` and `example_data/10x_arc_2.1.0` must point to
them; see `validation/README.md`. The six left-ventricle GEM wells also appear
in `configuration/` for `ENCODE_heart_LV_6x`; keep their rows identical, since
both directories may share one targets store.

`crew_controllers.R` runs every worker as a Slurm job on the node of the
allocation that runs the pipeline, so start the pipeline inside an allocation
on a node with free capacity. Runtimes are then measured on one CPU model, and
the workers' Slurm logs record crew's resource metrics for each target.

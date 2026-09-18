# Configuration

This directory contains the public demo settings and is used by default.

To keep a separate set of settings:

```bash
cp -r configuration configuration_my_project
printf '%s\n' 'configuration_my_project' > configuration.local
```

Edit the files in your chosen directory. Keep `cfg_GEM_wells.tsv`,
`cfg_aggregations.yaml`, and any enabled module settings together. Module files
use flat names such as `cfg_module_differential_analyses.yaml`.
Disabled modules do not need configuration files.

`configuration.local` is ignored by Git and contains one directory path.
Removing it selects `configuration/` again. The pipeline never falls back to
files in another directory. Data paths retain their existing interpretation,
and selecting configuration does not change the targets store. Do not change
the selection during a run.

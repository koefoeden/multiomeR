# Configuration

This directory contains the public example settings. It is selected when the
repository root has no `configuration.local`. To keep a separate set of
settings, copy it and select the copy:

```bash
cp -r configuration configuration_my_project
printf '%s\n' 'configuration_my_project' > configuration.local
```

Edit the files in the copy, and commit them if you share the settings.
`configuration.local` holds one directory path, relative to the checkout root
or absolute. It is local to the checkout and ignored by Git, and removing it
selects `configuration/` again.

All active configuration comes from the selected directory; missing files are
never taken from another directory. Core files are `cfg_GEM_wells.tsv`,
`cfg_aggregations.yaml` and `crew_controllers.R`; module files are flat
`cfg_module_<module>.yaml` files, needed only when that module is enabled.
Schemas live in `manifests/`. In R, `configuration_path("cfg_GEM_wells.tsv")`
shows the resolved file. Data paths keep their existing interpretation, so
moving a configuration does not make them relative to it, and `_targets.yaml`,
not the selection, sets the targets store. Do not change the selection during a
run.

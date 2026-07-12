# Manuscript Figures

These scripts document how selected manuscript and preprint figures are assembled
from multiomeR target outputs. They are included for figure provenance and
methodological transparency, not as guaranteed standalone public demos.

Some scripts depend on target stores, benchmark metadata, raw Cell Ranger output,
or in-house CBMR aggregations that are not distributed with the public repository.
Those scripts should fail loudly when required target outputs or plot-object
sidecars are unavailable.

Generated outputs are written below `manuscript_figures/outputs/` and are not
tracked.

The benchmark renderer caches the portable wall-time result below the active
targets store and writes its plot artifacts below `manuscript_figures/outputs/`:

```bash
pixi run Rscript manuscript_figures/render_benchmark_walltime_plot.R
```

Run that command first from the checkout that owns the targets store and its
private configuration. A different checkout can then render the cached result
by pointing at that store:

```bash
MULTIOMER_TARGETS_STORE=/path/to/targets/store \
  pixi run Rscript manuscript_figures/render_benchmark_walltime_plot.R
```

The combined figure renderer reads the resulting local plot object from
`manuscript_figures/outputs/benchmark/`.

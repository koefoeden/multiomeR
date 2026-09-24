# Local incremental validation

Run on an allocated compute node with the usual Pixi environment:

```sh
pixi run --use-environment-activation-cache validate-local
```

This checks local inputs and external sources, runs unqualified `tar_make()`,
and verifies outputs. It uses the same checkout and configured targets store;
valid targets are reused. Do not run it alongside another pipeline in that store.

The command sets `MULTIOMER_VALIDATION=1` for itself and its children. This selects
exactly `immune_human_2x`, `brain_mouse`, and `mixed_human_31x` plus their
required GEM wells, regardless of their activation flags. Other entries stay
excluded. Without the variable (or with `0`), normal flags apply: only the human
demo is enabled by default. Graph construction prints the validation selection.

For an interactive validation run, set the variable before `tar_make()`:

```r
Sys.setenv(MULTIOMER_VALIDATION = "1")
targets::tar_make()
```

That interactive command runs the pipeline only; use `validate-local` for the
preflight and output checks. Unset the variable to restore normal selection.
No files, branches, or stores are swapped or cleared.

## Inputs

Download the public demos if needed:

```sh
bash example_data/download_10X_cellranger_count_data.sh
```

The ignored link `example_data/encode` must point to a directory containing the
24 ENCODE count directories named by accession. They contain locally processed
public data: ARC 2.1.0, GRCh38 2024-A, and each count's
`outs/cellbender-output_gex_bent_only_filtered.h5`. Raw-read processing is outside
this test. Public donor metadata and ENCODE source URLs are committed.

`mixed_human_31x` also needs the ignored link `example_data/10x_arc_2.1.0`, a
directory with the seven public 10x Genomics datasets it names, reprocessed
from their raw reads in the same way as the ENCODE counts. Its PBMC well
`healthy_PBMC_human_2024A` is therefore separate from the demo's downloaded
`healthy_PBMC_human`, which uses the 10x Genomics outputs as published.

`mixed_human_31x` is the aggregation behind the output gallery and manuscript
figures, with all optional modules enabled. Its models compare male versus female
cardiomyocytes and left-ventricle cell composition, using all eligible features.
Only the four female and two male heart donors have sex recorded, with
imbalanced ages and health status; this is a software test and an unadjusted
association, not a causal sex-effect estimate.

## External checks and results

Preflight makes new, bounded network requests on every invocation:

- Pinned GitHub packages: installed commit identity, namespace loading, and the
  source DESCRIPTION at the pinned commit. This does not reinstall packages.
- Each configured Open Targets dataset: parse the live release listing with the
  same helper as the downloader; retrieve four bytes of one shard and require
  the Parquet header and HTTP partial-content response.
- CollecTRI: download into a temporary directory, verify the pinned checksum,
  and parse it with the pipeline reader.
- AnnotationHub: retrieve SQLite headers from its live metadata database and
  the four configured Ensembl resources. This checks transport and file identity,
  not a complete annotation database download or query. The Hub metadata
  database may be downloaded or refreshed to resolve accessions correctly
  (currently about 131 MB); later checks reuse it when unchanged.

Requests have time and size limits. Failed or unsupported range requests fail
preflight instead of falling back to a large download. Only AnnotationHub
metadata may be refreshed; large analysis resources are
untouched by these probes; the pipeline then downloads the GWAS configured for
`mixed_human_31x`'s genetic enrichment.

Reports under `outputs/validation/<timestamp>/` record the Git commit, dirty
working-tree status, session, preflight results and output summary. `PASS` is
written only after all checks succeed, including graph freshness. A dirty-tree
run is development evidence, not validation of the recorded commit alone.

For a quick external-only check, without running targets:

```sh
pixi run --use-environment-activation-cache Rscript -e 'source("validation/preflight.R"); checks <- run_validation_preflight(); print(checks); stopifnot(all(checks$result == "PASS"))'
```

Use a clean working tree for release evidence. Reuse relies on dependencies being
tracked correctly; these checks do not prove a fresh environment can be installed
or every remote shard is intact. Test clean installation separately when changing
dependencies. Never use `shortcut = TRUE` or force targets current for validation.

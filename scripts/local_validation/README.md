# Local release-candidate validation

An ordinary Esrum host polls public `vX.Y.Z-rc.N` tags and submits two isolated
Slurm jobs per commit. Both execute unqualified `targets::tar_make()`:

- The public demo uses the standard configs, with human and mouse enabled.
- An ENCODE differential-analysis profile replaces the configs only in its
  disposable checkout and provides additional result assertions.

The watcher accepts RC commits reachable from public `main`, records submissions
before queueing, and never automatically retries failures. RC tags are immutable.
No historical RC tags existed when this setup was prepared.

## Install on the host

Use a stable local copy of this directory. Set `MULTIOMER_VALIDATION_PROFILE_REPO`
to a separate, clean, locally committed private Git repository containing:

```
cfg_GEM_wells.tsv
cfg_aggregations.yaml
module_differential_analyses/cfg.yaml
validate.R
```

Configuration inheritance uses the normal pipeline readers. The private
`validate.R` must assert that expected module results exist and contain actual
results. Keep all private configuration and data outside the public checkout.
An optional `example_data/` directory is copied alongside the configs for donor
metadata. The installed profile uses six public ENCODE left-ventricle donors;
its processed ARC and CellBender inputs must already exist on the host.

Run `watch.sh` every five minutes using a user-systemd timer. On Thomas's host,
the declaration belongs in `~/personal/recurring_tasks.yml`; install it with
`sync-recurring-tasks`. The host needs Git, Python 3, gh, and authenticated Slurm
clients. Compute nodes need Pixi and access to the local run directory and data.
The host's user-systemd manager must remain running between logins; enable user
lingering if it is not already configured.
Authenticate gh with repository commit-status write permission. The watcher
keeps this authentication on the host; job submission uses `--export=NONE`.

The default state directory is `~/.local/state/multiomer-validation`. Keep it
private and on storage visible to compute nodes. Each attempt records the exact
public code SHA and private profile SHA, snapshots the harness, and has fresh
pipeline outputs. Pixi environments are cached by dependency specification;
public inputs are cached by download-manifest hash, and AnnotationHub downloads
persist locally. Jobs request
16 CPUs, 256 GB RAM, and 24 hours, matching the standard local Crew setup.

## Results and release gate

The host publishes only fixed, generic status messages under:

- `full-demo/human`
- `full-demo/mouse`
- `local-validation/differential-analysis`

Both demo statuses describe the same complete demo run. Full logs, profile
identifiers, and detailed results remain local. A zero exit code alone cannot
pass: the output assertions and freshness check must also succeed. Cancelled or
failed jobs without a result are reconciled through host-side Slurm accounting.

Before publishing a release, run:

```sh
python3 scripts/local_validation/check_release.py EXACT_40_CHARACTER_COMMIT_SHA
```

This fails unless the latest status for all three contexts is successful.
It does not prevent someone from manually publishing outside this procedure.

To retry a known terminal job, run `watch.sh --retry SHA/demo` or
`watch.sh --retry SHA/private`. A private retry snapshots the current committed
profile. Interrupted submissions with no recorded job ID require manual Slurm
reconciliation before changing local state; they are never silently resubmitted.

Test the coordinator without Slurm or GitHub writes:

```sh
pixi run --use-environment-activation-cache python -m unittest discover -s tests/local_validation -v
```

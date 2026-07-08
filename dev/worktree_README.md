# Pipeline Experiment Worktrees

Use these scripts when you want a sibling git worktree that can run the full
`targets` pipeline against a populated writable output store without risking the
base checkout's `outputs/` during experimentation.

## Motivation

The repository's `outputs/` store is large and contains the current target cache.
A plain fresh worktree is fast for code edits but starts with an empty store, so
full pipeline experiments would recompute too much. Symlinking base `outputs/`
into a worktree is unsafe because `tar_file()` targets and plotting helpers may
overwrite files or directories.

This workflow keeps base outputs isolated:

1. `dev/worktree_setup` refreshes a sibling scratch output store from
   `base/outputs/` with `rsync`.
2. The worktree gets `outputs/` as hardlinks to that scratch store, not to base.
3. The worktree gets `.pixi` as a symlink to the base environment.
4. The worktree gets local runtime paths from the base checkout:
   `example_data` is hardlinked and `crew_controllers.R` is copied as-is.
5. The worktree can run `targets::tar_make()` normally and only rerun outdated
   targets according to its own code and output metadata.
6. Output promotion back to base is a separate guarded step.

If a worktree mutates hardlinked files in place, only the scratch store can be
affected. The base store is touched only by `dev/worktree_promote --apply`.

The runtime paths are marked `skip-worktree` in the experiment worktree when
they are tracked by Git. This keeps local base wiring such as
`crew_controllers.R` faithful without making the experiment branch dirty.

## Create A Worktree

From the base checkout:

```bash
dev/worktree_setup marker-tfs
```

This creates sibling paths such as:

```text
../multiomeR_worktree_marker-tfs
../multiomeR_outputs_scratch_1
```

The label names the worktree and default branch only. Scratch stores are a
reusable pool: setup uses the first existing scratch directory whose active
marker does not point to a live git worktree, regardless of the label. If no
free scratch store exists, it creates the first free numeric scratch path.

Without a label, the first free numeric worktree label is used.

Useful overrides:

```bash
WORKTREE_BRANCH=codex/atac-marker-tf-invalidation \
WORKTREE_PATH=/path/to/worktree \
WORKTREE_SCRATCH=/path/to/scratch \
dev/worktree_setup marker-tfs
```

The scratch refresh runs locally by default with a bandwidth cap:

```bash
rsync --bwlimit=50M ...
```

Tune the cap when needed:

```bash
WORKTREE_RSYNC_BWLIMIT=25M dev/worktree_setup marker-tfs
```

To run the refresh through Slurm instead:

```bash
WORKTREE_RSYNC_RUNNER=srun dev/worktree_setup marker-tfs
```

This uses:

```bash
srun --cpus-per-task=1 --mem=8G --time=04:00:00 rsync --bwlimit=50M ...
```

Change the allocation flags with `WORKTREE_SRUN_ARGS`:

```bash
WORKTREE_RSYNC_RUNNER=srun \
WORKTREE_SRUN_ARGS="--partition=compute --cpus-per-task=1 --mem=16G --time=08:00:00" \
dev/worktree_setup marker-tfs
```

## Run The Pipeline

From the worktree:

```bash
pixi run Rscript - <<'RSCRIPT'
targets::tar_make()
RSCRIPT
```

The worktree uses its own `outputs/` store. The base store is not modified.

## Merge Code

From the worktree:

```bash
dev/worktree_merge --no-ff
```

Extra arguments are passed to `git merge`. This script only handles code. It
does not touch `base/outputs/`.

## Promote Outputs

Promotion is dry-run by default:

```bash
dev/worktree_promote
```

Apply the output sync only after reviewing the dry run:

```bash
dev/worktree_promote --apply
```

Promotion refuses to run if:

- base `outputs/meta/meta` no longer matches the checksum recorded at setup
- base `outputs/meta/process` appears to point to a live target process
- worktree `outputs/meta/process` appears to point to a live target process

The promotion sync excludes `meta/process` and `outputs/.worktree/`. It also
does not promote permissions or directory mtimes, so scratch/worktree directory
metadata does not overwrite shared base-store directory metadata.
After a successful `--apply`, the recorded base metadata checksum is refreshed
so the worktree can be reused for another promotion cycle.

## Convenience Wrapper

From the worktree:

```bash
dev/worktree_finish --apply --no-ff
```

This calls `dev/worktree_merge --no-ff`, then `dev/worktree_promote --apply`.

## Tear Down

After code has been merged and outputs have been promoted, check whether the
linked worktree can be removed:

```bash
dev/worktree_teardown
```

Remove it only after reviewing the dry run:

```bash
dev/worktree_teardown --apply
```

Teardown refuses to run unless:

- the worktree has no uncommitted Git changes
- the experiment branch is already merged into the base checkout
- base and worktree `outputs/meta/meta` match
- the recorded base metadata checksum matches the current base store
- no base or worktree target process appears active

This removes the linked git worktree and releases the scratch active marker. It
does not delete the scratch store itself.

## Safety Notes

- Do not promote outputs while a base pipeline run is active.
- Treat `dev/worktree_promote --apply` and `dev/worktree_teardown --apply` as
  the destructive steps.
- If promotion refuses because base metadata changed, refresh or create a new
  experiment worktree from the new base state.
- `tar_read()` always reads stored state, even when targets are outdated. Use
  `targets::tar_outdated()` or `targets::tar_make()` to check currency.

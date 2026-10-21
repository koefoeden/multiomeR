# multiomeR repository instructions

## Personal instructions

If `$HOME/AGENTS.md` exists, read it before working in this repository. It may
add user- or environment-specific instructions that are intentionally not
tracked here. A missing home-level file is normal.

## Documentation

The manual's sources are the Markdown files under `website/`; the generated
`website/multiomeR-manual-llm.md` gives compact, repository-wide context. Edit
documentation with
[multiomer-edit-documentation](.agents/skills/multiomer-edit-documentation/SKILL.md).

## Repository-specific agent workflows

Task-specific instructions live under `.agents/skills/`, which `.claude/skills`
links for Claude Code. Use the matching `SKILL.md`, especially for:

- running R code or the targets pipeline;
- debugging failed targets;
- validating R helper or target-graph changes;
- estimating which targets a change or run will rebuild;
- editing documentation, including pipeline diagrams;
- creating commits, pull requests or releases.

Always pass `--use-environment-activation-cache` when invoking `pixi run`. A
worktree uses this checkout's environment: run
`pixi run --manifest-path <this checkout>/pixi.toml --as-is <command>` from the
worktree, without the cache flag, which a worktree run would rewrite with its
own paths for every checkout. Pass a task's command, such as
`Rscript tests/testthat.R`, not its name, since a task runs in this checkout.
Pixi 0.82 rejects an environment reached through a worktree's `.pixi` link as
moved. Never delete it as Pixi suggests: the link leads to the shared
environment.

Treat repository-local skills as living instructions. Unless a task is
explicitly read-only, update the skill you are using in the same task whenever
its use reveals stale, changed, ambiguous or missing guidance. Keep revisions
succinct and high-level; omit details reliably implied by the code or existing
instructions.

## Configuration

Configuration follows [configuration/README.md](configuration/README.md); read
settings through `configuration_path()` and resolve the selection during graph
construction.

## Reusable helper source

Generally reusable process, resource, and structured-output helpers live under
`packages/multiomeRCore/R`; [its README](packages/multiomeRCore/README.md)
describes how the pipeline and standalone repositories use them. Do not
duplicate these implementations under the root `R/` directory or in downstream
runtime packages.

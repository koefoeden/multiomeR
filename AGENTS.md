# multiomeR repository instructions

## Personal instructions

If `$HOME/AGENTS.md` exists, read it before working in this repository. It may
add user- or environment-specific instructions that are intentionally not
tracked here. A missing home-level file is normal.

## Documentation entry points

- Start with `website/multiomeR-manual-llm.md` when compact, repository-wide
  documentation context is useful. It combines the user and implementation
  books in authored order and retains source-file provenance comments.
- The canonical documentation sources are the Quarto files under `website/`
  and `website/implementation/`; the LLM-oriented Markdown file is generated.
- Use `.github/CONTRIBUTING.md` for contribution scope and validation guidance.

After changing either documentation book, render both books and refresh the
LLM-oriented export:

```bash
pixi run --use-environment-activation-cache quarto render website
pixi run --use-environment-activation-cache quarto render website/implementation
pixi run --use-environment-activation-cache -e dev export-website-llm-markdown
```

## Repository-specific agent workflows

Task-specific instructions live under `.agents/skills/`. Agents that support
repository skills should use the matching `SKILL.md`, especially for:

- running R code or the targets pipeline;
- debugging failed targets;
- validating R helper or target-graph changes;
- creating commits, pull requests, releases, or pipeline diagrams.

Always pass `--use-environment-activation-cache` when invoking `pixi run`.

Treat repository-local skills as living workflow documentation. Unless a task
is explicitly read-only, revise the relevant skill when its use reveals
materially stale, ambiguous, or missing guidance. Keep revisions succinct and
high-level; omit details reliably implied by the code or existing instructions.

## Reusable helper source

The generally reusable process, resource, and structured-output helpers live
under `packages/multiomeRCore/R`. The multiomeR pipeline sources those files
directly, so ordinary users can edit them as pipeline code without installing
the nested package. Standalone repositories may install the same directory as
the `multiomeRCore` package from a pinned multiomeR Git commit. Do not duplicate
these implementations under the root `R/` directory or in downstream runtime
packages.

---
name: multiomer-create-mermaid-diagram
description: Generate or revise Mermaid dependency diagrams for the active multiomeR targets graph. Use for implementation-book graph updates, new part_of_graph views, or explicitly requested hand-curated workflow overviews.
---

# multiomeR Create Mermaid Diagram

## Manifest-derived diagrams

Read `website/figures/human_curated/README.md` for the tag contract and canonical
regeneration command. Inspect the owning target fragments and keep graph
membership focused on reader-relevant inputs, transformations, checkpoints,
and outputs.

Published website diagrams use an isolated public-demo configuration with the
documented optional examples enabled. Synchronize those artifacts between both
repositories; never overwrite them with diagrams from a private analysis.
Keep analysis-specific graph exports outside `website/`.

After changing membership, regenerate every tagged view and review the changed
`*_v2.mmd` files. Edge order can differ between checkouts: compare sorted lines,
discard reorder-only views, and keep a changed view's committed order so both
repositories hold the same file. Follow the root `AGENTS.md` workflow to render both books and
refresh the Markdown export. Graph generation must not execute targets.

## Hand-curated overviews

Only construct a separate overview manually when the user explicitly wants an
abstraction that the dependency graph cannot express. Use `flowchart TB`, keep
the node set small, and use the canonical theme and legend under
`website/figures/`. Use a filename distinct from generated `*_v2.mmd` files so
regeneration cannot overwrite it.

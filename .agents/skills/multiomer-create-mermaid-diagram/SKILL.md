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

After changing membership, regenerate every tagged view and review the changed
`*_v2.mmd` files. Follow the root `AGENTS.md` workflow to render both books and
refresh the Markdown export. Graph generation must not execute targets.

## Hand-curated overviews

Only construct a separate overview manually when the user explicitly wants an
abstraction that the dependency graph cannot express. Use `flowchart TB`, keep
the node set small, and use the canonical theme and legend under
`website/figures/`. Use a filename distinct from generated `*_v2.mmd` files so
regeneration cannot overwrite it.

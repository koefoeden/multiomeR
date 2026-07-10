---
name: multiomer-create-mermaid-diagram
description: Generate or revise Mermaid dependency diagrams for the active multiomeR targets graph. Use for implementation-book graph updates, new part_of_graph views, or explicitly requested hand-curated workflow overviews.
---

# multiomeR Create Mermaid Diagram


## Manifest-derived diagrams

Implementation-book graphs are generated from target descriptions tagged with
`[part_of_graph:<graph_id>]`.

1. Inspect existing graph IDs and the relevant target fragment:

```bash
rg -n "part_of_graph:" _targets.R extra_targets module_*
```

2. Add or remove graph tags in target descriptions. Keep membership focused on
   reader-relevant inputs, transformations, checkpoints, and outputs.
3. Regenerate every tagged view through the live manifest:

```bash
pixi run Rscript website/figures/human_curated/graphs_v2.R
```

4. Review the changed `website/figures/human_curated/<graph_id>_v2.mmd` files
   and render the implementation book:

```bash
pixi run quarto render website/implementation
```

Do not depend on `website/cache/targets_graphs`; that cache is optional and is
not the source for current implementation-book diagrams.

## Hand-curated overviews

Only construct a separate overview manually when the user explicitly wants an
abstraction that the dependency graph cannot express. Use `flowchart TB`, keep
the node set small, and use the canonical theme and legend under
`website/figures/`. Use a filename distinct from generated `*_v2.mmd` files so
regeneration cannot overwrite it.

## Output Expectations

- Updated target graph tags and generated `*_v2.mmd` files, or a separately
  named hand-curated `.mmd` file.
- A successful implementation-book render.

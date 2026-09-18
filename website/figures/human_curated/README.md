# Implementation graph views

`graphs_v2.R` generates every `*_v2.mmd` view from the active manifest and `[part_of_graph:<graph_id>]` tags in target descriptions. For the shared website, run it in an isolated public checkout configured with the public demo plus the documented optional examples: genetic enrichment and peak–gene correlation on `immune_human_2x`, and differential analyses on `ENCODE_heart_LV_6x` with its GEM wells active. Supply the matching module settings in that isolated configuration. Never switch a running analysis's configuration to make diagrams. From that checkout's root:

```bash
pixi run --use-environment-activation-cache Rscript website/figures/human_curated/graphs_v2.R
```

Edit membership in the owning target fragment, then regenerate all views. Graph IDs contain letters, numbers, and underscores; repeat the tag to include a target in multiple views. A tag selects visible nodes, not targets to execute.

The generator bypasses untagged intermediate nodes, normalizes configured suffixes, and merges duplicate labels. These diagrams explain dependencies; use `targets::tar_network()` for the exact configured graph.

Keep manually drawn conceptual overviews in separately named files. Do not hand-edit generated `*_v2.mmd` files. After changes, render both documentation books and refresh the Markdown export as described in the root `AGENTS.md`.

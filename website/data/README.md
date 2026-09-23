# Shared documentation inputs

`public_defaults/` is a generated snapshot of the public runtime schema and example YAML files. Both repositories render from these same inputs; personal configuration and deployment-specific schemas are not documentation sources.

After changing public defaults, refresh this snapshot from the reviewed public checkout, then synchronize `website/` and the documentation scripts between repositories:

```bash
pixi run --use-environment-activation-cache -e dev python dev/refresh_website_inputs.py --public-repo /path/to/public-checkout
pixi run --use-environment-activation-cache -e dev render-website
pixi run --use-environment-activation-cache -e dev export-website-llm-markdown
```

Graph diagrams are generated from the public demo with the documented optional module examples enabled, using `website/figures/human_curated/graphs_v2.R`. Generate them in an isolated public checkout/configuration and synchronize the resulting diagrams. Do not generate shared diagrams from a personal analysis.

The [output gallery](../gallery.md) is generated from one built aggregation. After its plot targets are rebuilt, refresh every preview and `output_gallery.yaml`:

```bash
pixi run --use-environment-activation-cache refresh-output-gallery mixed_human_31x
```

The command writes one WebP preview per plot target, replaces `gallery_assets/`, and keeps a hand-edited `source_file` while that file still exists. Refresh it instead of editing previews by hand.

## Parameter browser

`website/parameters.html` is a generated, standalone browser with embedded styles, data and JavaScript. Do not edit it directly. After changing the public parameter snapshot, `helpers/parameter_overview.R`, `helpers/parameter_overview.js`, or the stylesheet, regenerate it:

```bash
pixi run --use-environment-activation-cache -e dev render-parameter-overview
```

Link parameter names as `` [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes) `` from pages directly under `website/`. The browser selects the matching workflow tab, clears its search, and opens the linked row. Workflow links use `parameters.html#workflow=aggregation` (or a module scope). The website builder refreshes and publishes the browser automatically; ordinary prose edits do not require regeneration or a book render.

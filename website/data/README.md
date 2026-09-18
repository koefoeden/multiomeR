# Shared documentation inputs

`public_defaults/` is a generated snapshot of the public runtime schema and
example YAML files. Both repositories render from these same inputs; personal
configuration and deployment-specific schemas are not documentation sources.

After changing public defaults, refresh this snapshot from the reviewed public
checkout, then synchronize `website/` and the documentation scripts between
repositories:

```bash
pixi run --use-environment-activation-cache -e dev python dev/refresh_website_inputs.py --public-repo /path/to/public-checkout
pixi run --use-environment-activation-cache -e dev render-website
pixi run --use-environment-activation-cache -e dev export-website-llm-markdown
```

Graph diagrams are generated from the public demo with the documented optional
module examples enabled, using `website/figures/human_curated/graphs_v2.R`.
Generate them in an isolated public checkout/configuration and synchronize the
resulting diagrams. Do not generate shared diagrams from a personal analysis.

Gallery snapshots are checked-in artifacts. When a target or plot changes,
refresh the image from that public target or mark its preview pending; changing
a target label does not make an old image current.

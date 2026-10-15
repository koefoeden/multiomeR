---
name: multiomer-edit-documentation
description: Edit the multiomeR manual, including its methods pages, documentation figures, Mermaid graph views, or generated documentation inputs. Use for changes under website/, the output gallery, the parameter browser, the public-defaults snapshot, or the part_of_graph views of the targets graph.
---

# multiomeR Edit Documentation

## Sources and shared files

- Edit the Markdown sources under `website/`, never the temporary `.qmd` inputs
  the builder prepares; `website/multiomeR-manual-llm.md` is generated.
- The complete `website/` tree, `pipeline_manuscript_figures/` (whose results
  the private manuscript includes) and `dev/render_website.py`,
  `dev/export_website_llm_markdown.py`, `dev/refresh_website_inputs.py`,
  `dev/refresh_output_gallery.py`, `dev/refresh_demo_outputs.R`,
  `dev/check_website_links.py` and `dev/render_parameter_overview.R` must
  remain identical in the public and private repositories.
- Keep private setup and benchmark notes in repository-local guides, never
  under `website/`.

## Write

- Call the core workflow the **primary module** (`primary-module` as a
  modifier); the others are optional modules. Keep page file names and anchors
  stable when only wording changes.
- Use Configure → Run → Review with output-file trees; interpretation belongs
  in the plots.
- Describe what an analysis does in its `website/_shared_methods/` fragment,
  citing `website/references.bib` and stating only key values that are not
  configurable; configurable settings belong to the parameter browser. Keep
  fragments free of links: source links and graphs belong on the
  `website/methods_*.md` pages that include them.
- The manuscript supplement groups the fragments differently from the methods
  pages, so split or regroup them only at paragraph boundaries.

## Regenerate

- Render from the shared public-default snapshot, never a selected personal
  analysis.
- Follow `website/data/README.md` for the defaults snapshot, output gallery,
  demo outputs and parameter browser.
- Generate graph views as `website/figures/human_curated/README.md` describes,
  in an isolated public configuration. Keep their membership focused on
  reader-relevant inputs, transformations, checkpoints and outputs, and keep
  analysis-specific graph exports outside `website/`. Edge order can differ
  between checkouts: compare sorted lines, discard reorder-only views, and keep
  a changed view's committed order so both repositories hold the same file.
- `website/figures/*.drawio.svg` embed their draw.io model. Edit them in draw.io
  and export SVG with the diagram included; a hand edit must change the model
  and the rendered SVG together.

## Validate and commit

1. Check source links, render the manual with the builder, never
   `quarto render website` directly, so executable R chunks still generate
   galleries and diagrams, then check its links. In a worktree, run these
   through `pixi run --as-is`, not an activated shell hook.
   ```bash
   pixi run --use-environment-activation-cache -e dev check-source-links
   pixi run --use-environment-activation-cache -e dev render-website
   pixi run --use-environment-activation-cache -e dev check-website-links
   ```
2. Refresh the LLM export with
   `pixi run --use-environment-activation-cache -e dev export-website-llm-markdown`
   and commit it with the source change that produced it, so later commits do
   not pick up unrelated export changes.
3. Before publishing a shared change, confirm that the shared files are
   identical in both repositories and validate both renders.

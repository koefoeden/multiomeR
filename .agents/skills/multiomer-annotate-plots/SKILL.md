---
name: multiomer-annotate-plots
description: Write or review multiomeR plot titles, subtitles, and captions using three levels of progressive detail. Use when adding plots, revising plot annotations, or checking whether figures explain their interpretation and methods.
---

# Plot annotation specification

Pipeline plots should stand on their own when opened outside the workflow.
Use three levels of progressive detail: identify the subject, guide the reader,
then explain the relevant implementation. This is the repository specification
for new or revised plot annotations; existing plots are examples, not proof of
complete adoption.

## The three levels

| Level | Reader's question | Content |
| --- | --- | --- |
| **Title** | What am I looking at? | A succinct, human-friendly description of the measured quantity and comparison or grouping. Include modality or population when needed to distinguish the output. |
| **Subtitle** | How should I read this, and why does it matter? | The main reading instruction: which pattern to inspect, how to interpret it, and its relevance to the biological question or QC decision. Explain unfamiliar visual encodings when needed. |
| **Caption** | Exactly how was this produced? | Implementation and methodological details that affect meaning: input population, exclusions, aggregation and weighting, transformations, denominators, statistical procedures, thresholds, ordering, and limitations. Include only what is relevant to this plot. |

Titles should usually be a short phrase; subtitles one or two short sentences.
Captions can be longer, but should remain readable at the exported figure size.
These are writing targets, not fixed word limits. Avoid repeating the same
information at all three levels or filling a level with boilerplate.

## Writing decisions

- Use descriptive titles for reusable diagnostics. Assert a finding only when
  it is supported by the actual plotted data; avoid hard-coded conclusions that
  may become false for another aggregation.
- Give the subtitle an interpretation task, such as checking coherent marker
  support or a component dominated by a technical covariate. Explain why that
  pattern matters without implying it automatically proves identity, causality,
  or a filtering decision.
- Keep essential qualifications visible in the subtitle or beside the affected
  encoding when their absence would make the main reading misleading. The
  caption can explain them further. For example, marker-derived labels are not
  independent validation of the same markers.
- Keep units and basic encoding labels on axes and legends. Use annotation text
  for non-obvious meaning, such as what zero represents, whether scaling permits
  comparisons across rows, or which population supplies a percentage denominator.
- In captions, distinguish cells from independent samples, raw from adjusted
  p-values, and descriptive scores from probabilities. State the multiple-testing
  family, weighting, or reference population when relevant. Derive changing
  thresholds and counts from the values used by the plot.
- Prefer meaningful method names over internal object names, target IDs, file
  paths, or configuration keys. Include software/version provenance only when
  it helps explain the output; leave execution logs and exhaustive provenance
  to accompanying reports.
- Use explicit empty-result text when a plot has no eligible data. Distinguish
  insufficient data from an analysis that ran and found no qualifying result.

## Example and existing implementations

An improved wording example inspired by the marker dot plots:

**Title:** Marker expression by assigned cell type

**Subtitle:** Look for several markers supporting each matching cell type and
for broad signal elsewhere that may indicate weak specificity. Labels assigned
from these markers are not independent validation.

**Caption:** Colour shows per-gene scaled mean signal across the plotted groups;
dot area shows the percentage of cells with detected signal within each group.
Shaded boxes mark matches between assigned cell types and marker sets. Only
available positive markers are shown; shared genes can appear in multiple sets.

Before reusing this wording, verify the actual input population, scaling, and
denominator in the owning helper and target.

Useful source examples, each with aspects to learn from rather than copy blindly:

- `style_marker_dot_plot()` in [R/QC_helpers.R](../../../R/QC_helpers.R)
  separates reading instructions from scoring and ordering details.
- `plot_motif_family_accessibility_heatmap()` in the same file documents the
  input population, score aggregation, clustering, and biological limitations.
- PCA diagnostics in
  [extra_targets/GEX_merge_and_dim_reduc_targets.R](../../../extra_targets/GEX_merge_and_dim_reduc_targets.R)
  explain which patterns to inspect and their QC relevance.

## Layout and review

Use the plotting system's native title, subtitle, and caption fields. For
composite plots, share annotations when they apply to all panels and identify
panel-specific differences clearly. Manuscript assemblies may move subtitle
and caption content into the figure legend: preserve its meaning there. Check the assembled legend explicitly when panel annotations are hidden.

When changing annotations, trace the owning helper and target to verify every
claim. Render a representative output at its intended dimensions and check
wrapping, clipping, text size, and agreement with axes and legends. Include an
empty or alternative branch when wording depends on it. Use the repository's
R execution and validation skills for code changes; a specification-only edit
does not require a pipeline run or new tests. Review existing plots as they are
touched unless a broader annotation audit is requested.

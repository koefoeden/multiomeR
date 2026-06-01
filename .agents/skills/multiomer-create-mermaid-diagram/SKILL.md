---
name: multiomer-create-mermaid-diagram
description: Create or update a Mermaid diagram for a multiomeR targets workflow. Use when the user asks to create or revise a pipeline overview diagram for multiomeR.
---

# multiomeR Create Mermaid Diagram


## Diagram formatting rules

- Use `flowchart TB` for a top-to-bottom layout.
- Organize the graph into subgraphs where meaningful, i.e. in case of multiple similar objects of the same type. However, do not over-do this, or create subgraphs for categories that are handled by node styling.
- Node labels should consist of a short, reasonably descriptive label on the first line, and a representative target name (without suffixes) on the second line in parentheses.
- Processing steps (e.g. transformations, tool calls) should always be placed on arrow paths instead of as nodes.
- If multiple upstream objects feed the same tool step, prefer a merge or junction pattern so the tool label appears once. For this, use headless incoming connectors (`---`), so the visual emphasis stays on the shared downstream step rather than on separate directed edges.
- If a step and its associated nodes are optional, use dashed connectors and dashed node borders.

## Workflow to create/update the diagram

1. Read the specific workflow structure by inspecting the relevant file at `website/cache/targets_graphs/<graph_name>_<target_graph_match>_targets_graph.csv`, which contains a table of the edges of the targets graph. For current multiomeR work, prefer graph names based on the root workflow or module name rather than old targets project names.
2. Determine the high-level, abstract phases of the pipeline that help reader orientation, i.e. you shouldn't make a detailed replication of the target graph. In this process, create a node-inventory list based on all the relevant input, intermediate and output-objects that will exist in this graph. Assign each one a simple, alphanumeric ID and a clear descriptive label, e.g. ID [Descriptive Label]
3. Using plain English, list out the exact chronological flow of how these nodes connect, from start to finish. Include the text that should appear on the connection arrows, if any, e.g. [ID] connects to [ID] via "Arrow Label"
4. Briefly review your Edge Mapping and ensure there are no contradictory directional flows that would break a strict Top-Down layout.
5. Finally, create or update the `website/figures/<graph_name>_<target_graph_match>_overview.mmd` file with the workflow-specific structure only, while using the standard color scheme for node styling as shown in the `assets/standard_node_color_legend.mmd` file. Note that a helper file provides shared Mermaid theme settings, so do not add these.

## Output Expectations

When using this skill, produce:
- A new or updated Mermaid `.mmd` diagram that is renderable as-is.

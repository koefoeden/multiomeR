#!/usr/bin/env python3
"""Render and audit simplified Mermaid diagrams from YAML target-graph recipes."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict, deque
from pathlib import Path

import yaml


GRAPH_DIR = Path("website/cache/targets_graphs")
FIGURE_DIR = Path("website/figures")
RECIPE_DIR = Path("website/diagram_recipes")

CLASS_DEFS = [
    "  classDef regular_input fill:#dbeafe,stroke:#3b82f6,color:#1e3a8a;",
    "  classDef CFG_input fill:#f3e8ff,stroke:#8b5cf6,color:#6d28d9;",
    "  classDef intermediate_object fill:#f8fafc,stroke:#64748b,color:#1f2937;",
    "  classDef regular_output fill:#ecfdf5,stroke:#10b981,color:#065f46;",
    "  classDef QC_output fill:#fff7ed,stroke:#f97316,color:#9a3412;",
    "  classDef optional_regular_input fill:#dbeafe,stroke:#3b82f6,stroke-dasharray: 5 5,color:#1e3a8a;",
    "  classDef optional_intermediate_object fill:#f8fafc,stroke:#64748b,stroke-dasharray: 5 5,color:#1f2937;",
    "  classDef optional_regular_output fill:#ecfdf5,stroke:#10b981,stroke-dasharray: 5 5,color:#065f46;",
    "  classDef junction_node fill:none,stroke:none,color:#334155,font-weight:bold;",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_recipe(path: Path) -> dict:
    with path.open() as handle:
        recipe = yaml.safe_load(handle)
    recipe["recipe_path"] = str(path)
    return recipe


def matches_any(text: str, patterns: list[str]) -> bool:
    return any(re.search(pattern, text) for pattern in patterns)


def targets_for_node(names: list[str], node: dict) -> list[str]:
    return [name for name in names if matches_any(name, node["target_patterns"])]


def has_path(start_targets: list[str], end_targets: list[str], adjacency: dict[str, list[str]]) -> bool:
    if not start_targets or not end_targets:
        return False
    ends = set(end_targets)
    seen = set(start_targets)
    queue = deque(start_targets)
    while queue:
        current = queue.popleft()
        if current in ends:
            return True
        for child in adjacency.get(current, []):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return False


def mermaid_label(label: str, representative: str) -> str:
    return f"{label}<br/>({representative})"


def write_mermaid(recipe: dict) -> Path:
    output = FIGURE_DIR / f"{recipe['output_stem']}_overview.mmd"
    lines = ["%%{init: {'layout': 'elk'}}%%", "flowchart TB", *CLASS_DEFS, ""]
    for node in recipe["nodes"]:
        label = mermaid_label(node["label"], node["representative"])
        lines.append(f"  {node['id']}(\"{label}\")")
    lines.append("")

    for edge in recipe["edges"]:
        style = edge.get("style", "-->")
        label = edge.get("label", "")
        if label:
            lines.append(f"  {edge['from']} {style}|{label}| {edge['to']}")
        else:
            lines.append(f"  {edge['from']} {style} {edge['to']}")
    lines.append("")

    by_class = defaultdict(list)
    for node in recipe["nodes"]:
        by_class[node["class"]].append(node["id"])
    for class_name, node_ids in by_class.items():
        lines.append(f"  class {','.join(node_ids)} {class_name};")

    output.write_text("\n".join(lines) + "\n")
    return output


def write_mapping(recipe: dict, target_map: dict[str, list[str]]) -> Path:
    output = GRAPH_DIR / f"{recipe['graph_stem']}_{recipe['slug']}_collapsed_mapping.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("diagram", "diagram_node_id", "diagram_label", "represented_target"))
        writer.writeheader()
        for node in recipe["nodes"]:
            targets = target_map[node["id"]]
            if not targets:
                writer.writerow(
                    {
                        "diagram": recipe["slug"],
                        "diagram_node_id": node["id"],
                        "diagram_label": node["label"],
                        "represented_target": "",
                    }
                )
            for target in targets:
                writer.writerow(
                    {
                        "diagram": recipe["slug"],
                        "diagram_node_id": node["id"],
                        "diagram_label": node["label"],
                        "represented_target": target,
                    }
                )
    return output


def render_recipe(recipe: dict, graph_cache: dict[str, dict]) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict]]:
    if recipe["graph_stem"] not in graph_cache:
        graph_cache[recipe["graph_stem"]] = load_graph(recipe)
    graph = graph_cache[recipe["graph_stem"]]
    names = graph["names"]
    adjacency = graph["adjacency"]

    target_map = {node["id"]: targets_for_node(names, node) for node in recipe["nodes"]}
    assigned_nodes = defaultdict(list)
    for node_id, targets in target_map.items():
        for target in targets:
            assigned_nodes[target].append(node_id)
    represented = set(target for targets in target_map.values() for target in targets)
    scoped = {name for name in names if matches_any(name, recipe["scope_patterns"])}
    omitted = sorted(scoped - represented)
    extra = sorted(represented - scoped)
    duplicated = {target: node_ids for target, node_ids in sorted(assigned_nodes.items()) if len(node_ids) > 1}
    unsupported = []
    edge_rows = []

    for edge in recipe["edges"]:
        supported = has_path(target_map[edge["from"]], target_map[edge["to"]], adjacency)
        if not supported:
            unsupported.append(f"{edge['from']}->{edge['to']}")
        edge_rows.append(
            {
                "diagram": recipe["slug"],
                "from": edge["from"],
                "to": edge["to"],
                "label": edge.get("label", ""),
                "supported_by_original_path": supported,
            }
        )

    write_mermaid(recipe)
    write_mapping(recipe, target_map)

    open_questions = ""
    if omitted or extra or duplicated or unsupported:
        open_questions = "Review omitted targets, extra captures, duplicate assignments, and unsupported conceptual edges before treating as final documentation."
    review_rows = [
        {
            "diagram": recipe["slug"],
            "recipe_path": recipe["recipe_path"],
            "represented_targets_n": len(represented),
            "scope_targets_n": len(scoped),
            "omitted_targets_n": len(omitted),
            "extra_targets_n": len(extra),
            "duplicated_targets_n": len(duplicated),
            "unsupported_edges": ";".join(unsupported),
            "open_questions": open_questions,
        }
    ]
    omitted_rows = [{"diagram": recipe["slug"], "omitted_target": target} for target in omitted]
    extra_rows = [{"diagram": recipe["slug"], "extra_target": target} for target in extra]
    duplicate_rows = [
        {"diagram": recipe["slug"], "duplicated_target": target, "diagram_node_ids": ";".join(node_ids)}
        for target, node_ids in duplicated.items()
    ]
    return review_rows, edge_rows, omitted_rows, extra_rows, duplicate_rows


def graph_paths(graph_stem: str) -> dict[str, Path]:
    return {
        "nodes": GRAPH_DIR / f"{graph_stem}_targets_nodes.csv",
        "edges": GRAPH_DIR / f"{graph_stem}_targets_graph.csv",
    }


def r_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def ensure_graph_assets(recipe: dict) -> None:
    paths = graph_paths(recipe["graph_stem"])
    if paths["nodes"].exists() and paths["edges"].exists():
        return

    graph_build = recipe.get("graph_build")
    if not graph_build:
        raise FileNotFoundError(f"Missing graph cache for {recipe['graph_stem']} and no graph_build block in {recipe['recipe_path']}")

    allow_match = graph_build.get("allow_match", graph_build["target_graph_match"])
    force = "TRUE" if graph_build.get("force", False) else "FALSE"
    r_code = "\n".join(
        [
            'source("R/targets_graph_helpers.R")',
            "invisible(build_targets_graph_assets(",
            f"  pipeline_name = {r_string(graph_build['pipeline_name'])},",
            f"  target_graph_match = {r_string(graph_build['target_graph_match'])},",
            f"  allow_match = {r_string(allow_match)},",
            f"  force = {force}",
            "))",
        ]
    )
    subprocess.run(["pixi", "run", "Rscript", "-e", r_code], check=True)

    if not paths["nodes"].exists() or not paths["edges"].exists():
        raise FileNotFoundError(f"graph_build did not create expected graph cache for {recipe['graph_stem']}")


def load_graph(recipe: dict) -> dict:
    ensure_graph_assets(recipe)
    paths = graph_paths(recipe["graph_stem"])
    nodes = read_csv(paths["nodes"])
    edges = read_csv(paths["edges"])
    adjacency: dict[str, list[str]] = defaultdict(list)
    for edge in edges:
        adjacency[edge["from"]].append(edge["to"])
    return {"names": [row["name"] for row in nodes], "adjacency": adjacency}


def write_review_outputs(
    graph_stem: str,
    review_rows: list[dict],
    edge_rows: list[dict],
    omitted_rows: list[dict],
    extra_rows: list[dict],
    duplicate_rows: list[dict],
) -> None:
    with (GRAPH_DIR / f"{graph_stem}_simplified_diagram_review.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=review_rows[0].keys())
        writer.writeheader()
        writer.writerows(review_rows)
    with (GRAPH_DIR / f"{graph_stem}_simplified_diagram_edge_support.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=edge_rows[0].keys())
        writer.writeheader()
        writer.writerows(edge_rows)
    with (GRAPH_DIR / f"{graph_stem}_simplified_diagram_omitted_targets.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("diagram", "omitted_target"))
        writer.writeheader()
        writer.writerows(omitted_rows)
    with (GRAPH_DIR / f"{graph_stem}_simplified_diagram_extra_targets.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("diagram", "extra_target"))
        writer.writeheader()
        writer.writerows(extra_rows)
    with (GRAPH_DIR / f"{graph_stem}_simplified_diagram_duplicate_targets.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("diagram", "duplicated_target", "diagram_node_ids"))
        writer.writeheader()
        writer.writerows(duplicate_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recipes", nargs="*", type=Path, help="Recipe YAML files. Defaults to website/diagram_recipes/*.yaml.")
    parser.add_argument("--strict", action="store_true", help="Fail if any recipe has omitted targets, extra captures, duplicate assignments, or unsupported edges.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    recipe_paths = args.recipes or sorted(RECIPE_DIR.glob("*.yaml"))
    recipes = [read_recipe(path) for path in recipe_paths]
    graph_cache: dict[str, dict] = {}
    rows_by_graph: dict[str, dict[str, list[dict]]] = defaultdict(lambda: {"review": [], "edges": [], "omitted": [], "extra": [], "duplicates": []})

    for recipe in recipes:
        review_rows, edge_rows, omitted_rows, extra_rows, duplicate_rows = render_recipe(recipe, graph_cache)
        rows = rows_by_graph[recipe["graph_stem"]]
        rows["review"].extend(review_rows)
        rows["edges"].extend(edge_rows)
        rows["omitted"].extend(omitted_rows)
        rows["extra"].extend(extra_rows)
        rows["duplicates"].extend(duplicate_rows)

    for graph_stem, rows in rows_by_graph.items():
        write_review_outputs(graph_stem, rows["review"], rows["edges"], rows["omitted"], rows["extra"], rows["duplicates"])

    if args.strict:
        failed = [
            row["diagram"]
            for rows in rows_by_graph.values()
            for row in rows["review"]
            if row["omitted_targets_n"] or row["extra_targets_n"] or row["duplicated_targets_n"] or row["unsupported_edges"]
        ]
        if failed:
            raise SystemExit("Diagram recipe audit failed for: " + ", ".join(failed))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Export the Quarto website books as one LLM-friendly Markdown file."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml


INCLUDE_RE = re.compile(r"\{\{<\s*include\s+([^>\s]+)\s*>\}\}")
CHUNK_RE = re.compile(r"```{r[^}]*}\n.*?\n```", re.DOTALL)
IMAGE_RE = re.compile(r"!\[([^\]]*)\]\(([^)]+)\)(\{[^}]*\})?")
MERMAID_RE = re.compile(r'emit_mermaid\("([^"]+)"\)')
INLINE_R_LINK_RE = re.compile(r"\[([^\]]+)\]\(`r [^)]+`\)")
MAX_PLACEHOLDER_TEXT = 120


@dataclass(frozen=True)
class Page:
    path: Path


@dataclass(frozen=True)
class ExportOptions:
    include_mermaid: bool
    include_generated_chunks: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect website/ and website/implementation/ Quarto pages into one Markdown file."
    )
    parser.add_argument("--site-dir", type=Path, default=Path("website"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("website/multiomeR-manual-llm.md"),
    )
    parser.add_argument(
        "--no-orphans",
        action="store_true",
        help="Do not append a list of tracked QMD pages outside the book graph.",
    )
    parser.add_argument(
        "--include-mermaid",
        action="store_true",
        help="Inline Mermaid graph files instead of compact graph placeholders.",
    )
    parser.add_argument(
        "--include-generated-chunks",
        action="store_true",
        help="Keep unevaluated Quarto helper chunks that generate HTML/Markdown.",
    )
    return parser.parse_args()


def read_book_config(book_dir: Path) -> dict:
    config_path = book_dir / "_quarto.yml"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing Quarto config: {config_path}")
    return yaml.safe_load(config_path.read_text(encoding="utf-8"))


def walk_chapters(book_dir: Path, chapters: list) -> list[tuple[str, str | Page | Path]]:
    entries: list[tuple[str, str | Page | Path]] = []
    for item in chapters:
        if isinstance(item, str):
            entries.append(("page", Page(book_dir / item)))
            continue

        if not isinstance(item, dict):
            continue

        part = item.get("part")
        if part:
            entries.append(("part", str(part)))
            entries.extend(walk_chapters(book_dir, item.get("chapters", [])))
            continue

        href = item.get("href")
        if href:
            target = book_dir / href
            if target.is_dir():
                entries.append(("book", target))
            else:
                entries.append(("page", Page(target)))
            continue

        entries.extend(walk_chapters(book_dir, item.get("chapters", [])))

    return entries


def render_book(
    book_dir: Path,
    repo_root: Path,
    used_sources: set[Path],
    options: ExportOptions,
) -> str:
    config = read_book_config(book_dir)
    title = config.get("book", {}).get("title", book_dir.name)
    chapters = config.get("book", {}).get("chapters", [])
    out = [f"\n# Book: {title}\n"]

    for kind, value in walk_chapters(book_dir, chapters):
        if kind == "part":
            out.append(f"\n## Part: {value}\n")
        elif kind == "book":
            out.append(render_book(Path(value), repo_root, used_sources, options))
        elif kind == "page":
            page = value
            assert isinstance(page, Page)
            out.append(render_page(page.path, repo_root, used_sources, options))

    return "\n".join(out)


def render_page(
    path: Path,
    repo_root: Path,
    used_sources: set[Path],
    options: ExportOptions,
) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Book references missing page: {path}")

    rel = path.relative_to(repo_root)
    used_sources.add(path.resolve())
    body = expand_includes(path, repo_root, used_sources, ())
    body = transform_chunks(body, path, repo_root, options)
    body = transform_images(body)
    body = transform_inline_r_links(body)
    body = body.strip()
    return f"\n<!-- source: {rel} -->\n\n{body}\n"


def expand_includes(
    path: Path,
    repo_root: Path,
    used_sources: set[Path],
    stack: tuple[Path, ...],
) -> str:
    if path in stack:
        chain = " -> ".join(str(p) for p in (*stack, path))
        raise RuntimeError(f"Recursive Quarto include detected: {chain}")

    text = path.read_text(encoding="utf-8")

    def replace(match: re.Match[str]) -> str:
        include_path = (path.parent / match.group(1)).resolve()
        rel = include_path.relative_to(repo_root)
        used_sources.add(include_path)
        included = expand_includes(include_path, repo_root, used_sources, (*stack, path))
        return f"\n<!-- begin include: {rel} -->\n\n{included.strip()}\n\n<!-- end include: {rel} -->\n"

    return INCLUDE_RE.sub(replace, text)


def transform_chunks(
    markdown: str,
    page_path: Path,
    repo_root: Path,
    options: ExportOptions,
) -> str:
    def replace(match: re.Match[str]) -> str:
        chunk = match.group(0)
        header = chunk.split("\n", 1)[0]
        mermaid = MERMAID_RE.search(chunk)

        if mermaid:
            mermaid_path = resolve_asset_path(mermaid.group(1), page_path, repo_root)
            rel = mermaid_path.relative_to(repo_root)
            if not options.include_mermaid:
                return f"[Mermaid graph omitted; source: `{rel}`]"
            code = mermaid_path.read_text(encoding="utf-8").strip()
            return f"<!-- mermaid source: {rel} -->\n\n```mermaid\n{code}\n```"

        if re.search(r"\binclude\s*=\s*FALSE\b", header, re.IGNORECASE):
            return ""

        if is_generated_chunk(header) and not options.include_generated_chunks:
            return f"[Generated Quarto chunk omitted: `{summarize_chunk(chunk)}`]"

        return chunk

    return CHUNK_RE.sub(replace, markdown)


def is_generated_chunk(header: str) -> bool:
    return bool(
        re.search(r"\becho\s*=\s*FALSE\b", header, re.IGNORECASE)
        and re.search(r'\bresults\s*=\s*["\']asis["\']', header, re.IGNORECASE)
    )


def summarize_chunk(chunk: str) -> str:
    lines = chunk.splitlines()[1:-1]
    body = " ".join(line.strip() for line in lines if line.strip())
    return truncate(body)


def resolve_asset_path(asset: str, page_path: Path, repo_root: Path) -> Path:
    candidates = [
        repo_root / asset,
        page_path.parent / asset,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"Could not resolve asset from {page_path}: {asset}")


def transform_images(markdown: str) -> str:
    def replace(match: re.Match[str]) -> str:
        alt, target, _attrs = match.groups()
        label = truncate(alt) if alt else "Image"
        return f"[Image omitted; source: `{target}`; alt: {label}]"

    return IMAGE_RE.sub(replace, markdown)


def transform_inline_r_links(markdown: str) -> str:
    return INLINE_R_LINK_RE.sub(r"\1", markdown)


def truncate(text: str) -> str:
    if len(text) <= MAX_PLACEHOLDER_TEXT:
        return text
    return text[: MAX_PLACEHOLDER_TEXT - 3].rstrip() + "..."


def tracked_qmd_files(site_dir: Path) -> set[Path]:
    try:
        result = subprocess.run(
            ["git", "ls-files", str(site_dir)],
            check=True,
            text=True,
            capture_output=True,
        )
        paths = [Path(line) for line in result.stdout.splitlines()]
    except (subprocess.CalledProcessError, FileNotFoundError):
        paths = list(site_dir.rglob("*.qmd"))

    return {path.resolve() for path in paths if path.suffix == ".qmd"}


def render_orphans(site_dir: Path, repo_root: Path, used_sources: set[Path]) -> str:
    orphans = sorted(tracked_qmd_files(site_dir) - used_sources)
    if not orphans:
        return "\n# Orphaned QMD Pages\n\nNone.\n"

    lines = ["\n# Orphaned QMD Pages\n"]
    lines.append("Tracked QMD files not reached from the Quarto book graph or include graph.\n")
    for path in orphans:
        lines.append(f"- `{path.relative_to(repo_root)}`")
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd().resolve()
    site_dir = (repo_root / args.site_dir).resolve()
    output = (repo_root / args.output).resolve()
    used_sources: set[Path] = set()
    options = ExportOptions(
        include_mermaid=args.include_mermaid,
        include_generated_chunks=args.include_generated_chunks,
    )

    text = [
        "# multiomeR Website LLM Export\n",
        "This file is generated from the Quarto book outlines and resolves Quarto include shortcodes.",
        "Hidden setup chunks, generated helper chunks, Mermaid graph bodies, and verbose image metadata are omitted by default.\n",
        render_book(site_dir, repo_root, used_sources, options),
    ]

    if not args.no_orphans:
        text.append(render_orphans(site_dir, repo_root, used_sources))

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(text).rstrip() + "\n", encoding="utf-8")
    print(output.relative_to(repo_root))
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Render Markdown book sources through temporary Quarto execution inputs."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    subprocess.run(["Rscript", "dev/render_parameter_overview.R"], cwd=root, check=True,
                   env={**os.environ, "R_PROFILE_USER": os.devnull})
    with tempfile.TemporaryDirectory(prefix="multiomer-docs-") as directory:
        stage = Path(directory)
        for path in root.iterdir():
            if path.name not in {"website", "docs", "configuration.local"}:
                (stage / path.name).symlink_to(path, target_is_directory=path.is_dir())
        site = stage / "website"
        shutil.copytree(root / "website", site,
                        ignore=shutil.ignore_patterns(".quarto", "*_files", "*.knit.md"))
        sources = [path for path in site.rglob("*.md")
                   if path.name != "multiomeR-manual-llm.md"]
        names = {path.name for path in sources}
        for path in [*sources, *site.rglob("_quarto.yml")]:
            text = path.read_text()
            for name in names:
                text = text.replace(name, name.removesuffix(".md") + ".qmd")
            path.write_text(text)
        for path in sources:
            path.rename(path.with_suffix(".qmd"))
        for book in ("website", "website/implementation"):
            subprocess.run(
                ["quarto", "render", book], cwd=stage, check=True,
                env={**os.environ, "R_PROFILE_USER": os.devnull},
            )
        # Source links must point to the editable Markdown files on GitHub.
        for path in (stage / "docs").rglob("*.html"):
            path.write_text(path.read_text().replace(".qmd", ".md"))
        shutil.rmtree(root / "docs", ignore_errors=True)
        shutil.copytree(stage / "docs", root / "docs")


if __name__ == "__main__":
    main()

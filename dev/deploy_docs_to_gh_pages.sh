#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

docs_dir="$repo_root/docs"

if [[ ! -d "$docs_dir" ]]; then
  echo "docs/ does not exist. Build the site before deploying." >&2
  exit 1
fi

if [[ ! -f "$docs_dir/index.html" ]]; then
  echo "docs/ does not look like a built site (missing docs/index.html)." >&2
  exit 1
fi

current_branch="$(git branch --show-current)"
if [[ "$current_branch" == "gh-pages" ]]; then
  echo "Refusing to run from the gh-pages branch." >&2
  exit 1
fi

tmp_dir="$(mktemp -d)"
cleanup() {
  git worktree remove --force "$tmp_dir" >/dev/null 2>&1 || true
}
trap cleanup EXIT

if git ls-remote --exit-code --heads origin gh-pages >/dev/null 2>&1; then
  git fetch origin gh-pages
  git worktree add "$tmp_dir" --track -B gh-pages origin/gh-pages
else
  git worktree add "$tmp_dir" --detach
  (
    cd "$tmp_dir"
    git checkout --orphan gh-pages
    git rm -rf . >/dev/null 2>&1 || true
  )
fi

find "$tmp_dir" -mindepth 1 -maxdepth 1 \
  ! -name '.git' \
  -exec rm -rf {} +

rsync -a --delete \
  --exclude '.git' \
  --exclude '.git/' \
  "$docs_dir"/ "$tmp_dir"/

touch "$tmp_dir/.nojekyll"

(
  cd "$tmp_dir"
  git add -A

  if git diff --cached --quiet; then
    echo "No changes to deploy."
    exit 0
  fi

  git commit -m "Deploy docs site"
  git push origin gh-pages
)

echo "Deployed docs/ to gh-pages."

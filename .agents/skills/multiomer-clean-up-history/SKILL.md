---
name: multiomer-clean-up-history
description: Clean up and rewrite feature branch history in the multiomeR repository. Use when the user asks to regroup commits, squash noisy history, rewrite a branch before opening a pull request, or make commit structure/messages more coherent. Especially useful when targets changes, docs changes, generated files, and workflow tweaks were mixed together and the branch should be rebuilt into a smaller set of atomic commits.
---

# multiomeR Clean Up History

## Workflow

Use this skill when the goal is to improve branch history, not to change the final branch content.

Before rewriting:
- Ensure the worktree is clean with `git status --short --branch`.
- If there are unrelated tracked changes, stop and ask the user how to proceed.
- Identify untracked artifacts separately; do not let them influence commit
  grouping unless the branch intentionally adds them.

## Inspect The Branch

Start from the merge-base against the intended target branch, usually `origin/main`:

```bash
git merge-base HEAD origin/main
git log --oneline --reverse <base>..HEAD
git log --reverse --format='COMMIT %h%nSUBJECT %s%n' <base>..HEAD --name-only
git diff --stat <base>..HEAD
```

Use that survey to identify:
- coherent topics that should become separate commits
- files touched by multiple topics
- derived outputs that should be dropped or regenerated
- whether the branch is simple enough for `git rebase -i`

## Choose A Rewrite Strategy

Use interactive rebase only when:
- the existing commits are already mostly atomic
- the task is mostly `reword`, `fixup`, `squash`, or small reordering
- the same files are not heavily entangled across unrelated commits

Prefer rebuild-from-base when:
- the branch mixes unrelated topics
- many commits are tiny fixups
- a few files span several topics
- generated docs or outputs were committed and later removed

Choose from the branch evidence; do not rebuild a branch that only needs a few
fixups or rewords.

## Rebuild From Base

Only after the user has requested a rewrite, the worktree is clean, and the
rebuild strategy is selected, create a safety ref before resetting:

```bash
git branch backup/<branch>-pre-rewrite-<yyyymmdd> HEAD
git reset --hard <base>
```

Then reconstruct a smaller commit stack from the saved branch state:

```bash
git restore --source=backup/<branch>-pre-rewrite-<yyyymmdd> --staged --worktree -- <paths...>
git commit -F- <<'EOF'
<message>
EOF
```

Guidelines:
- Group by user-facing change or durable topic, not by original commit order.
- Keep workflow/skill/docs metadata separate from pipeline behavior changes.
- Keep `website/multiomeR-manual-llm.md` with the documentation source changes
  that regenerate it. Do not commit ignored `docs/` render output.
- Keep deployment helpers separate from the content they deploy.
- Keep plot-only targets as `contained_target_breaking` commits when possible.

If a file belongs to more than one new commit:
- Restore an intermediate revision for the earlier commit.
- Restore the final branch revision for the later commit.

Example:

```bash
git restore --source=<intermediate-commit> --staged --worktree -- R/general_helpers.R
git commit ...
git restore --source=<backup-branch> --staged --worktree -- R/general_helpers.R
git commit ...
```

## Write Correct Commit Messages

For every rebuilt commit:
- Use the format in the `multiomer-git-commit-format` skill.
- Add the mandatory impact keyword line described in `multiomer-impact-keyword-lines`.
- Prefer one of:
  - `non_target_breaking`
  - `contained_target_breaking: <target>`
  - `cascading_target_breaking: <target> [<dataset-scope>]`

Do not guess lazily. If the commit changes pipeline behavior, inspect the relevant target names first.

Helpful checks:

```bash
rg -n "tar_target\\(|tarchetypes::tar_" _targets.R extra_targets module_*
rg -n "<target-or-helper-name>" _targets.R R extra_targets module_*
```

## Verify The Rewrite

After rebuilding the branch:

```bash
git diff --stat backup/<branch>-pre-rewrite-<yyyymmdd>..HEAD
git log --oneline --decorate --graph <base>..HEAD
git status --short --branch
```

The diff against the backup branch should be empty. If it is not empty, the rewrite changed content rather than just history.

## Finish

Report:
- the new commit stack
- whether the final tree matches the original branch
- whether tests were run
- the backup branch name

Only suggest force-pushing after verification:

```bash
git push --force-with-lease origin <branch>
```

## Heuristics

- Prefer the smallest commit stack that still separates durable topics.
- Separate pipeline code, docs cleanup, deployment tooling, and skill/workflow docs unless they are inseparable.
- Prefer a slightly larger coherent commit over a misleadingly narrow one.
- Never delete the backup branch until the user is satisfied with the rewritten history.

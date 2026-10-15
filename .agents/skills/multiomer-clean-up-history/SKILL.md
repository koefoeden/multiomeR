---
name: multiomer-clean-up-history
description: Clean up and rewrite feature branch history in the multiomeR repository. Use when the user asks to regroup commits, squash noisy history, rewrite a branch before opening a pull request, or make commit structure/messages more coherent. Especially useful when targets changes, docs changes, generated files, and workflow tweaks were mixed together and the branch should be rebuilt into a smaller set of atomic commits.
---

# multiomeR Clean Up History

Rewrite only an unpublished feature branch, and only to improve its history:
the final tree must not change. Never rewrite `main` or commits another
repository's trailers name.

## Safety

- Start from a clean worktree; ask the user about unrelated tracked changes.
  Keep untracked artifacts out of the grouping unless the branch adds them.
- Before rewriting, save a backup ref:
  `git branch backup/<branch>-pre-rewrite-<yyyymmdd> HEAD`.
- Afterwards, `git diff --stat backup/<branch>-pre-rewrite-<yyyymmdd>..HEAD`
  must be empty. Keep the backup until the user is satisfied, and suggest
  `git push --force-with-lease origin <branch>` only after this check.
- Agent shells cannot run interactive commands. For fixups and rewords, use
  `git commit --fixup` and `GIT_SEQUENCE_EDITOR=true git rebase -i --autosquash <base>`.
  For entangled branches, reset to the merge base and rebuild commits with
  `git restore --source=<backup> --staged --worktree -- <paths>`, restoring an
  intermediate revision for a file shared by an earlier commit.

## Grouping

- Group by user-facing change or durable topic, not by original commit order;
  prefer the smallest stack that still separates durable topics.
- Keep pipeline behavior, documentation, deployment tooling and skill or
  workflow changes in separate commits unless they are inseparable.
- Keep `website/multiomeR-manual-llm.md` with the documentation source changes
  that regenerate it. Do not commit ignored `docs/` render output.
- Keep plot-only target changes as `contained_target_breaking` commits when
  possible.

Write every rebuilt commit with `multiomer-git-commit-format`, plus any
classification the repository's AGENTS.md requires. Classify the impact of
commits that change pipeline behavior with `multiomer-estimate-rebuild-impact`.

Report the new commit stack, whether the final tree matches the backup, which
tests ran, and the backup branch name.

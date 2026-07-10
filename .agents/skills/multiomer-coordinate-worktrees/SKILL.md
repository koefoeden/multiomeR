---
name: multiomer-coordinate-worktrees
description: Create and coordinate multiple user-owned Codex App tasks for multiomeR in separate App-managed worktrees. Use when the user wants parallel, sidebar-visible worker tasks with isolated codex branches, selectable read-only or full multiomeR environments, task links, progress checks, or worktree handoffs.
---

# Coordinate multiomeR Worktrees

Create top-level Codex App tasks that the user can track individually. Do not use
collaboration subagents, fork the coordinator task, or create Git worktrees
manually for this workflow.

Accept a numbered `Tasks:` list plus any shared or task-specific constraints. If
no concrete tasks are supplied, request the task list before dispatching.

## Select the worktree mode

Support these canonical modes per task:

- `read-only` (default): link `.pixi` and a read-only targets output snapshot for
  pipeline introspection without output writes.
- `full`: link `.pixi` and provision an independent writable scratch `outputs`
  store for pipeline execution.

Accept clear aliases such as `introspection` for `read-only` and `regular` or
`multiomeR` for `full`, but report the canonical mode. Apply a shared mode to all
tasks unless a task overrides it. Use `read-only` when no mode is specified. Ask
for clarification rather than inventing behavior for an unknown mode.

Recommended input forms are:

```text
Worktree mode: full

Tasks:
1. <task using the shared full mode>
2. [read-only] <task overriding the shared mode>
```

## Dispatch workers

1. Parse the requested tasks and identify dependencies. Dispatch only independent
   tasks in parallel; sequence dependent tasks or ask for the missing ordering
   decision.
2. List current App projects and select the remote `multiomeR` project on host
   `esrum`. Verify its path is
   `/maps/projects/cbmr_shared/people/tqb695/non-GDPR/multiomeR`; do not rely on a
   remembered project ID.
3. Create one separate user-owned App task per item with environment type
   `worktree`. Omit `startingState` to start from the project's default branch
   unless the user explicitly names an existing branch or requests the current
   working tree.
4. Give each task a short unique title and plan a unique branch named
   `codex/<task-slug>`. Check for an existing branch and choose an unambiguous
   suffix unless the user requires an exact name.
5. Record the canonical worktree mode in each worker prompt.
6. Start independent task-creation calls concurrently.

Put the complete task scope and validation requirements in each worker's initial
prompt. Require the worker to leave unrelated files, the base checkout, and all
other worktrees untouched and to report progress and blockers in its own task.

## Bootstrap the selected mode

Require every worker to verify that `worktree` is not the base checkout, then
classify its initial environment before creating its branch or doing task work:

- `read-only`: `outputs` is a symlink.
- `full`: `outputs` is a directory containing `.worktree/setup.tsv`.
- `bare`: neither `.pixi` nor `outputs` exists.
- anything else: partial or unknown setup; stop and report a blocker.

If the existing mode matches the requested mode, skip setup and run only that
mode's validation. If the worktree is bare, run the matching setup and validate
it. If the existing mode differs, stop rather than converting or deleting it.
This preserves compatibility if the App starts applying its declared setup hook
automatically while retaining the explicit workaround for coordinator-created
bare worktrees.

### Read-only mode

```bash
worktree=$(git rev-parse --show-toplevel)
CODEX_WORKTREE_PATH="$worktree" dev/worktree_codex_read_only

test -L .pixi
test -d .pixi
test -L outputs
test -d outputs
test -f outputs/meta/meta
readlink -f .pixi
readlink -f outputs
```

Report both resolved symlink targets. Treat `outputs` as read-only and do not
perform a write probe inside it.

### Full mode

```bash
worktree=$(git rev-parse --show-toplevel)
CODEX_WORKTREE_PATH="$worktree" dev/worktree_codex_setup

test -L .pixi
test -d .pixi
test -d outputs
test ! -L outputs
test -w outputs
test -f outputs/meta/meta
test -f outputs/.worktree/setup.tsv
readlink -f .pixi
sed -n '1,8p' outputs/.worktree/setup.tsv
```

Report the resolved `.pixi` target plus the recorded scratch path. Treat this
task's `outputs` as writable and isolated from the base store. Do not share its
scratch store with another active worktree.

For either mode, stop and report a blocker if setup or any validation command
fails. Do not rerun full setup on an already configured full worktree, and do not
silently fall back from `full` to `read-only` or vice versa.

## Create the branch and perform the task

After bootstrap succeeds, require the worker to:

1. Create and switch to its planned `codex/<task-slug>` branch.
2. Confirm the current worktree and branch before editing.
3. Make only the requested changes.
4. Run validation proportional to the task and preserve the selected `outputs`
   contract.
5. Report the final branch, changed files, validation, and any blocker.

Do not commit, push, open a pull request, merge, or promote outputs unless the
user explicitly requests that action.

For a full worktree, retain `outputs/.worktree/setup.tsv` while work is active.
When the user explicitly requests teardown, run `dev/worktree_codex_delete`
before removing the worktree so its scratch marker is released.

Branch creation happens inside the worker because App task creation can select a
starting ref but cannot name a new branch.

## Return trackable tasks

Worktree creation may initially return only a client task ID. Wait for setup to
finish, find each permanent task ID, and set the intended short title. Return a
table with title, permanent task ID, planned branch, canonical worktree mode, and
current status. Emit one clickable task directive per successfully created task:

```text
::created-thread{threadId="<permanent-task-id>"}
```

Do not claim that a task is ready until it has a permanent ID. If provisioning
fails, report the client task ID and failure separately.

## Monitor and hand off

- For progress requests, read the worker tasks and summarize title, link, App
  status, branch, worktree mode, latest milestone, blocker, and next action
  without interrupting them.
- Send follow-up prompts to individual workers when the user asks for additional
  checks or changes.
- Use the App handoff operation only on worker task IDs, never on the calling
  coordinator task. Wait for handoff completion and report the resulting state.
- Hand workers back to the base checkout one at a time because the base checkout
  cannot host several worker branches simultaneously.

Keep the coordinator in the base checkout and free of worker edits.

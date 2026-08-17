---
name: multiomer-cancel-run
description: Cancel a currently running multiomeR targets run by finding its PID from the process file and sending SIGINT to stop it. Use when the user asks to stop, cancel, or kill a running multiomeR run, or when a run needs to be interrupted.
---

# multiomeR Cancel Run

The active root workflow records its orchestrator PID in the configured targets
store. Resolve the store from `_targets.yaml`; do not use the retained legacy
store under `pipelines/`.

## Steps

### 1. Resolve the store and PID

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
store <- targets::tar_config_get("store")
cat("store:", store, "\n")
cat("recorded pid:", targets::tar_pid(store = store), "\n")
EOF
```

The process record remains after a completed run, so a recorded PID is not proof
that the workflow is active.

### 2. Verify the process is running

```bash
ps -p <pid> -o pid,stat,etime,cmd
```

Confirm that the command is the expected R/targets orchestrator for this
checkout. If the PID is absent or belongs to another process, do not signal it.
If the user asked to clear a stale targets lock, use
`targets::tar_unblock_process()` only after this verification.

### 3. Send SIGINT

```bash
kill -INT <pid>
```

SIGINT requests graceful targets and crew shutdown. After signalling, poll the
orchestrator and report whether it exited; do not escalate to SIGKILL unless the
user explicitly requests forced termination.

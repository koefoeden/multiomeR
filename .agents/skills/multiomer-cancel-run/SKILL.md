---
name: multiomer-cancel-run
description: Cancel a currently running multiomeR targets run by finding its PID from the process file and sending SIGINT to stop it. Use when the user asks to stop, cancel, or kill a running multiomeR run, or when a run needs to be interrupted.
---

# multiomeR Cancel Run

The active root workflow stores its orchestrator PID in `<store>/meta/process`. Look up the root store path, then send SIGINT.

## Pipeline Store Paths

| Workflow | Store |
|---|---|
| active root `_targets.R` workflow | `pipelines/processing_and_aggregation/outputs` |

`differential_analyses` and `genetic_enrichment` are optional modules inside the same root workflow, so cancel their runs through this store. Extracted workflows such as cellranger, GWAS, and genotype processing live in sister repos and should be handled from those repos.

## Steps

### 1. Read the PID

```bash
awk -F'|' '$1=="pid"{print $2}' <store>/meta/process
```

### 2. Verify the process is running

```bash
ps -p <pid> -o pid,stat,etime,cmd
```

### 3. Send SIGINT

```bash
kill -2 <pid>
```

SIGINT (equivalent to Ctrl+C) stops the orchestrator and triggers graceful crew shutdown. If scheduler-backed workers are active, controller shutdown should cancel them.

## Quick One-Liner

```bash
kill -2 $(awk -F'|' '$1=="pid"{print $2}' <store>/meta/process)
```

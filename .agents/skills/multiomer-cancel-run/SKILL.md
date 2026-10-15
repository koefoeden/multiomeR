---
name: multiomer-cancel-run
description: Cancel a currently running multiomeR targets run by finding its driver from the process file and sending SIGINT to stop it. Use when the user asks to stop, cancel, or kill a running multiomeR run, or when a run needs to be interrupted.
---

# multiomeR Cancel Run

## 1. Find the driver

Run the active-run check in `multiomer-run-pipeline`. The recorded PID is the
process that runs the pipeline for the configured store.

- **Running here:** its working directory is this checkout. Confirm with
  `ps -p <pid> -o pid,etime,cmd` that it is an R process, then send SIGINT.
  Never pick a process by name, such as the first `callr` process on the node:
  other sessions run their own pipelines.
- **Running on another node:** the PID is not here, but progress was written
  within the last few minutes with targets dispatched. If the driver is its own
  scheduler job, cancel that job by its ID, such as `scancel <jobid>` on Slurm.
  Otherwise ask the session or user that started it to stop it.
- **Neither:** no run is active; report that and signal nothing.

## 2. Send SIGINT

```bash
kill -INT <pid>
```

SIGINT requests graceful targets and crew shutdown. After signalling, poll the
driver and report whether it exited; do not escalate to SIGKILL unless the user
explicitly requests forced termination. A few minutes later, check that the
run's worker jobs have left the queue, and cancel any that remain by their
explicit job IDs.

# Choose where the analysis runs

A **worker** is an R process that builds targets. A **controller** starts and manages workers, either on your machine or through a cluster scheduler. multiomeR defines its controllers with `crew` in `crew_controllers.R` in the repository root; the [targets distributed-computing guide](https://books.ropensci.org/targets/crew.html) explains the general setup.

Independent GEM wells, modalities, and analysis branches run in parallel, so more workers shorten a run until the longest chain of dependent targets limits it. Use local workers on a workstation that meets the [system requirements](demo_installation.md#system-requirements). For large datasets, use a cluster scheduler and ask your support team which account and resource limits to use.

## What to expect

The table shows two recorded runs to `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`, estimated from the recorded runtime of each target:

| Aggregation | GEM wells | Nuclei called by Cell Ranger ARC | Estimated critical path | Sum if targets ran one at a time |
|---|---:|---:|---:|---:|
| `immune_human_2x` | 2 | 17,277 | 17.3 minutes | 31.4 minutes |
| `PBMC_human_6x` | 6 | 51,291 | 23.8 minutes | 52.4 minutes |

`immune_human_2x` is the demo aggregation. `PBMC_human_6x` combines the six PBMC GEM wells in the public `cfg_GEM_wells.tsv`; the public `cfg_aggregations.yaml` does not define it.

The critical path is the longest chain of dependent targets: the shortest possible run time when enough workers are available. Fewer workers and scheduler queue time make real runs longer, as do more nuclei, peaks, plots, or modules. Treat the numbers as examples, not predictions for another machine or configuration, and time one representative aggregation before sizing a large run.

## Local execution

The committed `crew_controllers.R` is sized for a 16-CPU, 256-GB workstation: four `local-light` workers (1 core, 16 GB each) and two `local-heavy` workers (6 cores, 60 GB each). On a machine near the 60-GB minimum, lower the two `workers` values so that only one heavy target runs at a time:

```{.r filename="crew_controllers.R"}
controller_list <- list(
  crew::crew_controller_local(
    name = "local-light",
    workers = 2
  ),
  crew::crew_controller_local(
    name = "local-heavy",
    workers = 1,
    crashes_max = 1
  )
)
```

The `RAM_GB` values in `controller_resources_tibble` route each target to a controller with enough declared memory; they do not limit memory use. Choose `workers` values so that the targets that can run at once fit in physical memory, and raise them only when you know the memory headroom.

After editing the file, restart R or reload the runtime, which also checks the file:

```{.r filename="R"}
load_project_runtime()
```

## Scheduler execution

For SLURM, PBS, SGE, or LSF, replace the local controllers with the matching `crew.cluster` controllers. The commented SLURM example in `crew_controllers.R` defines light, heavy, and GPU tiers. For each tier:

1. Request, in the controller's scheduler options, the CPUs, memory, and GPUs that its `controller_resources_tibble` row declares.
2. Set the queue, account, wall time, modules, and worker start-up commands your cluster requires.
3. Give GPUs their own tier: only targets that request GPUs are routed to it.

Test a small target selection before increasing worker counts. For start-up and routing errors, see [Troubleshooting](troubleshooting.md#controller-and-scheduler-failures).

## Rules for `crew_controllers.R` {#controller-rules}

The last expression in the file must be a list with `controller_list` and `controller_resources_tibble`:

- Each `controller_name` in the table names one controller in `controller_list`. Names are unique.
- The table has exactly the columns `controller_name`, `cores`, `RAM_GB`, and `gpus`, in that order, with numeric, non-missing resources.
- Each target runs on the smallest tier that meets its CPU, RAM, and GPU request; targets without a request run on the smallest tier. Tiers are compared by GPUs, then cores, then RAM, so row order does not matter.

[Implementation conventions](implementation/implementation_conventions.html#runtime-bootstrap) describe how the runtime loads this file and how targets request resources.

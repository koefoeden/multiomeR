# Install and prepare the demo

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

## System requirements

- Linux system with `git` and `curl`
- At least 60 GB of RAM. This is enough for one heavy target at a time; machines near the minimum should reduce concurrent workers in `crew_controllers.R`.
- Multiple CPU cores are strongly recommended. The timing quoted in the next chapter was measured with 16 logical threads.

::: {.callout-tip title="Machines with less than 256 GB of RAM"}
The committed `crew_controllers.R` is sized for a 16-CPU, 256-GB workstation: four light workers and two heavy workers that may each use 60 GB. On a machine near the 60-GB minimum, edit the two `workers` values in `crew_controllers.R` before running the demo so that only one heavy target runs at a time:

``` {.r filename="crew_controllers.R"}
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

The `RAM_GB` values in the same file describe routing capacity, not enforced limits, so the workers that can run at once must fit in physical memory. To run on a SLURM or other scheduler instead, see [Choose where the analysis runs](performance_distributed_computing.md).
:::

## Set up the demo

Run this block from the directory where you want to clone multiomeR to. The single `pixi run` setup command installs the locked environment before its `setup-demo` task downloads the two configured public inputs and installs the pinned GitHub-only R packages.

``` {.bash filename="Bash"}
# Clone the repository and enter its root directory.
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR

# Skip these two lines when pixi is already available on PATH.
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"

# Install the locked environment, download 3.9 GB of demo inputs, and install
# the pinned GitHub versions of BPCells, Signac, and betterChromVAR.
pixi run --use-environment-activation-cache --locked --run-post-link-scripts setup-demo

# Start R in the configured environment for the commands in the next chapter.
pixi run --use-environment-activation-cache --locked R
```

The download task is restart-safe: non-empty files already present under `example_data` are skipped.

Continue to [Run the demo](demo_running.md) from the R prompt.
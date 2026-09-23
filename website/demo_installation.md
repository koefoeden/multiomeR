# Install and prepare the demo

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

The demo analyzes two public 10x Genomics GEM wells: peripheral blood mononuclear cells from a healthy donor and a lymph node with lymphoma, 17,277 nuclei in total.

## System requirements

- 64-bit x86 Linux with `git` and `curl`.
- At least 60 GB of RAM, enough for one memory-intensive target at a time.
- About 30 GB of free disk space. This covers 3.9 GB of demo inputs, about 6 GB of results and up to 15 GB for the Pixi environment and its package cache, which Pixi keeps in your home directory by default.
- Multiple CPU cores. The run time in the next chapter was measured with 16 logical threads.

The committed `crew_controllers.R` suits a 16-CPU, 256-GB workstation. On a machine near the 60-GB minimum, lower its worker counts before running the demo, starting with the heavy workers; see [Local execution](performance_distributed_computing.md#local-execution).

## Set up the demo

Run these commands in the folder where the clone should be created.

```{.bash filename="Bash"}
# Clone the repository and enter its root folder.
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR

# Install Pixi. Skip these two lines if `pixi` is already on your PATH.
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"

# Install the locked environment, download the demo inputs into example_data/,
# and install the pinned GitHub versions of BPCells, Signac and betterChromVAR.
pixi run --use-environment-activation-cache --locked --run-post-link-scripts setup-demo

# Start R in that environment.
pixi run --use-environment-activation-cache --locked R
```

If the download is interrupted, run the `setup-demo` command again; it skips files that are already downloaded.

Keep the R session open and continue to [Run the demo](demo_running.md).
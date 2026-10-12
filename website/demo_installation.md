# Install and prepare the demo

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

The demo analyzes two public 10x Genomics GEM wells of peripheral blood mononuclear cells (PBMCs) from two healthy donors: one with granulocytes removed by cell sorting (2,711 nuclei) and 10x Genomics' 10k PBMC set run on the Chromium Controller (10,691 nuclei), 13,402 nuclei in total.

## System requirements

- 64-bit x86 Linux with `git` and `curl`.
- At least 32 GB of RAM. The most memory-intensive targets need about 7 GB plus 2.5 GB per 100,000 nuclei: about 12 GB for 213,000 nuclei and 25 GB for 713,000. The optional `Seurat_SCT` GEX backend needs more, about 0.25 GB per 1,000 nuclei to normalize gene expression.
- About 30 GB of free disk space. This covers 4.6 GB of demo inputs, about 6.6 GB of results and up to 15 GB for the Pixi environment and its package cache, which Pixi keeps in your home directory by default.
- Multiple CPU cores. The run time in the next chapter was measured with 16 logical threads.

The committed `configuration/crew_controllers.R` runs the demo on a 16-CPU, 32-GB machine. For larger aggregations on such a machine, lower its worker counts, starting with the heavy workers; see [Local execution](performance_distributed_computing.md#local-execution).

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
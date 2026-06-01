#!/bin/bash

for version in 4.4.2 4.5.1; do
    module load R/$version --auto
    Rscript - <<'EOF'
    renv::init(bioconductor = TRUE, bare = TRUE)
EOF
# we need to restart R here (maybe?) to ensure repos are set correctly.
    Rscript - <<'EOF'
    Sys.setenv("RENV_CONFIG_PAK_ENABLED" = "TRUE")
    options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
    renv::install("tidyverse")
    source("R/shared/dev_helpers.R")
    snapshot_packages_for_R_version()
EOF
    module unload R/$version
done
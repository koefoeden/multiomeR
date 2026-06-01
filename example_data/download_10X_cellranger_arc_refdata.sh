#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

wget -c https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
tar -xzf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

wget -c https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -xzf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

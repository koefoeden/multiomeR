#!/bin/bash
set -euo pipefail

# Script to download public cellranger-arc output datasets from the official 10x site to serve as test data in the pipeline
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

# Embryonic mouse data from https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0
mkdir -p embryonic_brain_mouse/outs
wget -O embryonic_brain_mouse/outs/web_summary.html https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_web_summary.html
wget -O embryonic_brain_mouse/outs/summary.csv https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_summary.csv
wget -O embryonic_brain_mouse/outs/per_barcode_metrics.csv https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_per_barcode_metrics.csv
wget -O embryonic_brain_mouse/outs/filtered_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5
wget -O embryonic_brain_mouse/outs/raw_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_raw_feature_bc_matrix.h5
wget -O embryonic_brain_mouse/outs/atac_fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz
wget -O embryonic_brain_mouse/outs/atac_fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz.tbi

# PBMC from healthy human data from https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0
mkdir -p healthy_PBMC_human/outs
wget -O healthy_PBMC_human/outs/web_summary.html https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_web_summary.html
wget -O healthy_PBMC_human/outs/summary.csv https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_summary.csv
wget -O healthy_PBMC_human/outs/per_barcode_metrics.csv https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv
wget -O healthy_PBMC_human/outs/filtered_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5
wget -O healthy_PBMC_human/outs/raw_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5
wget -O healthy_PBMC_human/outs/atac_fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz
wget -O healthy_PBMC_human/outs/atac_fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi


# Flash-Frozen Lymph Node with B Cell Lymphoma from https://www.10xgenomics.com/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-2-0-0
# Output Files
mkdir -p lymphoma_lymph_human/outs
wget -O lymphoma_lymph_human/outs/web_summary.html https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_web_summary.html
wget -O lymphoma_lymph_human/outs/summary.csv https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_summary.csv
wget -O lymphoma_lymph_human/outs/per_barcode_metrics.csv https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_per_barcode_metrics.csv
wget -O lymphoma_lymph_human/outs/filtered_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5
wget -O lymphoma_lymph_human/outs/raw_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_raw_feature_bc_matrix.h5
wget -O lymphoma_lymph_human/outs/atac_fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_atac_fragments.tsv.gz
wget -O lymphoma_lymph_human/outs/atac_fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_atac_fragments.tsv.gz.tbi

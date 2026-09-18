# Plan your analysis

Before configuring anything, please check that you have the following inputs available:

**cellranger-arc count directories:** The output directories produced by `cellranger-arc count` containing the `outs/`-folder with `summary.csv`, `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`, its `.tbi` index, and `per_barcode_metrics.csv` inside. All count directories must have been generated using the same reference if you want to combine later in the pipeline.

**donor metadata TSV** with one row per donor and the phenotypes or covariates you will use; see the [donor metadata table](reference_donor_metadata.md).

**VCF files for demultiplexing by genotype (optional):** If you have multiplexed several donors on one or more GEM-wells, and you wish to demultiplex them using reference genotypes, please prepare a VCF file for each unique combination/batch of donors as described in the vireo documentation \<format as link to the right documentation\> . This functionality also requires `atac_possorted_bam.bam` in the cellragner-arc count directories.

**CellBender H5 files (optional)**: If you wish the pipeline to use gene-expression data that has been filtered for ambient RNA, please run cellbender \<format as link to right place\> beforehand.

**Compute setup:** If your dataset is large, it is highly recommended to run the pipeline on a compute cluster with a job-scheduler available. See [Choose where the analysis runs](performance_distributed_computing.md) for more info.

When these things are in order, please continue to [Run your own analysis](main_running.md#steps).
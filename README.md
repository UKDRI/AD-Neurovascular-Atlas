# Project overview

This project is using single-nuclear RNA sequencing data to investigate the blood brain barrier in Alzheimer's disease.

It employs a novel method to physically separate the vasculature from the parenchymal tissue from human post-mortem prefrontal cortex samples.

There are 20 cases and 20 controls each with a parenchyma and vascular fraction sequenced for a total of 80 samples via 10X.

Data can be found on GEO [here]()

## Processing of data

FastQC/MultiQC were run to assess sequencing quality.
Subsequently standard processing was done with [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), using an updated reference with Ensembl version 109 and Gencode version 43.
The output of this was run through [scFlow](https://github.com/combiz/nf-core-scflow/tree/dev-nf) with standard parameters which was then read into R via the [Seurat](https://satijalab.org/seurat/) package.

Looking at the read counts and number of features output from scFlow revealed that 3 of the donors, and the vascular fraction from a fourth donor, were of poor quality and so they were excluded from downstream analysis.
PCA was performed and 35 dimensions were used to find nearest-neighbors and a resolution of 0.6 was used to find clusters via Seurat.
UMAP was computed with the same 35 dimensions.
Cluster markers were found using genes with at least 20% of cells in a cluster expressing them, and a minimum log foldchange of 0.5, again via Seurat.

Differences in celltype proportions between cases and controls were examined with the [speckle](https://github.com/phipsonlab/speckle) R package.
It leverages biological replication to obtain measures of variability of cell type proportion estimates and uses empirical Bayes to stabilize variance estimates by borrowing information across cell types.

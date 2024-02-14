# Project overview

This project is using single-nuclear RNA sequencing data to investigate the blood brain barrier in Alzheimer's disease.

There are 20 cases and 20 controls with two cell populations, from the parenchyma and vascular, sequenced for a total of 80 samples.

## Data

The are three sets of data.

Set 1 contains the first 12 donors, 24 samples (12x P - parenchymal and 12x V - vascular fraction), which are independently sequenced of Set_2 (the remaining 16 samples from 8 donors).
The sequencing was performed for set 1 initially in Oxford, but there were some issues with index hopping and cleaning the data did not seem to work too well.
So Set 1 was re-sequenced in Cardiff and the folder is: 211008_A00748_0157_AHT7TJDSX2_fastq_L2_3_4

For set 2 the data was sequenced initially on one Illumina lane, but then more lanes were added for better coverage.

The last set is another bath of 40 samples from 20 donors.

## Processing of data

FastQC/MultiQC were run to assess sequencing quality, no issues were found with the latest batch of samples.
Subsequently standard processing was done with [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), using an updated reference with Ensembl version 109 and Gencode version 43.
The output of this was run through [scFlow](https://github.com/combiz/nf-core-scflow/tree/dev-nf) with standard parameters which was then read into R via the [Seurat](https://satijalab.org/seurat/) package.

Looking at the read counts and number of features output from scFlow revealed that 3 of the donors, and the vascular fraction from a fourth donor, were of poor quality and so they were excluded from downstream analysis. 
PCA was performed and 35 dimensions were used to find nearest-neighbors and a resolution of 0.6 was used to find clusters via Seurat.
UMAP was computed with the same 35 dimensions.
Cluster markers were found using genes with at least 20% of cells in a cluster expressing them, and a minimum log foldchange of 0.5, again via Seurat.

Differences in celltype proportions between cases and controls were examined with the [speckle](https://github.com/phipsonlab/speckle) R package.
It leverages biological replication to obtain measures of variability of cell type proportion estimates and uses empirical Bayes to stabilize variance estimates by borrowing information across cell types. 

## Contributors

- Frank Wessely
- Mateus Harrington
- Caleb Webber

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

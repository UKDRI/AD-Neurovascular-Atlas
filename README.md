# Neurovascular Atlas of Alzheimer's Disease

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI](https://zenodo.org/badge/1105323397.svg)](https://doi.org/10.5281/zenodo.17737331)

## Overview

This repository contains the analysis code for a single-nucleus RNA sequencing study of the neurovascular unit (NVU) in Alzheimer's disease (AD).
We profiled prefrontal cortex samples from AD patients and age-matched controls to generate a high-resolution atlas of vascular and parenchymal cell types.

**Associated manuscript:** *In preparation*

## Processing of data

In brief: standard processing was done with [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest), using an updated reference with Ensembl version 109 and Gencode version 43.
The output of this was run through [scFlow](https://github.com/combiz/nf-core-scflow/tree/dev-nf) with standard parameters which was then read into R via the [Seurat](https://satijalab.org/seurat/) package.

Looking at the read counts and number of features output from scFlow revealed that 3 of the donors, and the vascular fraction from a fourth donor, were of poor quality and so they were excluded from downstream analysis.

## Data Availability

The raw sequencing data (10X Chromium single-nucleus RNA-seq) and processed count matrices are available from GEO, accessions GSE310554 & GSE222007.
Note that GSE310554 contains the processed data of all samples from both datasets except for those donors mentioned above, and is a Seurat object.

## Computational environment

There's a few options if you want to replicate the R environment used for this analysis.

### Option 1: Using Docker (Recommended)

The easiest way to reproduce our analysis environment is using Docker.

#### Pull from Docker Hub

The image is on Docker Hub [here](https://hub.docker.com/r/multitude5286/ad-bbb-analysis)

```bash
# Pull the image from Docker Hub
docker pull multitude5286/ad-bbb-analysis
```

#### Run the Docker container

Having got the Docker image ready you can either enter it as an interactive R session

```bash
# Interactive R session
docker run --rm -ti ad-bbb-analysis R
```

Or use it via [RStudio](https://posit.co/downloads/) like so

```bash
# Use RStudio Server (optional)
docker run --rm -ti \
  -e USER=rstudio \
  -e PASSWORD=yourpassword \
  -p 8787:8787 \
  ad-bbb-analysis \
  /init
```

Having run this you can visit [http://localhost:8787/](http://localhost:8787/) and use username: `rstudio`, password: `yourpassword`

#### Build the Docker image yourself

Alternatively you can build it yourself with the Dockerfile in this repo

```bash
# Clone this repository
git clone https://github.com/UKDRI/AD-Neurovascular-Atlas.git
cd ad-neurovascular-atlas

# Build the Docker image (takes ~30-60 minutes)
bash build-docker.sh

# Or build manually
docker build -t ad-bbb-analysis:1.0.0 .
```

### Option 2: Using renv (Local Installation)

If you prefer to install packages locally:

```bash
# Clone repository
git clone https://github.com/UKDRI/AD-Neurovascular-Atlas.git
cd ad-neurovascular-atlas

# Install renv if you don't have it
R -e "install.packages('renv')"

# Restore the R environment (this will take a while)
R -e "renv::restore()"
```

Note: You'll need to install system dependencies manually.
See Dockerfile for the required system libraries.

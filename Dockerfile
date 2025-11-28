FROM rocker/rstudio:4.5.1

# Build arguments
ARG VERSION="1.0.0"
ARG BUILD_DATE
ARG GIT_COMMIT
ARG GITHUB_REPO
ARG ZENODO_DOI

# Metadata labels
LABEL maintainer="Bernardo-harringtong@cardiff.ac.uk"
LABEL org.opencontainers.image.title="AD-BBB Analysis Environment"
LABEL org.opencontainers.image.description="R environment for AD blood-brain barrier single-cell RNA-seq analysis"
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.created="${BUILD_DATE}"
LABEL org.opencontainers.image.revision="${GIT_COMMIT}"
LABEL org.opencontainers.image.url="https://github.com/${GITHUB_REPO}"
LABEL org.opencontainers.image.source="https://github.com/${GITHUB_REPO}"
LABEL org.opencontainers.image.documentation="https://doi.org/${ZENODO_DOI}"
LABEL org.opencontainers.image.licenses="CC-BY"
LABEL org.opencontainers.image.base.name="rocker/rstudio:4.5.1"
LABEL r.version="4.5.1"
LABEL zenodo.doi="${ZENODO_DOI}"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    make libcurl4-openssl-dev libnode-dev libxml2-dev libx11-dev \
    libssl-dev libcairo2-dev libfontconfig1-dev libfreetype6-dev \
    git zlib1g-dev libglpk-dev pandoc libmagick++-dev gsfonts cmake \
    libpng-dev libjpeg-dev libtiff-dev python3 libabsl-dev libgdal-dev \
    gdal-bin libgeos-dev libproj-dev libsqlite3-dev libicu-dev \
    libfribidi-dev libharfbuzz-dev libudunits2-dev default-jdk \
    libfftw3-dev libgit2-dev libhdf5-dev libhwloc-dev libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Use Posit Public Package Manager for faster binary installation
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' >> /usr/local/lib/R/etc/Rprofile.site

WORKDIR /home/ad-bbb

# Install CRAN packages using install2.r (parallel installation)
RUN install2.r --error --skipinstalled --ncpus -1 \
    stringr \
    DT \
    reshape \
    dplyr \
    tidyr \
    tibble \
    purrr \
    here \
    readr \
    kableExtra \
    janitor \
    labelled \
    BiocManager \
    data.table \
    gt \
    europepmc \
    readxl \
    furrr \
    future.apply \
    parallel \
    targets \
    forcats \
    msigdbr \
    downlit \
    xml2 \
    corrgram \
    RColorBrewer \
    pheatmap \
    ggplot2 \
    ggsci \
    cowplot \
    scales \
    ggvenn \
    ggpubr \
    gridExtra \
    tidytext \
    viridis \
    wesanderson \
    pals \
    ggalluvial \
    GGally \
    ggrastr \
    ggcorrplot \
    circlize \
    ggpointdensity \
    && rm -rf /tmp/downloaded_packages

RUN R -e "options(timeout = 600); BiocManager::install(c( \
    'Biobase' \
), ask = FALSE, update = FALSE)" \
    && rm -rf /tmp/downloaded_packages

RUN install2.r --error --skipinstalled --ncpus -1 \
    NMF \
    && rm -rf /tmp/downloaded_packages

# Install remotes for GitHub packages
RUN install2.r --error --skipinstalled remotes && rm -rf /tmp/downloaded_packages

# Install Bioconductor packages (with explicit timeout in the R command as backup)
RUN R -e "options(timeout = 600); BiocManager::install(c( \
    'cellity', \
    'scater', \
    'Seurat', \
    'SC3', \
    'DESeq2', \
    'edgeR', \
    'BiocParallel', \
    'org.Hs.eg.db', \
    'pathview', \
    'gage', \
    'gageData', \
    'DropletUtils', \
    'ensembldb', \
    'AnnotationHub', \
    'patchwork', \
    'scran', \
    'scuttle', \
    'PCAtools', \
    'batchelor', \
    'bluster', \
    'igraph', \
    'SingleCellExperiment', \
    'scRNAseq', \
    'biomaRt', \
    'GEOquery', \
    'AnnotationDbi', \
    'clusterProfiler', \
    'DOSE', \
    'enrichplot', \
    'EnsDb.Hsapiens.v79', \
    'fgsea', \
    'MAST', \
    'speckle', \
    'rrvgo', \
    'MungeSumstats', \
    'BSgenome.Hsapiens.1000genomes.hs37d5', \
    'BSgenome.Hsapiens.NCBI.GRCh38', \
    'SNPlocs.Hsapiens.dbSNP155.GRCh37', \
    'qusage', \
    'SNPlocs.Hsapiens.dbSNP155.GRCh38', \
    'GenomicFiles', \
    'simona', \
    'ComplexHeatmap', \
    'switchde', \
    'slingshot', \
    'tradeSeq' \
), ask = FALSE, update = FALSE)" \
    && rm -rf /tmp/downloaded_packages

# Install GitHub packages with timeout
RUN R -e "remotes::install_github('jinworks/CellChat', upgrade = 'never')"
RUN R -e "remotes::install_github('rpolicastro/scProportionTest', upgrade = 'never')"
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade = 'never')"
RUN R -e "remotes::install_github('neurogenomics/scFlow', upgrade = 'never')"
RUN R -e "remotes::install_github('immunogenomics/presto', upgrade = 'never')"
RUN R -e "remotes::install_github('cole-trapnell-lab/monocle3', upgrade = 'never')"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers', upgrade = 'never')"

# Clean up
RUN rm -rf /tmp/downloaded_packages /tmp/*.rds

# Set working directory for your project
WORKDIR /home/ad-bbb


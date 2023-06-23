# load packages
if (!require("Seurat", character.only = TRUE)) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("future", character.only = TRUE)) {
  install.packages("future")
  library(future)
}

# read data
seurat <- readr::read_rds(here::here("03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/scflow-seurat-preprocessing.rds"))

## number of cores
n_workers <- 30
## 200 GB
options(future.globals.maxSize = 200000 * 1024^2)
plan(strategy = "multicore", workers = n_workers)

## normalisation
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

## variable genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

## scaling
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

## PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

# Save
readr::write_rds(seurat, here::here("03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/scflow-seurat-postprocessing.rds"))

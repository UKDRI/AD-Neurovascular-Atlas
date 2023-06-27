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
seurat <- readr::read_rds(here::here("03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/scflow-seurat-preprocessing.rds/scflow-seurat-preprocessing.rds"))

## number of cores
#n_workers <- 2
### 365 GB
#options(future.globals.maxSize = 369000 * 1024^2)
#plan(strategy = "multicore", workers = n_workers)

## normalisation
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

## variable genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

## scaling
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

## PCA
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

print("Finished running PCA")

# save pca plots
pca_plot <- DimPlot(seurat, reduction = "pca", group.by = "diagnosis")
elbow_plot <- ElbowPlot(seurat, ndims = 50)
# save plots
ggplot2::ggsave(
  here::here(
    "03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/pca_plot.png"
  ),
  plot = pca_plot
)
ggplot2::ggsave(
  here::here(
    "03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/elbow_plot.png"
  ),
  plot = elbow_plot
)

# remove plots to save memory
rm(pca_plot, elbow_plot)

# find clusters/neighbors
n_dims <- 20
## n_dims <- 30
## resolution_value <- 0.5
resolution_value <- 0.3

cat("Dimensions of reduction used: ", n_dims, fill=T)
cat("Resolution: ", resolution_value, fill=T)

seurat <- FindNeighbors(seurat, dims = 1:n_dims)
seurat <- FindClusters(seurat, resolution = resolution_value)

rm(resolution_value)

print("Running UMAP")

seurat <- RunUMAP(seurat, dims = 1:n_dims)

# The scale data slot is huge, too much to read in on my laptop
# It's only needed to compute dimensional reductions (PCA/UMAP) so 
# we can remove it now
print("Remove scaled data")

seurat <- DietSeurat(seurat)

print("Saving data")

# Save
readr::write_rds(seurat, here::here("03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/scflow-seurat-postprocessing.rds"))

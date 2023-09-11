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
n_workers <- 2
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
#n_dims <- 20
n_dims <- 35
## resolution_value <- 0.3
resolution_value <- 0.6

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

#seurat@assays$RNA@scale.data <- NULL
# keep reductions and graphs
reductions <- names(seurat@reductions)
graphs <- names(seurat@graphs)

seurat <- DietSeurat(seurat, dimreducs = reductions, graphs = graphs)

rm(reductions, graphs)

# Since the find markers takes a while I'll compute this here too
print("Find markers")

## options(future.globals.maxSize = 4000 * 1024^2) ## 4 GB
options(future.globals.maxSize = 360000 * 1024^2) ## 360 GB
plan(strategy = "multicore", workers = n_workers)

cluster_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

print("Saving data")

# Save
readr::write_rds(seurat, here::here("03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/scflow-seurat-postprocessing_35-dims.rds"))
readr::write_rds(cluster_markers, here::here("03_data/990_processed_data/001_snrnaseq/07_scflow_analysis/cluter_markers_35-dims.rds"))

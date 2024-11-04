# load packages
if (!require("Seurat", character.only = TRUE)) {
  install.packages("Seurat")
  library(Seurat)
}
if (!require("future", character.only = TRUE)) {
  install.packages("future")
  library(future)
}
if (!require("qs", character.only = TRUE)) {
  install.packages("qs")
  library(qs)
}
library(slingshot)
library(tradeSeq)

# read data
sce_slingshot <- qs::qread(here::here("03_data/990_processed_data/008_pseudotime",
                              "slingshot_obj.qs"))

# subsample to a very small number of cells to test code is working
subsample_cells <- function(sce, celltype_col, n = 5000) {
  cell_types <- unique(colData(sce)[[celltype_col]])
  sampled_cells <- unlist(lapply(cell_types, function(ct) {
    cells <- which(colData(sce)[[celltype_col]] == ct)
    if (length(cells) > n) {
      sample(cells, n)
    } else {
      cells
    }
  }))
  return(sampled_cells)
}

celltype_col <- "level1_celltypes_with_endomt_subclusters"
set.seed(1234)
sampled_cells <- subsample_cells(sce_slingshot, celltype_col)
sce_subsampled <- sce_slingshot[, sampled_cells]

sce_subsampled <- slingshot(
  sce_subsampled,
  clusterLabels = 'level1_celltypes_with_endomt_subclusters',
  start.clus = "Endothelial",
  reducedDim = 'PCA'
)

sds <- SlingshotDataSet(sce_subsampled)

# Annoyingly this doesn't subset all the elements of the sds object, 
# so I'll need to do that manually
original_curves <- slingCurves(sds)

# Create a subset version of the curves
subset_curves <- lapply(original_curves, function(curve) {
  curve$lambda <- curve$lambda[sampled_cells]
  curve$dist_ind <- curve$dist_ind[sampled_cells]
  curve$w<- curve$w[sampled_cells]
  
  return(curve)
})

sds@curves <- subset_curves

# Set up parallel processing
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8
BPPARAM

# Evaluate the optimal number of knots
file <- here::here("03_data/990_processed_data/008_pseudotime/slingshot_kvalues.qs")
if(!file.exists(file)) {
  
set.seed(123)  # For reproducibility
print("Running evaluateK...")
k_values <- evaluateK(counts = assays(sce_subsampled)$counts, 
                      sds = sds, nGenes = 500, k = 3:10, 
                      parallel = TRUE, BPPARAM = BPPARAM)

qs::qsave(k_values, file)
print("Finished evaluateK!")
} else {
  k_values <- qs::qread(file)
}


#optimal_k <- which.min(k_values$ic$BIC)
# Fit GAMs for each gene along the pseudotime
print("Running fitGAM...")
sce_subsampled <- fitGAM(
  counts = assays(sce_subsampled)$counts,
  pseudotime = slingPseudotime(sds, na = FALSE),
  cellWeights = slingCurveWeights(sds),
  nknots = 6,
  parallel = TRUE,
  BPPARAM = BPPARAM
)
print("Finished fitGAM!")

# Extract the model results
gam_results <- rowData(sce_slingshot)$tradeSeq

# Save
qs::qsave(sce_slingshot, here::here("03_data/990_processed_data/008_pseudotime/slingshot_tradeseq.qs"))

print("gam results head:")
print(head(gam_results))


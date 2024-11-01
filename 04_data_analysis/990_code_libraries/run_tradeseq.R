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


sds <- SlingshotDataSet(sce_subsampled)
# Set up parallel processing
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- parallel::detectCores()
BPPARAM

# Evaluate the optimal number of knots
set.seed(123)  # For reproducibility
print("Running evaluateK...")
k_values <- evaluateK(counts = assays(sce_subsampled)$counts, 
                      sds = test$slingshot, nGenes = 500, k = 3:10, 
                      parallel = TRUE, BPPARAM = BPPARAM)
print("Finished evaluateK!")

qs::qsave(k_values, here::here("03_data/990_processed_data/008_pseudotime/slingshot_kvalues.qs"))

k_values
optimal_k <- which.min(k_values$ic$BIC)
# Fit GAMs for each gene along the pseudotime
print("Running fitGAM...")
sce_subsampled <- fitGAM(
  counts = assays(sce_subsample)$counts,
  pseudotime = slingPseudotime(sds, na = FALSE),
  cellWeights = slingCurveWeights(sds),
  nknots = optimal_k,
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


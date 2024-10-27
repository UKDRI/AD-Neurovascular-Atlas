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
sce_slingshot <- qs::qread(sce_slingshot, here("03_data/990_processed_data/008_pseudotime",
                              "slingshot_obj.qs"))

sds <- SlingshotDataSet(sce_slingshot)
# Set up parallel processing
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- parallel::detectCores()
BPPARAM

# Evaluate the optimal number of knots
set.seed(123)  # For reproducibility
k_values <- evaluateK(counts = assays(sce_slingshot)$counts, sds = sds, 
                      nGenes = 500, k = 3:10, parallel = TRUE, BPPARAM = BPPARAM)
k_values

optimal_k <- which.min(k_values$ic$BIC)
# Having a hard time running evaluateK, it's slow and seems to crash, I'll just try progressing with nknots of 6 (the default)
# Fit GAMs for each gene along the pseudotime
sce_slingshot <- fitGAM(counts = assays(sce_slingshot)$counts, pseudotime = slingPseudotime(sds, na = FALSE), cellWeights = slingCurveWeights(sds), nknots = optimal_k, parallel = TRUE, BPPARAM = BPPARAM)

# Extract the model results
gam_results <- rowData(sce_slingshot)$tradeSeq

print("gam results head:")
print(head(gam_results))

# Save
qs::qsave(sce_slingshot, here::here("03_data/990_processed_data/008_pseudotime/slingshot_tradeseq.qs"))
qs::qsave(k_values, here::here("03_data/990_processed_data/008_pseudotime/slingshot_kvalues.qs"))

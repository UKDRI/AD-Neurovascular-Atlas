# load packages
if (!require("qs", character.only = TRUE)) {
  install.packages("qs")
  library(qs)
}
if (!require("tidyr", character.only = TRUE)) {
  install.packages("tidyr")
  library(tidyr)
}
if (!require("dplyr", character.only = TRUE)) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require("stringr", character.only = TRUE)) {
  install.packages("stringr")
  library(stringr)
}
if (!require("scran", character.only = TRUE)) {
  install.packages("scran")
  library(scran)
}
library(slingshot)
library(tradeSeq)

# read data 
# NOTE: for some reason the data was behaving differently on the HPC than 
# locally, some data was lost when the sds object is created, so I've run this 
# part locally and just read in the file on the HPC

# set subset
subset = "pericyte"

if(subset == "pericyte"){
  sce_slingshot <- qs::qload(here::here("03_data/990_processed_data/008_pseudotime",
                              "slingshot_obj_pericyte.qs"))
} else {
  sce_slingshot <- qs::qload(here::here("03_data/990_processed_data/008_pseudotime",
                              "slingshot_obj_smc.qs"))
}

# Function to create balanced subsample for each celltype - vectorized version
create_balanced_subsample <- function(sce,
                                      celltype_col = "celltype",
                                      condition_col = "diagnosis",
                                      cells_per_group = 1500,
                                      seed = 42) {
  
  set.seed(seed)
  
  # Convert colData to tibble for easier manipulation
  cell_data <- colData(sce)[,c(celltype_col, condition_col)] |> 
    as_tibble() |>
    dplyr::mutate(cell_id = colnames(sce))
    
  # Create balanced sample using dplyr operations
  if(cases_or_controls == "controls") {
      selected_cells <- cell_data |>
        dplyr::filter(base::get(condition_col) == "Control") |>
        # Group by both celltype and condition
        group_by(across(all_of(c(celltype_col, condition_col)))) |>
        # Sample cells from each group
        slice_sample(n = cells_per_group * 2) |>
        ungroup()
  } else if (cases_or_controls == "cases") {
      selected_cells <- cell_data |>
        dplyr::filter(base::get(condition_col) == "Case") |>
        # Group by both celltype and condition
        group_by(across(all_of(c(celltype_col, condition_col)))) |>
        # Sample cells from each group
        slice_sample(n = cells_per_group * 2) |>
        ungroup()
  } else {
      # Create balanced sample using dplyr operations
      selected_cells <- cell_data |>
        # Group by both celltype and condition
        group_by(across(all_of(c(celltype_col, condition_col)))) |>
        # Sample cells from each group
        slice_sample(n = cells_per_group) |>
        ungroup()
  }
  
  
  # Print summary of balanced dataset
  summary_table <- selected_cells |>
    dplyr::count(across(all_of(c(celltype_col, condition_col)))) |>
    pivot_wider(names_from = all_of(condition_col),
                values_from = n)
  
  cat("Summary of balanced dataset for ", cases_or_controls, ":", sep = "")
  print(summary_table)
  
  return(selected_cells)
}

celltype_col <- "level1_celltypes_with_endomt_subclusters"
selected_cells <- create_balanced_subsample(sce_slingshot, celltype_col)
selected_cells <- selected_cells$cell_id
# Subset original SCE
sce_subsampled <- sce_slingshot[, selected_cells]

#--------- Fiilter out genes that're lowly expressed
# Identify cell types in your data
cell_types <- unique(sce_subsampled$level1_celltypes_with_endomt_subclusters)

# Create a list to store genes per cell type
genes_per_celltype <- lapply(cell_types, function(ct) {
  # Subset to specific cell type
  ct_sce <- sce_subsampled[, sce_subsampled$level1_celltypes_with_endomt_subclusters == ct]
  
  # Calculate gene expression proportion in this cell type
  ct_counts <- counts(ct_sce)
  ct_gene_prop <- rowMeans(ct_counts > 0)
  
  # Keep genes expressed in at least 10% of cells in this cell type
  names(ct_gene_prop)[ct_gene_prop >= 0.1]
})

# Find genes that are expressed in at least one cell type
genes_to_keep <- unique(unlist(genes_per_celltype))

# Filter SCE object
sce_subsampled <- sce_subsampled[genes_to_keep, ]

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
  curve$lambda <- curve$lambda[selected_cells]
  curve$dist_ind <- curve$dist_ind[selected_cells]
  curve$w <- curve$w[selected_cells]
  
  return(curve)
})


# Set up parallel processing
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 12
BPPARAM

# It gave an error due to NAs in pseudotime, which should be impossible, 
# but perhaps due to aspect not being subset equally
pseudotime <- slingPseudotime(sds, na = FALSE)
na_counts <- sum(is.na(pseudotime))
if (na_counts > 0) {
  stop("There are still NA values in the pseudotime matrix.")
}

valid_cells <- complete.cases(pseudotime)
#pseudotime <- pseudotime[valid_cells, ]
counts <- assays(sce_subsampled)$counts
cellWeights <- slingCurveWeights(sds)

# Ensure that age and sex are in the colData
if (!all(c("age", "sex") %in% colnames(colData(sce_subsampled)))) {
  stop("The colData must contain 'age' and 'sex' columns.")
}


# Create the design matrix
col_data <- colData(sce_subsampled)[,c("age", "sex")] 

design <- model.matrix(~ age + sex, data = col_data)
design <- design[,-1]

# Fit GAMs for each gene along the pseudotime
print("Running fitGAM...")
  sce_subsampled <- fitGAM(
    counts = counts,
    # genes = final_gene_subset,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    conditions = as.factor(sce_subsampled$diagnosis),
    U = design,
    nknots = 6,
    parallel = TRUE,
    BPPARAM = BPPARAM
  )
print("Finished fitGAM!")

# Save
file <- here::here("03_data/990_processed_data/008_pseudotime",
                   "slingshot_tradeseq_3k_pericyte.qs")
if(subset == "pericyte") {
  file <- str_replace(file, "_pericyte", "_smc")
} 

qs::qsave(sce_subsampled, file)

# Extract the model results
gam_results <- rowData(sce_subsampled)$tradeSeq

print("gam results head:")
print(head(gam_results))

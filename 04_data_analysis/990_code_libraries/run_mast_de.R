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
library(tidyverse)

# read data
sce <- qs::qread(here::here("03_data/990_processed_data/001_snrnaseq",
    "11_cell_network_interactions",
    "scflow-sce-annotated.qs"))

gene_symbols <- rownames(sce@assays$RNA@data)
rownames(sce@assays$RNA@data) <- rownames(sce@assays$RNA@counts)
sce <- subset(sce,
              idents = c("ambiguous", "low-feature-cells"),
              invert = TRUE
)

# Define a function to perform differential expression analysis for a given cell type
perform_DE <- function(seruat_obj, cell_type, level = "level2") {
  file <- here::here("03_data/990_processed_data/001_snrnaseq/13_mast_de",
               level, paste0(cell_type, "mast_de.qs"))
  if (!exists(file)) {
    # Perform differential expression analysis using MAST
    df <- FindMarkers(
      object = seruat_obj,
      ident.1 = "AD",
      ident.2 = "Control",
      group.by = "diagnosis",
      test.use = "MAST",
      latent.vars = c("donor", "age", "sex"),
      subset.ident = cell_type
    )
    df$celltype <- cell_type
    qs::qsave(df, file)
  } else {
    df <- qs::qread(file)
  }
  return(df)
  print(paste0("Finished: ", cell_type))
}

# setup columns to aggregate by
sce$celltype <- Idents(sce)
# remove the "_" in the celltype names and replace with "-"
sce$celltype <- gsub("_", "-", sce$celltype)
sce$diagnosis <- forcats::fct_recode(sce$diagnosis, AD = "Case")
# Use purrr::map to apply the perform_DE function to each cell type
# Set up parallel processing
plan(strategy = "multicore", workers = parallel::detectCores())
de_results <- purrr::map(unique(sce$celltype), perform_DE, seruat_obj = sce)

# Name the list elements by cell type
names(de_results) <- unique(sce$celltype)

# Save
qs::qsave(de_results, here::here("03_data/990_processed_data/001_snrnaseq/13_mast_de/level2/mast_de_results_list.qs"))

# level 1
Idents(sce) <- sce$highlevel_manual_annotations
de_results <- purrr::map(unique(sce$highlevel_manual_annotations), perform_DE, 
                         seruat_obj = sce, level = "leve1")

# Name the list elements by cell type
names(de_results) <- unique(sce$highlevel_manual_annotations)

# Save
qs::qsave(de_results, here::here("03_data/990_processed_data/001_snrnaseq/13_mast_de/level1/mast_de_results_list.qs"))
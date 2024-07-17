library(targets)
library(tarchetypes)
# library(crew)

# Set target options
tar_option_set(
  packages = c("tidyverse", "Seurat", "CellChat", "future", "biomaRt", "DoubletFinder"),
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE # ,
  # controller = crew_controller_local(workers = parallel::detectCores() - 2,
  #                                    seconds_idle = 10)
)

# This is an example _targets.R file. Every
# {targets} pipeline needs one.
# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(data_summary) to view the results.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
read_sce <- function(file) {
  sce <- qs::qread(file)
  gene_symbols <- rownames(sce@assays$RNA@data)
  rownames(sce@assays$RNA@data) <- rownames(sce@assays$RNA@counts)
  sce <- subset(sce,
    idents = c("ambiguous", "low-feature-cells"),
    invert = TRUE
  )
  rownames(sce@assays$RNA@data) <- gene_symbols
  return(sce)
}
get_ids <- function(sce) {
  data.frame(
    gene = rownames(sce@assays$RNA@data),
    ensembl = rownames(sce@assays$RNA@counts)
  )
}
translate_ids <- function(gene_ids) {
  # Get entrez IDs
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  entrez_ids <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = gene_ids$ensembl,
    mart = ensembl
  )
  return(entrez_ids)
}
make_cellchat <- function(sce, cellchat_db, case_or_control) {
  # Add a "samples" column that cellChat expects when making the object
  sce$samples <- sce$orig.ident

  df <- data.frame(
    gene = rownames(sce@assays$RNA@data),
    ensembl = rownames(sce@assays$RNA@counts)
  )
  # Subset to cases or controls
  rownames(sce@assays$RNA@data) <- df$ensembl
  sce <- subset(sce, subset = diagnosis == case_or_control)
  rownames(sce@assays$RNA@data) <- df$gene

  cellchat <- createCellChat(
    object = sce,
    group.by = "ident",
    assay = "RNA"
  )
  rm(sce, df)

  # use a subset of CellChatDB for cell-cell communication analysis
  cellchat@DB <- subsetDB(cellchat_db,
    search = "Secreted Signaling",
    key = "annotation"
  ) # use Secreted Signaling

  options(future.globals.maxSize = 3145728000)

  future::plan("sequential", workers = parallel::detectCores() - 2)

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}
cellchat_merge <- function(cellchat_list) {
  mergeCellChat(cellchat_list, add.names = names(cellchat_list))
}
seurat_list_by_sample <- function(seurat_file) {
  seurat <- readr::read_rds(seurat_file)
  seurat_list <- SplitObject(seurat, split.by = "manifest")
  return(seurat_list)
}
# Estimate best pK value for the sample
estimate_best_pk <- function(seurat_list) {

  estimate_best_pk_per_sample <- function(seurat, pc_num = 1:30) {
    sweep_res_list <- paramSweep(seurat,
                                 PCs = pc_num,
                                 sct = FALSE,
                                 num.cores = detectCores() - 1)
    sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    
    # get best pK value
    pK <- bcmvn %>%
      dplyr::filter(BCmetric == max(BCmetric)) %>%
      dplyr::select(pK)
    pK <- as.numeric(as.character(pK[[1]]))
    
    return(pK)
  }
  pks <- map(seurat_list, estimate_best_pk_per_sample)
  return(pks)
}
estimated_doublet_rate <- function(seurat_file) {
  seurat <- readr::read_rds(seurat_file)
  num_cells <- table(seurat$manifest) |> as.data.frame()
  
  # Expected multiplet rates from 10x Genomics
  multiplet_rates <- data.frame(
    recovered_cells = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 
                        9000, 10000), # Example values
    rate = c(0.004, 0.008, 0.016, 0.024, 0.032, 0.04, 0.048, 0.056, 0.064, 
             0.072, 0.08) # Example rates
  )
  
  num_cells <- num_cells |>
    dplyr::mutate(recovered_cells = round(Freq, -3)) |>
    dplyr::mutate(recovered_cells = if_else(recovered_cells == 0, 500, 
                                            recovered_cells)) |>
    dplyr::mutate(recovered_cells = if_else(recovered_cells > 10000, 10000, 
                                            recovered_cells)) |>
    dplyr::left_join(multiplet_rates)  
  return(num_cells)
}
# Function to run DoubletFinder on a single Seurat object
run_doubletfinder_on_list <- function(seurat_obj, pks, estimated_doublets) {
  
  run_doubletfinder <- function(seurat_obj, pks, estimated_doublets, 
                              num_pcs = 30) {
  # Run PCA if not already done
  if (!("pca" %in% names(seurat_obj@reductions))) {
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  }
  
  # Filter to the samples expected doublet rate row
  doublet_rate <- estimated_doublets |>
    dplyr::filter(Var1 %in% seurat_obj$manifest)
  
  # There should only be one row after this filter
  assertthat::assert_that(nrow(doublet_rate) == 1)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(Idents(seurat_obj))  
  nExp_poi <- round(doublet_rate$rate*nrow(seurat_obj@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
  
  print(paste0(
    "Expected number of doublets for sample ",
    unique(seurat_obj$manifest),
    ": ",
    nExp_poi.adj
  ))

  # Run DoubletFinder
  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = 1:num_pcs,
    pN = 0.25,
    pK = pks,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  
  print("Finished!")
  
  colnames(seurat_obj@meta.data)[35] <- "doublet_estimation"
  return(seurat_obj)
  }
  
  # Set up a future plan
  seurat_list <- map2(seurat_obj, pks, run_doubletfinder, 
                             estimated_doublets = estimated_doublets)
  seurat <- merge(x = seurat_list[[1]], y = seurat_list[-1])
  return(seurat)
}
# Set target-specific options such as packages:
# tar_option_set(packages = "utils") # nolint

# End this file with a list of target objects.
list(
  tar_target(file, here::here(
    "03_data/990_processed_data/001_snrnaseq",
    "11_cell_network_interactions",
    "scflow-sce-annotated.qs"
  ), format = "file"),
  tar_target(sce, read_sce(file)),
  tar_target(gene_ids, get_ids(sce)),
  tar_target(translated_id, translate_ids(gene_ids)),
  tar_target(cellchat_case, make_cellchat(sce, CellChatDB.human, "Case")),
  tar_target(cellchat_control, make_cellchat(sce, CellChatDB.human, "Control")),
  tar_target(seruat_qc_file, here::here(
    "03_data/990_processed_data/001_snrnaseq/07_scflow_analysis",
    "scflow-seurat-postprocessing_donors-excluded-updated.rds"
  ), format = "file"), 
  tar_target(estimated_doublets, estimated_doublet_rate(seruat_qc_file)),
  tar_target(seurat_list, seurat_list_by_sample(seruat_qc_file)),
  tar_target(pks, estimate_best_pk(seurat_list)),
  tar_target(seurat_doublets, run_doubletfinder_on_list(seurat_list, pks, 
                                                estimated_doublets)),
  tar_target(cellchat, cellchat_merge(list(
    case = cellchat_case,
    control = cellchat_control
  ))) # ,
  # tar_quarto(paper, "index.qmd")
)

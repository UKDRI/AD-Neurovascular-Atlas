library("EWCE", lib.loc = "/scratch/c.mpmfw/Tools/R/Older_R_packages")
library(scFlow)
library(monocle3) ## needed during annotate_merged_sce()
library(Seurat)
library(future)
library(pryr) ## memory check
library(data.table)
library(magrittr)

#############################################
## use Cell Ranger (CR) count calls to be filtered here
## i.e. do NOT use find_cells() from scFlow
## remove multiplets with DoubletFinder
## note, some scFlow code changes: keep "^MT" genes, if not mito genes, wrongly excluded by current code
## filtering settings: 
##  max_mito = 0.05 (default = adaptive)
##  min_features = 200 (default = 100)
##  max_features = 10000 (default = adaptive)
##  max_library_size = 80000 (default = adaptive)
## use default scFlow filters of:
##  min_counts = 2
##  min_cells = 2

#############################################
## input/output files

## metadata
samples_meta_file <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Samples_info/samples_meta_merged_final.txt"

## CR folders depending on set/batch
set_1_path <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/CellRanger_count/CellRanger_count_v611_ensembl_106_set_1_R2_full/"
set_2_path <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/CellRanger_count/CellRanger_count_v611_ensembl_106_set_2_R2_full/"

## generated via scFlow::map_ensembl_gene_id() and slightly modified 
## see get_gene_mapping_scflow.Rmd
## needed for annotate_sce()
map_ensembl_scflow_file <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Gene_info/map_ensembl_scflow_v106.txt"

## folder for sce objects after scflow QC
out_folder_sce <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Results_scFlow/SCE_objects_after_QC_CRv611_ensembl_106/"

## number of CPU cores to use
n_cores <- 2
options(future.globals.maxSize = 16000 * 1024^2)
plan(strategy = "multicore", workers = n_cores)

#############################################
## meta data of all samples
samples_meta <- fread(samples_meta_file, colClasses = c(donor = "character"))

## process input arguments: row index in samples meta data
args <- commandArgs(TRUE)
sample_row_idx <- as.numeric(args[1])
sample_meta_now <- samples_meta[sample_row_idx, ] ## one row


#############################################
## change scFlow mito genes handling

my_annotate_sce_genes <- function (sce, drop_unmapped = TRUE, drop_mito = TRUE, drop_ribo = FALSE, 
                                   ensembl_mapping_file = NULL, species = getOption("scflow_species", 
                                                                                    default = "human")) 
{
  cat(cli::rule("Annotating SingleCellExperiment genes", line = 2), 
      "\r\n")
  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }
  if (!("ensembl_gene_id" %in% colnames(SummarizedExperiment::rowData(sce)))) {
    stop(cli::cli_alert_danger("The rowData is missing ensembl_gene_id."))
  }
  mapping_results <- map_ensembl_gene_id(SummarizedExperiment::rowData(sce)$ensembl_gene_id, 
                                         ensembl_mapping_file = ensembl_mapping_file, species = species)
  mapping_results$ensembl_gene_id <- as.character(mapping_results$ensembl_gene_id)
  SummarizedExperiment::rowData(sce) <- dplyr::full_join(data.frame(SummarizedExperiment::rowData(sce)), 
                                                         mapping_results, by = "ensembl_gene_id", all = TRUE) %>% 
    dplyr::rename(gene = external_gene_name)
  if (!all(SummarizedExperiment::rowData(sce)$ensembl_gene_id == 
           rownames(sce))) {
    stop(cli::cli_alert_danger("Fatal error: Misaligned new_rowdata."))
  }
  SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped <- !(is.na(SummarizedExperiment::rowData(sce)$gene)) + 
    0
  ####################
  ## CHANGE
  qc_metric_is_mito <- grepl("^mt-|^MT-", as.character(SummarizedExperiment::rowData(sce)$gene))
  ####################
  qc_metric_is_mito[is.na(qc_metric_is_mito)] <- 0
  SummarizedExperiment::rowData(sce)$qc_metric_is_mito <- qc_metric_is_mito
  qc_metric_is_ribo <- grepl("^RPS|^Rps|^RPL|^Rpl", as.character(SummarizedExperiment::rowData(sce)$gene))
  qc_metric_is_ribo[is.na(qc_metric_is_ribo)] <- 0
  SummarizedExperiment::rowData(sce)$qc_metric_is_ribo <- qc_metric_is_ribo
  SummarizedExperiment::rowData(sce)$qc_metric_mapped_keep <- (SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped | 
                                                                 !drop_unmapped)
  SummarizedExperiment::rowData(sce)$qc_metric_mito_keep <- !(qc_metric_is_mito & 
                                                                drop_mito)
  SummarizedExperiment::rowData(sce)$qc_metric_ribo_keep <- !(qc_metric_is_ribo & 
                                                                drop_ribo)
  sce@metadata$scflow_steps$genes_annotated <- 1
  return(sce)
}


my_annotate_sce <- function (sce, min_library_size = 300, max_library_size = "adaptive", 
                             min_features = 100, max_features = "adaptive", max_mito = "adaptive", 
                             min_ribo = 0, max_ribo = 1, min_counts = 2, min_cells = 2, 
                             drop_unmapped = TRUE, drop_mito = TRUE, drop_ribo = FALSE, 
                             annotate_genes = TRUE, annotate_cells = TRUE, nmads = 4, 
                             ensembl_mapping_file = NULL, species = getOption("scflow_species", 
                                                                              default = "human")) 
{
  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }
  cat(cli::rule("Annotating SingleCellExperiment", line = 2), 
      "\r\n")
  before_coldata_colnames <- colnames(sce@colData)
  before_rowdata_colnames <- colnames(SummarizedExperiment::rowData(sce))
  qc_params <- setdiff(names(formals(annotate_sce)), c("sce", 
                                                       "ensembl_mapping_file", "annotate_genes", "annotate_cells"))
  qc_params_l <- purrr::map(qc_params, ~get(.))
  qc_params_l <- purrr::set_names(qc_params_l, qc_params)
  sce@metadata[["qc_params"]] <- qc_params_l
  if (("adaptive" %in% sce@metadata$qc_params) & !sce@metadata$scflow_steps$emptydrops_annotated) {
    cli::cli_alert_warning("To improve adaptive thresholding, first run emptyDrops!")
  }
  if (annotate_genes) {
    sce <- my_annotate_sce_genes(sce, drop_unmapped = drop_unmapped, 
                                 drop_mito = drop_mito, drop_ribo = drop_ribo, ensembl_mapping_file = ensembl_mapping_file)
  }
  if (annotate_cells) {
    sce <- annotate_sce_cells(sce, min_library_size = min_library_size, 
                              max_library_size = max_library_size, min_features = min_features, 
                              max_features = max_features, max_mito = max_mito, 
                              min_ribo = min_ribo, max_ribo = max_ribo, min_counts = min_counts, 
                              min_cells = min_cells, nmads = nmads)
  }
  else {
    if (!annotate_genes) {
      stop(cli::cli_alert_danger("Nothing to do. Specify gene/cell/both."))
    }
  }
  if (annotate_cells) {
    cli::cli_alert_success("SingleCellExperiment cells were successfully annotated with: \r\n")
    cli::cli_ul(setdiff(colnames(sce@colData), before_coldata_colnames))
  }
  if (annotate_genes) {
    cli::cli_alert_success("SingleCellExperiment genes were successfully annotated with: \r\n")
    cli::cli_ul(setdiff(colnames(SummarizedExperiment::rowData(sce)), 
                        before_rowdata_colnames))
  }
  cat(cli::rule("Generating QC plots for SingleCellExperiment", 
                line = 2), "\r\n")
  cli::cli_text("Generating QC plots and appending to metadata.")
  all_scflow_fns <- ls(getNamespace("scFlow"), all.names = TRUE)
  qc_plot_fns <- all_scflow_fns[startsWith(all_scflow_fns, 
                                           "scFlow:::.qc_plot_")]
  for (fn in qc_plot_fns) {
    sce <- get(fn)(sce)
  }
  cli::cli_text("Generating QC summary table and appending to metadata.")
  sce <- scFlow:::.qc_append_summary_table(sce)
  return(sce)
}


#############################################
## run scFlow

## note,  barcodes will be re-named in generate_sce() with first 2 columns of metadata as prefix:
## barcode <- paste(metadata[[1]], metadata[[2]], colnames(mat), sep = "_")
## 
## note, ensembl gene IDs are used for feature IDs by scFlow
## whereas Seurat with Read10X() extracts gene symbols from features.tsv.gz in 2nd column as default (gene.column = 2)
## (and then makes them unique with make.unique())

if (sample_meta_now[, batch] == "set_1"){
  input_folder <- paste0(set_1_path, sample_meta_now[, sample], "/outs/filtered_feature_bc_matrix/")
} else {
  input_folder <- paste0(set_2_path, sample_meta_now[, sample], "/outs/filtered_feature_bc_matrix/")
}


## keep default min_counts = 2 and min_cells = 2, include MT genes that are not mitochondrial genes ("MT-")
sce <- read_sparse_matrix(input_folder) %>%
  generate_sce(metadata = sample_meta_now) %>%
  my_annotate_sce(species = "human", 
                  ensembl_mapping_file = map_ensembl_scflow_file,
                  max_mito = 0.05,
                  min_features = 200,
                  max_features = 10000,
                  max_library_size = 80000,
                  drop_unmapped = FALSE,
                  drop_mito = TRUE) %>%
  filter_sce() %>%
  find_singlets(singlet_find_method = "doubletfinder") %>%
  filter_sce()


#############################################
## write SCE after QC
sce %>% write_sce(folder_path = paste0(out_folder_sce, sample_meta_now[, sample]), write_metadata = TRUE)


#############################################
cat("", fill = T)
pryr::mem_used()
cat("", fill = T)
devtools::session_info()


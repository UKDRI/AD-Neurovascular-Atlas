library("EWCE", lib.loc = "/scratch/c.mpmfw/Tools/R/Older_R_packages")
library(scFlow)
library(data.table)
library(scran)          
library(Seurat)         ## also loads Rtsne
library(stringr)        ## str_split
library(dplyr)
library(future)         ## use multiple cores
library(monocle3)       ## needed during annotate_merged_sce()

## Endo 10X single nuclei: set 1 + set 2
## 40 samples: 20 V + 20 P samples
## 20 CT + 20 AD samples
## input is from modified scFlow filtering

############################
## files

## folder with all SCE objects from scFLow, already de-swapped
sce_folder <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Results_scFlow/SCE_objects_after_QC_CRv611_ensembl_106/"

## meta data of all samples
samples_meta_file <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Samples_info/samples_meta_merged_final.txt"

## generated via scFlow::map_ensembl_gene_id() and slightly modified 
## see get_gene_mapping_scflow.Rmd
## needed for annotate_sce()
## corresponds to features.tsv.gz from CR counts output
map_ensembl_scflow_file <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Gene_info/map_ensembl_scflow_v106.txt"

## number of cores
n_workers <- 40

## save seurat object as rds
out_seurat <- "/scratch/c.mpmfw/Endo_10X_main/Analysis/Processed_data_final/sdata_endo_40_samples.rds"


#####################################
## Import SCE objects from scFLow output
## NOTE, cannot convert merged SCE object from scFLow to Seurat object:
## Error in checkSlotAssignment(object, name, value) : assignment of an object of class “dgTMatrix” is not valid for slot ‘data’ in an object of class “Assay”; is(value, "AnyMatrix") is not TRUE
## so, instead load and merge individual objects
samples_meta <- fread(samples_meta_file, colClasses = c(donor = "character"))
sce_paths <- as.list(paste0(sce_folder, samples_meta$sample))
sce_merged <- merge_sce(sce_paths, ensembl_mapping_file = map_ensembl_scflow_file)


##############
## change feature names from scFLow ensembl ID to gene symbol
## some gene symbols (<20) are not unique
## use make.unique as in Seurat::Read10X()
map_ensembl_scflow <- fread(map_ensembl_scflow_file)
r_dt <- data.table(data.frame(rowData(sce_merged)))
stopifnot(length(setdiff(r_dt$ensembl_gene_id, map_ensembl_scflow$ensembl_gene_id)) == 0)
r_dt <- plyr::join(r_dt, map_ensembl_scflow[, .(ensembl_gene_id, external_gene_name)], by = "ensembl_gene_id") ## keep order of genes
stopifnot(all.equal(rownames(sce_merged), r_dt$ensembl_gene_id))
rownames(sce_merged) <- make.unique(names = r_dt$external_gene_name)


###########################
## create Seurat object
sdata_endo <- as.Seurat(sce_merged, data = "counts")

## note, active assay is called "originalexp" instead of "RNA"
sdata_endo <- RenameAssays(object = sdata_endo, originalexp = 'RNA')

## change orig.ident to sample ID
sdata_endo$orig.ident <- sdata_endo$sample

## nuclei per sample
## sdata_endo_meta <- sdata_endo@meta.data
## table(sdata_endo_meta$orig.ident)

#########################
## standard seurat processing

## 16 GB
options(future.globals.maxSize = 16000 * 1024^2)
plan(strategy = "multicore", workers = n_workers)

## normalisation
sdata_endo <- NormalizeData(object = sdata_endo, normalization.method = "LogNormalize", scale.factor = 10000)

## variable genes
sdata_endo <- FindVariableFeatures(sdata_endo, selection.method = "vst", nfeatures = 2000)

## scaling
all.genes <- rownames(sdata_endo)
sdata_endo <- ScaleData(sdata_endo, features = all.genes)

## PCA
sdata_endo <- RunPCA(sdata_endo, features = VariableFeatures(object = sdata_endo))


####################
## save object
saveRDS(sdata_endo, out_seurat, compress = T)


#############################################
cat("", fill = T)
pryr::mem_used()
cat("", fill = T)
devtools::session_info()
#| label: load libraries
#| include: false

## create vector of packages used
pkg <- c("helpers" = c("stringr", "DT", "reshape", "dplyr", "tidyr",
                       "tibble", "purrr", "here", "readr", "kableExtra", 
                       "janitor", "labelled", "BiocManager", "data.table",
                       "gt", "europepmc", "readxl", "NMF"),
         "website" = c("downlit", "xml2"),
         "plots" = c("corrgram", "RColorBrewer", "pheatmap", "ggplot2",
                     "ggsci", "cowplot", "scales", "ggvenn", "ggpubr", 
                     "gridExtra", "tidytext", "viridis", "wesanderson",
                     "pals", "ggalluvial"))

pkg_bioconductor <- c("cellity", "scater", "Seurat", "SC3", "DESeq2", "edgeR",
                      "BiocParallel", "org.Hs.eg.db", "pathview", "gage",
                      "gageData", "DropletUtils", "ensembldb", 
                      "AnnotationHub", "patchwork", "scran", "scuttle",  
                      "PCAtools", "batchelor", "bluster", "igraph", 
                      "SingleCellExperiment", "scRNAseq", "biomaRt",
                      "GEOquery", "AnnotationDbi", "clusterProfiler",
                      "DOSE", "enrichplot", "EnsDb.Hsapiens.v79", "fgsea")

pkg_devtools <- c("sqjin/CellChat")

## create function to installing any missing packages and then load them all
sort_packages <- function(packages, bioconductor = FALSE, devtools = FALSE) {
    ## Check if packages are not installed and assign the
    ## names of the packages not installed to the variable new.pkg
    new_pkg <- packages[!(packages %in% installed.packages())]
    ## If the packages are bioconductor, install via biocmanager
    ## If there are any packages in the list that aren't installed,
    ## install them
    if (length(new_pkg) & !bioconductor) {
        install.packages(new_pkg, repos = "https://cran.rstudio.com")
    } else if (length(new_pkg) & bioconductor) {
        BiocManager::install(new_pkg)
    } else if (length(new_pkg) & devtools) {
      devtools::install_github(new_pkg)
    }
    
    if (devtools) {
      packages <- sapply(strsplit(packages, "/"), tail, n = 1)
    }
    ## load packages
    lapply(packages, require, character.only = TRUE)
}
## use function
sort_packages(pkg)
sort_packages(pkg_bioconductor, bioconductor = TRUE)
sort_packages(pkg_devtools, devtools = TRUE)





rm(pkg, pkg_bioconductor, pkg_devtools, sort_packages)


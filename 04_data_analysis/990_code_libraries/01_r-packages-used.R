#| label: load libraries
#| include: false

## create vector of packages used
pkg <- c("helpers" = c("stringr", "DT", "reshape", "dplyr", "tidyr",
                       "tibble", "purrr", "here", "readr", "kableExtra", 
                       "janitor", "labelled", "BiocManager"),
         "website" = c("downlit", "xml2"),
         "plots" = c("corrgram", "RColorBrewer", "pheatmap", "ggplot2",
                     "ggsci", "cowplot", "scales", "ggvenn", "ggpubr", 
                     "gridExtra", "tidytext"))

pkg_bioconductor <- c("cellity", "scater", "Seurat", "SC3", "DESeq2", "edgeR",
                      "BiocParallel", "org.Mm.eg.db", "pathview", "gage",
                      "gageData", "glmGamPoi", "DropletUtils", "ensembldb", 
                      "AnnotationHub", "patchwork", "scran", "scuttle",  
                      "PCAtools", "batchelor", "bluster", "igraph")

## create function to installing any missing packages and then load them all
sort_packages <- function(packages, bioconductor = FALSE) {
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
    }
    ## load packages
    lapply(packages, require, character.only = TRUE)
}
## use function
sort_packages(pkg)
sort_packages(pkg_bioconductor, bioconductor = TRUE)





rm(pkg, pkg_bioconductor, sort_packages)


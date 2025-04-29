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

# Set to run on cases, controls, or both
cases_or_controls <- "both"

file <- here::here(
  "03_data/990_processed_data/008_pseudotime",
  "slingshot_sds_subset_3k_controls.qs"
)

if (cases_or_controls == "cases") {
  file <- str_replace(file, "_controls", "_cases")
} else if (cases_or_controls == "both") {
  file <- str_replace(file, "_controls", "_case_and_control")
}


# Set up parallel processing
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 12
BPPARAM

# Evaluate the optimal number of knots
# file <- here::here("03_data/990_processed_data/008_pseudotime/slingshot_kvalues_1k.qs")
#
# if(!file.exists(file)) {
#
# set.seed(123)  # For reproducibility
# print("Running evaluateK...")
# k_values <- evaluateK(counts = assays(sce_subsampled)$counts,
#                       sds = sds, nGenes = 500, k = 3:10,
#                       parallel = TRUE, BPPARAM = BPPARAM)
#
# qs::qsave(k_values, file)
# print("Finished evaluateK!")
# } else {
#   k_values <- qs::qread(file)
# }

# Assuming kResults is your matrix from evaluateK
# Rows are genes, columns are different values of K
# summarized_metrics <- apply(kResults, 2, mean)  # You can use mean, median, etc.

# If you want to use a different summary statistic, e.g., median:
# summarized_metrics <- apply(k_values, 2, median)
# summarized_metrics[summarized_metrics == max(summarized_metrics)] |> names()
#
# # Create a data frame for plotting
# ks <- as.numeric(colnames(k_values))
# metrics_df <- data.frame(
#   K = seq_along((colnames(k_values))) + 2,
#   Metric = summarized_metrics
# )
#
# # Plot the summarized metric
# ggplot(metrics_df, aes(x = K, y = Metric)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Summarized Metric for Different K", x = "Number of Clusters (K)", y = "Summarized Metric")

sce_slingshot <- qs::qread(here::here(
  "03_data/990_processed_data/008_pseudotime",
  "slingshot_obj.qs"
))

pseudotime <- slingPseudotime(sce_slingshot, na = FALSE)
na_counts <- sum(is.na(pseudotime))
if (na_counts > 0) {
  stop("There are still NA values in the pseudotime matrix.")
}

counts <- assays(sce_slingshot)$counts
cellWeights <- slingCurveWeights(sds)


# Ensure that age and sex are in the colData
if (!all(c("age", "sex") %in% colnames(colData(sce_slingshot)))) {
  stop("The colData must contain 'age' and 'sex' columns.")
}


# Create the design matrix
col_data <- colData(sce_slingshot)[, c("age", "sex")]

design <- model.matrix(~ age + sex, data = col_data)
design <- design[, -1]

# Fit GAMs for each gene along the pseudotime
print("Running fitGAM...")
if (cases_or_controls == "both") {
  sce_subsampled <- fitGAM(
    counts = counts,
    # genes = final_gene_subset,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    conditions = as.factor(sce_slingshot$diagnosis),
    U = design,
    nknots = 6,
    parallel = TRUE,
    BPPARAM = BPPARAM
  )
} else {
  sce_subsampled <- fitGAM(
    counts = counts,
    # genes = final_gene_subset,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    U = design,
    nknots = 6,
    parallel = TRUE,
    BPPARAM = BPPARAM
  )
}
print("Finished fitGAM!")

# Save
file <- here::here(
  "03_data/990_processed_data/008_pseudotime",
  "slingshot_tradeseq_3k_controls.qs"
)
if (cases_or_controls == "cases") {
  file <- str_replace(file, "_controls", "_cases")
} else if (cases_or_controls == "both") {
  file <- str_replace(file, "_controls", "_case_and_control")
}
qs::qsave(sce_subsampled, file)

# Extract the model results
gam_results <- rowData(sce_subsampled)$tradeSeq

print("gam results head:")
print(head(gam_results))

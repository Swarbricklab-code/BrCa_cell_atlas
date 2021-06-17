# ===================
# title: "dwls_deconvolute_bulk"
# description: runs DWLS cell-type deconvolution of bulk samples using DWLS cell-type signature matrix 
# project: "BrCa-atlas"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# created: 01/12/2020
# version: 0.1
# ===================
debug <- FALSE

file_log <- file(snakemake@log[[1]], open="wt")
sink(file_log)
sink(file_log, type="message")

message(paste0("LOG FILE: ", snakemake@log[[1]]))

start_time <- Sys.time()
message("START TIME: ", start_time)
message("\n----------------------------------") 
message("CONFIG")

# libraries

# ================
# DWLS
# ================
# Tsoucas, D. et al. Accurate estimation of cell-type composition from gene expression data. Nat Commun 10, 2975 (2019). 
# https://www.nature.com/articles/s41467-019-10802-z (DOI: 10.1038/s41467-019-10802-z)
# code obtained from here (30/11/2020): https://github.com/dtsoucas/DWLS 
# https://github.com/dtsoucas/DWLS/blob/master/Deconvolution_functions.R
source("code/external/DWLS/Deconvolution_functions.R")

message("INPUT")
#file_signature_matrix <- "testing/DWLS/results/signature_matrix.txt"
file_signature_matrix <- snakemake@input[["signature_matrix"]]
message("file_signature_matrix: ", file_signature_matrix)

# bulk dataset
file_bulk_data <- snakemake@input[["bulk_data"]]
message("file_bulk_data: ", file_bulk_data)
message("")

message("OUTPUT")
file_celltype_fractions <- snakemake@output[["celltype_fractions"]]
message("file_celltype_fractions: ", file_celltype_fractions)
message("")

message("WILDCARDS")
message("")

message("PARAMS")
message("")

message("--------------------------------")
message("ANALYSIS ......") # ----

# Signature Matrix
# NOTE: needs to be a matrix for solveDampenedWLS() function
df_signature_matrix <- as.matrix(read.delim(file_signature_matrix, row.names = 1))
celltypes <- colnames(df_signature_matrix)
message("loaded signature matrix with: ", ncol(df_signature_matrix), " cell-types and ", nrow(df_signature_matrix), " genes")

# Load Bulk data
df_bulk_data <- read.delim(file_bulk_data, row.names = 1)
message("loaded bulk dataset with: ", ncol(df_bulk_data), " samples and ", nrow(df_bulk_data), " genes")

message("RUNNING DWLS Deconvolution")
sample_ids <- colnames(df_bulk_data)

df_cell_fractions <- NULL
sample_count <- 1
for(sample_id in sample_ids)
{
  message("SAMPLE_ID = ", sample_id, " (", sample_count, "/", length(sample_ids), ")")
  bulk_sample <- df_bulk_data[, sample_id]
  names(bulk_sample) <- rownames(df_bulk_data)
  
  message("trim signature and bulk data to contain the same differentially expressed genes")
  signature_matrix_trimmed <- trimData(df_signature_matrix, bulk_sample)
  message("Trimmed data contains: ", length(signature_matrix_trimmed$bulk), " genes")
  
  # estimate using dampened weighted least squares
  solDWLS <- solveDampenedWLS(signature_matrix_trimmed$sig, signature_matrix_trimmed$bulk)
  # round results
  solDWLS <- round(solDWLS, 5)
  
  solDWLS <- data.frame(t(solDWLS))
  solDWLS[, "SAMPLE_ID"] <- sample_id
  df_cell_fractions <- rbind(df_cell_fractions, solDWLS)  
  
  sample_count <- sample_count + 1
}

message("re-arrange columns and save results")
df_cell_fractions <- df_cell_fractions[, c("SAMPLE_ID", celltypes)]
write.table(x = df_cell_fractions, file = file_celltype_fractions, sep = "\t", row.names = F)

message("--------------------------------")
message("FINALISE") # ----
finish_time <- Sys.time()
print(paste0("start time = ", start_time))
print(paste0("finish time = ", finish_time))

message("--------------------------------")
message("WARNINGS") # ----
warnings()

message("--------------------------------")
message("SESSION INFO") # ----
sessionInfo()

# ===================
# title: "dwls_generate_signature_matrix"
# description: runs generation of cell-type signature matrix for DWLS cell-type deconvolution 
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
library(Matrix)

# ================
# DWLS
# ================
# Tsoucas, D. et al. Accurate estimation of cell-type composition from gene expression data. Nat Commun 10, 2975 (2019). 
# https://www.nature.com/articles/s41467-019-10802-z (DOI: 10.1038/s41467-019-10802-z)
# code obtained from here (30/11/2020): https://github.com/dtsoucas/DWLS 
# https://github.com/dtsoucas/DWLS/blob/master/Deconvolution_functions.R
source("code/external/DWLS/Deconvolution_functions.R")

message("INPUT")
# single-cell data ("dense" counts matrix)
file_count_matrix <- snakemake@input[["count_matrix"]]
message("file_count_matrix: ", file_count_matrix)

# vector of cell-types
file_cell_ids <- snakemake@input[["cell_ids"]]
message("file_cell_ids: ", file_cell_ids)
message("")

message("OUTPUT")
file_signature_matrix <- snakemake@output[["signature_matrix"]]
message("file_signature_matrix: ", file_signature_matrix)

dir_out <- dirname(file_signature_matrix)
message("DWLS output directory: ", dir_out)

message("")

message("WILDCARDS")
message("")

message("PARAMS")
message("")

message("--------------------------------")
message("ANALYSIS ......") # ----

# ANALYSIS STEPS
message("LOAD COUNTS")
mat_counts <- read.delim(file_count_matrix, row.names = 1)
cell_ids <- scan(file_cell_ids, what = "character")

message("cell-types identified and frequency of coccurence:")
print(table(cell_ids))

message("CREATE SIGNATURE MATRIX")
signature_matrix <- buildSignatureMatrixMAST(mat_counts, cell_ids, dir_out)
celltypes <- colnames(signature_matrix)

message("signature matrix contains: ", ncol(signature_matrix), " cell-types and ", nrow(signature_matrix), " genes")

message("SAVING SIGNATURE MATRIX")
df_signature_matrix <- data.frame(signature_matrix)
df_signature_matrix[, "GENE_SYMBOL"] <- rownames(df_signature_matrix)
df_signature_matrix <- df_signature_matrix[, c("GENE_SYMBOL", celltypes)]

write.table(file = file_signature_matrix, x = df_signature_matrix, sep = "\t", row.names = F)

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

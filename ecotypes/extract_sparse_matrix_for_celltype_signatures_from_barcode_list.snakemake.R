# ==========================================
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# date: "`r Sys.Date()`"
# created: 04/08/2020
# version: 1.0
#
# Overview
# Extract sparse matrix from seurat (using cell barcode list)
# ==========================================

# Config
library(Matrix)
library(Seurat)
library(dplyr)

# PARAMS

# INPUT
file_seurat <- snakemake@input[["seurat"]]
file_barcode_list <- snakemake@input[["barcode_list"]]

# OUTPUT
file_cells <- snakemake@output[["cells"]]
file_genes <- snakemake@output[["genes"]]
file_sparse_matrix <- snakemake@output[["sparse_matrix"]]

# WILDCARDS
subsample_fraction <- snakemake@wildcards[["subsample_fraction"]]

# ================================
# read the integrated dataset
sData <- readRDS(file_seurat)

# and barcode list
df_barcode_list <- read.delim(file = file_barcode_list, sep = "\t", stringsAsFactors = F)

# -----------------------
# Randomly sample from each cell-type
# -----------------------
message("BEFORE SUB SAMPLE")
table(df_barcode_list[, 2])

# Fix the cell-type annotation names for CIBERSORT input
# remove "+"
df_barcode_list[, 2] <- gsub(pattern = "\\+", replacement = "", x = df_barcode_list[, 2])
# remove " "
df_barcode_list[, 2] <- gsub(pattern = "\\ ", replacement = "_", x = df_barcode_list[, 2])   
# remove "-"
df_barcode_list[, 2] <- gsub(pattern = "\\-", replacement = "", x = df_barcode_list[, 2])  
# append "x" to cell-types that end with numbers (R annoyingly appends .n to each one when reading if not - see: https://stackoverflow.com/questions/13235163/how-do-i-stop-r-from-prepending-x-to-my-numeric-column-labels)

celltype_sym <- rlang::sym(colnames(df_barcode_list)[2])

if(subsample_fraction != "all")
{
  subsample_fraction <- as.numeric(subsample_fraction)
  set.seed(123)
  df_sampled <- df_barcode_list %>% 
    group_by(!! celltype_sym) %>%   
    sample_frac(subsample_fraction)
  
  message("AFTER SUB SAMPLE")
  print(table(df_sampled[, 2]))
  
  cell_ids <- df_sampled %>% dplyr::pull("barcode")
  sData <- SubsetData(sData, cells = cell_ids)
  
  # Check we have all expected barcodes
  if(length(cell_ids) == ncol(sData)) {
    message(paste0("Extracted all ", ncol(sData), " barcodes"))
  } else {
    stop("ERROR: not all barcodes in barcode list are found in seurat object. Stopping.")
  }
  
} else {
  message("NO SUB-SAMPLING DONE.")  
}

# -----------------------
# generate count matrix
# -----------------------
sparseMatrix <- sData@assays$RNA@counts

# write as sparse matrix first (R generates errors when using as.matrix())
rownames(df_barcode_list) <- df_barcode_list$barcode
cell_ids <- df_barcode_list[colnames(sparseMatrix), 2]  
write(cell_ids, file = file_cells)
write(rownames(sparseMatrix), file = file_genes)
Matrix::writeMM(sparseMatrix, file = file_sparse_matrix)

# NOW, use python "convert_sparse_matrix_to_dense.py" script to write to standard matrix

# R Session Info
sessionInfo()

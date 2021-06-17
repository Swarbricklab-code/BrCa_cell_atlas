# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# date: "`r Sys.Date()`"
# created: 15/05/2020
# version: 1.0
#
# Overview
# Generate mixture matrix from the scRNA-Seq (i.e., psuedo-bulk) for CIBERSORTx analysis of brca-atlas

# Config
library(Matrix)
library(Seurat)
library(dplyr)

# PARAMS

# INPUT
file_seurat <- snakemake@input[["seurat"]]

file_avg_matching_pseudobulk_mixture_matrix <- snakemake@output[["avg_matching_pseudobulk_mixture_matrix"]]
file_sum_matching_pseudobulk_mixture_matrix <- snakemake@output[["sum_matching_pseudobulk_mixture_matrix"]]

# read the integrated dataset
sData <- readRDS(file_seurat)

# =============================
# Generate matching pseudo-bulk for each tumour
# -----------------------
# using Average expression in non-log count space
# -----------------------
sData <- SetIdent(sData, value = "orig.ident") 
matAvgMatchingPseudoBulk <- Seurat::AverageExpression(sData, assays = "RNA", slot = "counts", return.seurat = F)

# write mixture matrix
write.table(matAvgMatchingPseudoBulk, file = file_avg_matching_pseudobulk_mixture_matrix, sep = "\t", col.names = NA, row.names = TRUE)

# -----------------------
# using summed expression in non-log count space
# -----------------------
dfMeta <- sData@meta.data
dfMeta[, "barcode"] <- rownames(dfMeta)
sampleIDs <- unique(dfMeta$orig.ident)

dfAllSum <- NULL
for(sampleID in sampleIDs)
{
  cell_barcodes <- dfMeta[dfMeta[, "orig.ident"] == sampleID, "barcode"]
  sSample <- Seurat::SubsetData(sData, cells = cell_barcodes)
  message(paste0("Extracted ", length(cell_barcodes), " cells and ", nrow(sSample), " genes from SAMPLE: ", sampleID))
  
  # get summed gene counts
  dfSum <- as.data.frame(rowSums(sSample@assays$RNA@counts))
  colnames(dfSum) <- c(sampleID)
  dfSum[, "gene"] <- rownames(dfSum)
  if(is.null(dfAllSum)) {
    dfAllSum <- dfSum
  } else {
    dfAllSum <- dplyr::left_join(dfAllSum, dfSum, by = "gene")
  }
}
# rearrange columns to fit CIBERSORTx format
dfAllSum <- dfAllSum[, c("gene", sampleIDs)]

# write mixture matrix
write.table(dfAllSum, file = file_sum_matching_pseudobulk_mixture_matrix, sep = "\t", row.names = F)

# -----------------------
message("R Session Info")
# -----------------------
sessionInfo()


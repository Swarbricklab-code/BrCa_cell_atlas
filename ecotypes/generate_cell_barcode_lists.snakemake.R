# ===================
# title: "generate_cell_barcode_lists"
# project: "brca_atlas"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# date: "`r Sys.Date()`"
# version: 1.0
# ----------------------
# DESCRIPTION: generates the following cell barcode lists for use in DECON
# 
# CONFIG  ------------------------------------------------------------
message("CONFIG")
library(Seurat)

temp_start_time <- Sys.time()
print(temp_start_time)

# INPUT
seurat_file <- snakemake@input[["seurat"]]

# OUTPUT
barcode_list_file <- snakemake@output[["barcode_list"]]

# WILDCARDS ----
normal_epithelial_cells <- snakemake@wildcards[["normal_epithelial_cells"]]
celltype <- snakemake@wildcards[["cell_type"]]
cancer_cell_annotation <- snakemake@wildcards[["cancer_cell_annotation"]]
cycling_cell_annotation <- snakemake@wildcards[["cycling_cell_annotation"]]

# -----
message("READ SEURAT")
seurat_10X <- readRDS(seurat_file)
dfMeta <- seurat_10X@meta.data
dfMeta$barcode <- rownames(dfMeta)
message(paste0("PARSED OBJECT of ", nrow(dfMeta), " cells"))  

# remove spaces from cell-type names
dfMeta[, celltype] <- gsub(pattern = " ", replacement = "_", x = dfMeta[, celltype])

# -----
message("EXTRACT CELL TYPES")

# include normal epithelial cells?
if(normal_epithelial_cells == "without_normal_epithelial_cells") 
{
  dfMeta <- dfMeta[dfMeta$normal_cell_call != "normal", ]
  message(paste0("removed normal epithelial cells"))  
  message(paste0("EXTRACTED ", nrow(dfMeta), " cells"))  
} else if(normal_epithelial_cells == "with_normal_epithelial_cells") {
  message("no normal cell filtering done")
} else {
  stop(paste0("unknown normal_epithelial_cells wildcard: ", normal_epithelial_cells))
}

if(celltype == "celltype_major")
{
  # no further processing needed
} else {
  # cancer cell annotations
  if(cancer_cell_annotation == "SCTyper_Dec") {
    # NO change needed for now as the cancer cells are already annotated with SCTyper calls
    message("Annotated cancer cells with current SCTyper call")    
  } else if(cancer_cell_annotation == "ith_gene_modules") {
    cancer_cells <- !(grepl(pattern = "no_gene_module", x = dfMeta[, "gene_module"]))
    dfMeta[cancer_cells, celltype] <- paste0("GM", dfMeta[cancer_cells, "gene_module"])    
    message("Annotated cancer cells with ITH gene module annotaions")
  }
  
  # cycling cells
  if(cycling_cell_annotation == "separate") {
    message("no changes made to cycling cell annotations")
  } else if(cycling_cell_annotation == "combined") {
    # combine all cycling annotations into one (this can help prevent confounding in deconvolution)
    cycling_cells <- grepl(pattern = "Cycling", x = dfMeta[, celltype]) | grepl(pattern = "MKI67", x = dfMeta[, celltype])
    dfMeta[cycling_cells, celltype] <- "Cycling"
    message("combined all cycling cells")
  } else if(cycling_cell_annotation == "cancer_cell_annotation") {
    if(cancer_cell_annotation == "SCTyper_Dec")
    {
      # replace cancer cell cycling annotation with the original SCTyper call
      cycling_cancer_cells <- grepl(pattern = "Cancer_Cycling", x = dfMeta[, celltype])
      dfMeta[cycling_cancer_cells, celltype] <- paste0("Cancer_", dfMeta[cycling_cancer_cells, "Calls"])
      message("Annotated cycling cancer cells with SCTyper call")
    } else if(cancer_cell_annotation == "ith_gene_modules") {
      # Nothing further to be done
      message("Annotated cycling cancer cells with ITH gene modules")      
    }
  } else {
    stop(paste0("unknown cycling_cell_annotation wildcard: ", cycling_cell_annotation))    
  }
}

# write the barcode list
dfMeta <- dfMeta[, c("barcode", celltype), drop = F] 
colnames(dfMeta)[2] <- "annotation"
message(paste0("EXTRACTED ", nrow(dfMeta), " cells"))  
table(dfMeta[, "annotation"])

write.table(dfMeta, file = barcode_list_file, row.names = F, sep = "\t")

# SESSIONINFO ----------------------------------------------------------------------
message("SESSIONINFO")

temp_finish_time <- Sys.time()
print(paste0("start time = ", temp_start_time))
print(paste0("finish time = ", temp_finish_time))

sessionInfo()

# ==========================================
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# date: "`r Sys.Date()`"
# created: 15/05/2020
# version: 1.0
#
# Overview
# CIBERSORTx analysis
# ==========================================

# Config
library(Matrix)
library(Seurat)
source("code/R/plot_utils.R")
library(dplyr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(cola)
library(ComplexHeatmap)
library(randomcoloR)

# PARAMS
analysis_id <- "Jul2020"

celltype_classification <- "celltype_subset"

normal_epithelial_cells <- "with_normal_epithelial_cells"
cancer_cell_annotation <- "SCTyper_Dec"
cycling_cell_annotation <- "combined"
subsample_fraction <- 0.15
run_consensus_clustering <- TRUE

# ------------
# FUNCTIONS
# ------------
calc_cellfractions_persample_from_barcode_list <- function(df_barcode_list)
{
  # split into samples
  df_barcode_list <- df_barcode_list
  df_barcode_list[, "SampleID"] <- stringr::str_split_fixed(df_barcode_list$barcode, pattern = "_", n = 2)[, 1]
  
  # calc fraction per sample
  df_cell_fractions <- df_barcode_list %>% 
    dplyr::group_by(SampleID, annotation) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(fraction = n / sum(n))
  
  # wide form
  df_cell_fractions <- spread(df_cell_fractions[, c("SampleID", "annotation", "fraction")], key = annotation, value = fraction)
  df_cell_fractions <- as.data.frame(df_cell_fractions)
  df_cell_fractions[is.na(df_cell_fractions)] <- 0
  
  return(df_cell_fractions)
}

# ------------------
# minExp_0.75.rep_5
# ------------------
# INPUT
file_barcode_list <- paste0("analysis/CIBERSORTx/", 
                            analysis_id, "/",
                            celltype_classification, "/",
                            normal_epithelial_cells, "/",
                            cancer_cell_annotation, "/",
                            "cycling_", cycling_cell_annotation, "/",
                            "barcode_list.txt")

# pseudo-bulk
file_pseudobulk_results_none <- paste0("analysis/CIBERSORTx/", 
                                       analysis_id, "/",
                                       celltype_classification, "/",
                                       normal_epithelial_cells, "/",
                                       cancer_cell_annotation, "/",
                                       "cycling_", cycling_cell_annotation, "/",                                       
                                       "sampled_", subsample_fraction, "/",
                                       "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/pseudo_bulk/sum_matching/rmbatch_none/perm_0/relative/CIBERSORTx_Results.txt")
file_pseudobulk_results_smode <- paste0("analysis/CIBERSORTx/", 
                                        analysis_id, "/",
                                        celltype_classification, "/",
                                        normal_epithelial_cells, "/",
                                        cancer_cell_annotation, "/",
                                        "cycling_", cycling_cell_annotation, "/",                                       
                                        "sampled_", subsample_fraction, "/",                                        
                                        "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/pseudo_bulk/sum_matching/rmbatchSmode/refsample_counts/perm_0/relative/CIBERSORTx_Adjusted.txt")

# METABRIC
file_metabric_PAM50_info <- paste0("data/METABRIC/heloisa/METABRIC_ClassLabels_PAM50.txt")

# discovery
file_metabric_discovery_results_none <- paste0("analysis/CIBERSORTx/", 
                                               analysis_id, "/",
                                               celltype_classification, "/",
                                               normal_epithelial_cells, "/",
                                               cancer_cell_annotation, "/",
                                               "cycling_", cycling_cell_annotation, "/",                                       
                                               "sampled_", subsample_fraction, "/",                                               
                                               "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/discovery/rmbatch_none/QN_FALSE/perm_0/relative/CIBERSORTx_Results.txt")
file_metabric_discovery_results_smode <- paste0("analysis/CIBERSORTx/", 
                                                analysis_id, "/",
                                                celltype_classification, "/",
                                                normal_epithelial_cells, "/",
                                                cancer_cell_annotation, "/",
                                                "cycling_", cycling_cell_annotation, "/",                                       
                                                "sampled_", subsample_fraction, "/",                                                
                                                "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/discovery/rmbatchSmode/refsample_counts/QN_FALSE/perm_0/relative/CIBERSORTx_Adjusted.txt")

# validation
file_metabric_validation_results_none <- paste0("analysis/CIBERSORTx/", 
                                                analysis_id, "/",
                                                celltype_classification, "/",
                                                normal_epithelial_cells, "/",
                                                cancer_cell_annotation, "/",
                                                "cycling_", cycling_cell_annotation, "/",                                       
                                                "sampled_", subsample_fraction, "/",
                                                "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/validation/rmbatch_none/QN_FALSE/perm_0/relative/CIBERSORTx_Results.txt")
file_metabric_validation_results_smode <- paste0("analysis/CIBERSORTx/", 
                                                 analysis_id, "/",
                                                 celltype_classification, "/",
                                                 normal_epithelial_cells, "/",
                                                 cancer_cell_annotation, "/",
                                                 "cycling_", cycling_cell_annotation, "/",                                       
                                                 "sampled_", subsample_fraction, "/",                                                 
                                                 "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/validation/rmbatchSmode/refsample_counts/QN_FALSE/perm_0/relative/CIBERSORTx_Adjusted.txt")

# meta-data
file_ctp_clinical <- "data/sample_clinical_molecular_subtype_info.txt"
file_pseudobulk_PAM50 <- "data/pseudo_bulk/Allcells_PAM50files/SwarbrickAll_pam50scores.txt"

# OUTPUT
dirOutRoot <- paste0("analysis/CIBERSORTx/analysis/", 
                     analysis_id, "/",
                     celltype_classification, "/",
                     normal_epithelial_cells, "/",
                     cancer_cell_annotation, "/",
                     "cycling_", cycling_cell_annotation, "/",                                       
                     "sampled_", subsample_fraction, "/")
if(!dir.exists(dirOutRoot)) {
  dir.create(dirOutRoot, recursive = T, showWarnings = F)
}

# ======================
# PARSE CTP clinical
# ======================
df_ctp_clinical <- read.delim(file_ctp_clinical, stringsAsFactors = F)

df_pseudobulk_PAM50 <- read.delim(file_pseudobulk_PAM50, stringsAsFactors = F)
colnames(df_pseudobulk_PAM50)[1] <- "SampleID"
colnames(df_pseudobulk_PAM50)[7] <- "Allcells_PAM50"

# ======================
# Actual cell-fractions
# ======================
df_barcode_list <- read.delim(file_barcode_list, stringsAsFactors = F)
dfActual <- calc_cellfractions_persample_from_barcode_list(df_barcode_list)
colnames(dfActual) <- gsub(pattern = "-", replacement = "", colnames(dfActual))
colnames(dfActual) <- gsub(pattern = "\\+", replacement = "", colnames(dfActual))

cell_type_order <- NULL
celltype_name_fix <- NULL
# set the cell-type ordering for plot formatting and ordering
if(celltype_classification == "celltype_major") 
{
  if(normal_epithelial_cells == "with_normal_epithelial_cells") {
    cell_type_order <- c("Endothelial", "CAFs", "SMCs", "Normal_Epithelial", "Cancer_Epithelial", "Bcells", "Plasmablasts", "Tcells", "Myeloid")      
  } else if(normal_epithelial_cells == "without_normal_epithelial_cells") {
    cell_type_order <- c("Endothelial", "CAFs", "SMCs", "Cancer_Epithelial", "Bcells", "Plasmablasts", "Tcells", "Myeloid")
  } else {
    stop(paste0("ERROR: unknown normal_epithelial_cells: ", normal_epithelial_cells))
  }
    
} else if(celltype_classification == "celltype_minor") {
  if(normal_epithelial_cells == "with_normal_epithelial_cells") {
    cell_type_order <- c("Cancer_LumA_SC", "Cancer_LumB_SC", "Cancer_Her2_SC", "Cancer_Basal_SC", "Cycling",
                         "Myoepithelial", "Luminal_Progenitors", "Mature_Luminal",
                         "CAFs_MSC_iCAFlike", "CAFs_myCAFlike", 
                         "PVL_Differentiated", "PVL_Immature",
                         "Endothelial_ACKR1", "Endothelial_RGS5", "Endothelial_CXCL12",  "Endothelial_Lymphatic_LYVE1", 
                         "T_cells_CD8", "T_cells_CD4", 
                         "NK_cells", "NKT_cells", 
                         "B_cells_Memory", "B_cells_Naive", "Plasmablasts", 
                         "Macrophage", "Monocyte", "DCs")    
  } else if(normal_epithelial_cells == "without_normal_epithelial_cells") {
    cell_type_order <- c("Cancer_LumA_SC", "Cancer_LumB_SC", "Cancer_Her2_SC", "Cancer_Basal_SC", "Cycling",
                         "CAFs_MSC_iCAFlike", "CAFs_myCAFlike", 
                         "PVL_Differentiated", "PVL_Immature",
                         "Endothelial_ACKR1", "Endothelial_RGS5", "Endothelial_CXCL12",  "Endothelial_Lymphatic_LYVE1", 
                         "T_cells_CD8", "T_cells_CD4", 
                         "NK_cells", "NKT_cells", 
                         "B_cells_Memory", "B_cells_Naive", "Plasmablasts", 
                         "Macrophage", "Monocyte", "DCs")
  } else {
    stop(paste0("ERROR: unknown normal_epithelial_cells: ", normal_epithelial_cells))
  }
  
} else if(celltype_classification == "celltype_subset") {
  
  if(normal_epithelial_cells == "with_normal_epithelial_cells") {
    cell_type_order <- c("Cancer_LumA_SC", "Cancer_LumB_SC", "Cancer_Her2_SC", "Cancer_Basal_SC", 
                         "Cycling",
                         "Myoepithelial", "Luminal_Progenitors", "Mature_Luminal",                           
                         "CAFs_MSC_iCAFlike_s1", "CAFs_MSC_iCAFlike_s2", "CAFs_myCAF_like_s4", "CAFs_myCAF_like_s5", "CAFs_Transitioning_s3",
                         "Endothelial_ACKR1", "Endothelial_CXCL12", "Endothelial_Lymphatic_LYVE1", "Endothelial_RGS5",
                         "PVL_Immature_s1", "PVL_Immature_s2", "PVL_Differentiated_s3",
                         "T_cells_c0_CD4_CCR7", "T_cells_c1_CD4_IL7R", "T_cells_c2_CD4_Tregs_FOXP3", "T_cells_c3_CD4_Tfh_CXCL13",     
                         "T_cells_c4_CD8_ZFP36", "T_cells_c5_CD8_GZMK",  "T_cells_c7_CD8_IFNG", "T_cells_c8_CD8_LAG3",
                         "T_cells_c6_IFIT1", "T_cells_c9_NK_cells_AREG", "T_cells_c10_NKT_cells_FCGR3A",
                         "B_cells_Naive", "B_cells_Memory", "Plasmablasts",
                         "Myeloid_c0_DC_LAMP3", "Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE", "Myeloid_c3_cDC1_CLEC9A", "Myeloid_c4_DCs_pDC_IRF7", "Myeloid_c5_Macrophage_3_SIGLEC1", "Myeloid_c7_Monocyte_3_FCGR3A", "Myeloid_c8_Monocyte_2_S100A9", "Myeloid_c9_Macrophage_2_CXCL10", "Myeloid_c10_Macrophage_1_EGR1", "Myeloid_c11_cDC2_CD1C", "Myeloid_c12_Monocyte_1_IL1B")
    
  } else if(normal_epithelial_cells == "without_normal_epithelial_cells") {
    cell_type_order <- c("Cancer_LumA_SC", "Cancer_LumB_SC", "Cancer_Her2_SC", "Cancer_Basal_SC", 
                         "Cycling",
                         "CAFs_MSC_iCAFlike_s1", "CAFs_MSC_iCAFlike_s2", "CAFs_myCAF_like_s4", "CAFs_myCAF_like_s5", "CAFs_Transitioning_s3",
                         "Endothelial_ACKR1", "Endothelial_CXCL12", "Endothelial_Lymphatic_LYVE1", "Endothelial_RGS5",
                         "PVL_Immature_s1", "PVL_Immature_s2", "PVL_Differentiated_s3",
                         "T_cells_c0_CD4_CCR7", "T_cells_c1_CD4_IL7R", "T_cells_c2_CD4_Tregs_FOXP3", "T_cells_c3_CD4_Tfh_CXCL13",     
                         "T_cells_c4_CD8_ZFP36", "T_cells_c5_CD8_GZMK",  "T_cells_c7_CD8_IFNG", "T_cells_c8_CD8_LAG3",
                         "T_cells_c6_IFIT1", "T_cells_c9_NK_cells_AREG", "T_cells_c10_NKT_cells_FCGR3A",
                         "B_cells_Naive", "B_cells_Memory", "Plasmablasts",
                         "Myeloid_c0_DC_LAMP3", "Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE", "Myeloid_c3_cDC1_CLEC9A", "Myeloid_c4_DCs_pDC_IRF7", "Myeloid_c5_Macrophage_3_SIGLEC1", "Myeloid_c7_Monocyte_3_FCGR3A", "Myeloid_c8_Monocyte_2_S100A9", "Myeloid_c9_Macrophage_2_CXCL10", "Myeloid_c10_Macrophage_1_EGR1", "Myeloid_c11_cDC2_CD1C", "Myeloid_c12_Monocyte_1_IL1B")
  } else {
    stop(paste0("ERROR: unknown normal_epithelial_cells: ", normal_epithelial_cells))
  }
}

# ========================
# compare to pseudo-bulk
# ========================

# -----------------------
# PARSE
# -----------------------
df_pseudobulk_results_none <- read.delim(file_pseudobulk_results_none, stringsAsFactors = F)
colnames(df_pseudobulk_results_none)[1] <- "SampleID"

df_pseudobulk_results_smode <- read.delim(file_pseudobulk_results_smode, stringsAsFactors = F)
colnames(df_pseudobulk_results_smode)[1] <- "SampleID"

batch_correction <- c("none")

dirOut <- paste0(dirOutRoot, "compare_to_pseudobulk/", "rmbatch_", batch_correction, "/")
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T, showWarnings = F)
}

df_pseudobulk_results <- eval(sym(paste0("df_pseudobulk_results_", batch_correction)))

# ------------
# Correlate Identified vs actual cell fraction in each Sample
# ------------
sampleIDs <- unique(dfActual$SampleID)
cell_type_order <- colnames(dfActual)[-1]

dfSample_CellFraction_detail <- NULL
dfSample_CellFraction_summary <- NULL

for(sampleID in sampleIDs)
{
  actual <- t(dfActual[dfActual$SampleID == sampleID, cell_type_order])
  result <- t(df_pseudobulk_results[df_pseudobulk_results$SampleID == sampleID, cell_type_order])

  dfCompDetailed <- data.frame(actual=actual, result=result)
  colnames(dfCompDetailed) <- c("actual", "result")
  dfCompDetailed$SampleID <- sampleID
  dfCompDetailed$celltype <- rownames(dfCompDetailed)
  corr <- cor(actual, result, method = "pearson")
  
  dfCompDetailed$correlation <- rep(corr, nrow(dfCompDetailed))
  dfSample_CellFraction_detail <- rbind(dfSample_CellFraction_detail, dfCompDetailed)
  
  dfCompSummary <- data.frame(SampleID=sampleID, correlation=corr[1])
  dfSample_CellFraction_summary <- rbind(dfSample_CellFraction_summary, dfCompSummary)
  
  print(sampleID)
  print(paste0("Correlation", " = ", corr))
  print(dfCompDetailed)
}

# Save Tables
fOut <- paste0(dirOut, "sample_cellfraction_summary.txt")
write.table(dfSample_CellFraction_summary, file = fOut, row.names = F, sep = "\t")

fOut <- paste0(dirOut, "sample_cellfraction_detailed.txt")
write.table(dfSample_CellFraction_detail, file = fOut, row.names = F, sep = "\t")

# Plot

# barplot (summary)
dfPlot <- dfSample_CellFraction_summary
gg <- ggplot(dfPlot, aes(x = reorder(SampleID, -correlation), y = correlation)) +
  xlab("SampleID") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5))
#gg

fOut <- paste0(dirOut, "sample_cellfraction_correlation_summary.barplot.pdf")
ggsave(filename = fOut, plot = gg)

# boxplot (summary)
dfPlot <- dfSample_CellFraction_summary
gg <- ggplot(dfPlot, aes(x = "SampleIDs", y = correlation)) +
  geom_boxplot() +
  xlab("Samples") +    
  ylab("Correlation (per sample)") +  
#  geom_jitter(width = 0.1) +
  theme(axis.text.x = element_blank())
#gg

fOut <- paste0(dirOut, "sample_cellfraction_correlation_summary.boxplot.pdf")
ggsave(filename = fOut, plot = gg, width = 5, height = 10, units = "cm")

# ---------------------
# detailed
# ---------------------
# stacked barplot
dfPlot <- dfSample_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))
gg <- ggplot(dfPlot, aes(x = SampleID, fill = celltype)) +
  geom_bar(aes(y = value), stat="identity", position = "stack", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  xlab("SampleID") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 10)) +
  facet_grid(rows = vars(variable))
#gg

fOut <- paste0(dirOut, "comparison_between_actual_and_predicted_celltypes_for_each_sample.stacked_barplot.pdf")
ggsave(filename = fOut, plot = gg, width = 28, height = 17, units = "cm")

#
dfPlot <- dfSample_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))
gg <- ggplot(dfPlot, aes(x = variable, fill = celltype)) +
  geom_bar(aes(y = value), stat="identity", position = "stack", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  #xlab("CIBERSORTx comparison") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 7.2)) +
  facet_grid(cols = vars(SampleID))
#gg

fOut <- paste0(dirOut, "comparison_between_actual_and_predicted_celltypes_for_each_sample.facet_SampleID.stacked_barplot.pdf")
ggsave(filename = fOut, plot = gg, width = 48, height = 10, units = "cm")

# dodged barplot
dfPlot <- dfSample_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))
gg <- ggplot(dfPlot, aes(x = celltype, fill = celltype)) +
  geom_bar(aes(y = value), stat="identity", position = "dodge", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  xlab("cell-type") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  #theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   # major
  theme(axis.text.x = element_text(size = 6, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   # minor  
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 10)) +
  facet_wrap(SampleID~variable, nrow = 4) +
  guides(fill = FALSE)
#gg

fOut <- paste0(dirOut, "comparison_between_actual_and_predicted_celltypes_for_each_sample.facet_SampleID.dodged_barplot.pdf")
if(celltype_classification == "celltype_major") {
  ggsave(filename = fOut, plot = gg, width = 50, height = 25, units = "cm")
} else if(celltype_classification == "celltype_minor") {
  ggsave(filename = fOut, plot = gg, width = 50, height = 25, units = "cm")
} else if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 100, height = 30, units = "cm")
}

# ------------
# Correlate Identified vs actual cell fraction for each cell-type
# ------------
rownames(dfActual) <- dfActual$SampleID
dfActual_t <- dfActual[-1]
dfActual_t <- as.data.frame(t(dfActual_t))
dfActual_t[, "cell_type"] <- rownames(dfActual_t)

rownames(df_pseudobulk_results) <- df_pseudobulk_results$SampleID
df_pseudobulk_results_t <- df_pseudobulk_results[-1]
df_pseudobulk_results_t <- as.data.frame(t(df_pseudobulk_results_t))
df_pseudobulk_results_t[, "cell_type"] <- rownames(df_pseudobulk_results_t)

cell_types <- rownames(dfActual_t)
sample_order <- colnames(dfActual_t)[-length(colnames(dfActual_t))]

dfCellType_CellFraction_detail <- NULL
dfCellType_CellFraction_summary <- NULL
for(cell_type in cell_types)
{
  actual <- t(dfActual_t[dfActual_t$cell_type == cell_type, sample_order])
  result <- t(df_pseudobulk_results_t[df_pseudobulk_results_t$cell_type == cell_type, sample_order])
  
  dfCompDetailed <- data.frame(actual=actual, result=result)
  colnames(dfCompDetailed) <- c("actual", "result")
  dfCompDetailed$difference <- dfCompDetailed$result - dfCompDetailed$actual
  dfCompDetailed$celltype <- cell_type
  dfCompDetailed$SampleID <- rownames(dfCompDetailed) 
  corr <- cor(actual, result, method = "pearson")
  
  dfCompDetailed$correlation <- rep(corr, nrow(dfCompDetailed))
  dfCellType_CellFraction_detail <- rbind(dfCellType_CellFraction_detail, dfCompDetailed)
  
  dfCompSummary <- data.frame(celltype=cell_type, correlation=corr[1])
  dfCellType_CellFraction_summary <- rbind(dfCellType_CellFraction_summary, dfCompSummary)
  
  print(cell_type)
  print(paste0("Correlation", " = ", corr))
  print(dfCompDetailed)  
}

# save tables
fOut <- paste0(dirOut, "celltype_cellfraction_summary.txt")
write.table(dfCellType_CellFraction_summary, file = fOut, row.names = F, sep = "\t")

fOut <- paste0(dirOut, "celltype_cellfraction_detailed.txt")
write.table(dfCellType_CellFraction_detail, file = fOut, row.names = F, sep = "\t")

# Plot

# barplot (summary)
dfPlot <- dfCellType_CellFraction_summary
gg <- ggplot(dfPlot, aes(x = reorder(celltype, -correlation), y = correlation, fill = celltype)) +
  xlab("Cell-type") +
  ylim(c(NA,1)) +
  geom_bar(stat="identity", colour = "black") +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = FALSE)  
#gg

fOut <- paste0(dirOut, "celltype_cellfraction_correlation_summary.barplot.pdf")
ggsave(filename = fOut, plot = gg, width = 15, height = 15, units = "cm")

# boxplot (summary)
dfPlot <- dfCellType_CellFraction_summary
gg <- ggplot(dfPlot, aes(x = "Cell-types", y = correlation)) +
  geom_boxplot() +
  ylab("Correlation (per cell-type)") +  
  xlab("Cell-type") +
  ylim(c(NA,1)) +
  #geom_jitter(width = 0.1) +
  theme(axis.text.x = element_blank())
#gg

fOut <- paste0(dirOut, "celltype_cellfraction_correlation_summary.boxplot.pdf")
ggsave(filename = fOut, plot = gg, width = 5, height = 10, units = "cm")

# ---------------------
# detailed
# ---------------------
# boxplot comparing cell-type proportions between actual and predicted
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))

gg <- ggplot(dfPlot, aes(x = celltype, y = value, fill = variable)) +
  geom_boxplot() +
  xlab("cell-type") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold"))
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.boxplot.pdf")
if(celltype_classification == "celltype_major") {
  ggsave(filename = fOut, plot = gg, width = 20, height = 18, units = "cm")
} else if(celltype_classification == "celltype_minor") {
  ggsave(filename = fOut, plot = gg, width = 25, height = 18, units = "cm")
} else if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 30, height = 18, units = "cm")
}

#
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))

gg <- ggplot(dfPlot, aes(x = variable, fill = SampleID)) +
  geom_bar(aes(y = value), stat="identity", position = "stack", colour="black") +
#  scale_y_continuous(labels=scales::percent) +  
#  xlab("CIBERSORTx comparison") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 7.5)) +
  facet_grid(cols = vars(celltype))
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.stacked_barplot.pdf")
if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 50, height = 18, units = "cm")  
} else {
  ggsave(filename = fOut, plot = gg, width = 30, height = 18, units = "cm")
}

# dodged barplot
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))

gg <- ggplot(dfPlot, aes(x = SampleID, fill = SampleID)) +
  geom_bar(aes(y = value), stat="identity", position = "dodge", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  xlab("SampleID") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(strip.text = element_text(face = "bold", size = 10)) +
  facet_grid(celltype~variable) + 
  guides(fill = FALSE)
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.facet_celltype_vs_variable.dodged_barplot.pdf")
if(celltype_classification == "celltype_major") {
  ggsave(filename = fOut, plot = gg, width = 20, height = 28, units = "cm")
} else if(celltype_classification == "celltype_minor") {
  ggsave(filename = fOut, plot = gg, width = 20, height = 50, units = "cm")
} else if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 20, height = 130, units = "cm", limitsize = FALSE)
}

#
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("actual", "result"), id.vars = c("SampleID", "celltype"))

gg <- ggplot(dfPlot, aes(x = variable, fill = celltype)) +
  geom_bar(aes(y = value), stat="identity", position = "dodge", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  xlab("cell-type") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  theme(strip.text.x = element_text(face = "bold", size = 7)) +
  #theme(strip.text.y = element_text(face = "bold", size = 10)) + # major 
  theme(strip.text.y = element_text(face = "bold", size = 6)) + # minor   
  facet_grid(celltype~SampleID) +
  guides(fill = FALSE)
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.facet_celltype_vs_SampleID.dodged_barplot.pdf")
if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 40, height = 90, units = "cm")
} else {
  ggsave(filename = fOut, plot = gg, width = 40, height = 40, units = "cm")
}

# -------------------------------
# Difference in proportions
# Predicted - Actual 
#    => +ve means over predicted
#    => -ve means under-predicted
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("difference"), id.vars = c("SampleID", "celltype"))

gg <- ggplot(dfPlot, aes(x = SampleID, fill = celltype)) +
  geom_bar(aes(y = value), stat="identity", position = "dodge", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  xlab("cell-type") +
  ylab("Difference in cell type proportion (Predicted - Actual)") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  #theme(strip.text = element_text(face = "bold", size = 10)) + # major
  theme(strip.text = element_text(face = "bold", size = 6)) + # minor  
  facet_grid(celltype~.) + 
  guides(fill = FALSE)  
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.difference_of_proportions.facet_celltype.dodged_barplot.pdf")
if(celltype_classification == "celltype_major") {
  ggsave(filename = fOut, plot = gg, width = 30, height = 28, units = "cm")
} else if(celltype_classification == "celltype_minor") {
  ggsave(filename = fOut, plot = gg, width = 30, height = 50, units = "cm")
} else if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 30, height = 90, units = "cm") 
}

# Predicted - Actual 
#    => +ve means over predicted
#    => -ve means under-predicted
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("difference"), id.vars = c("SampleID", "celltype"))

gg <- ggplot(dfPlot, aes(x = SampleID, fill = celltype)) +
  geom_bar(aes(y = value), stat="identity", position = "dodge", colour="black") +
  scale_y_continuous(labels=scales::percent) +  
  xlab("cell-type") +
  ylab("Difference in cell type proportion (Predicted - Actual)") +
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 15, face = "bold")) +
  #theme(strip.text = element_text(face = "bold", size = 10)) + # major
  theme(strip.text = element_text(face = "bold", size = 6)) + # minor  
  facet_grid(.~celltype) + 
  guides(fill = FALSE)  
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.difference_of_proportions.facet_sample.dodged_barplot.pdf")
if(celltype_classification == "celltype_major") {
  ggsave(filename = fOut, plot = gg, width = 30, height = 28, units = "cm") 
} else if(celltype_classification == "celltype_minor") {
  ggsave(filename = fOut, plot = gg, width = 70, height = 20, units = "cm")
} else if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 100, height = 20, units = "cm", limitsize = F)
}

#
dfPlot <- dfCellType_CellFraction_detail
dfPlot <- reshape2::melt(dfPlot, measure.vars = c("difference"), id.vars = c("SampleID", "celltype"))
colnames(dfPlot) <- c("SampleID", "celltype", "variable", "difference")

gg <- ggplot(dfPlot, aes(x = celltype, y = difference, fill = celltype)) +
  geom_boxplot() +
  scale_y_continuous(labels=scales::percent) +  
  xlab("cell-type") +
  ylab("Difference in cell type proportion (Predicted - Actual)") +
  theme(axis.text = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +   
  theme(axis.title = element_text(size = 10, face = "bold")) +
  guides(fill = FALSE)  
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.difference_of_proportions_across_celltypes.boxplot.pdf")
if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 40, height = 15, units = "cm")
} else {
  ggsave(filename = fOut, plot = gg, width = 20, height = 15, units = "cm")
}


# ===========================
# METABRIC
# ===========================

# clinical info (e.g., molecular subtypes)
dfMETABRICClinicalInfo <- read.delim(file_metabric_PAM50_info, stringsAsFactors = FALSE, comment.char = "#")
colnames(dfMETABRICClinicalInfo)[2] <- "SUBTYPE"

# ===========================
# METABRIC (combined)
# ===========================

# PARSE:
# ------
df_metabric_results <- NULL

# Discovery cohort results (Smode batch correction)
df_metabric <- read.delim(file_metabric_discovery_results_smode, stringsAsFactors = F)
colnames(df_metabric)[1] <- "SampleID"

df_metabric[, "batch_correction"] <- "smode"
df_metabric[, "COHORT"] <- "METABRIC_discovery"

df_metabric_results <- bind_rows(df_metabric_results, df_metabric) 

# Validation cohort results (Smode batch correction)
df_metabric <- read.delim(file_metabric_validation_results_smode, stringsAsFactors = F)
colnames(df_metabric)[1] <- "SampleID"

df_metabric[, "batch_correction"] <- "smode"
df_metabric[, "COHORT"] <- "METABRIC_validation"

df_metabric_results <- bind_rows(df_metabric_results, df_metabric) 

df_metabric_results <- dplyr::left_join(df_metabric_results, 
                                        dfMETABRICClinicalInfo, 
                                        by = c("SampleID" = "METABRIC_ID"))

# Save Table
dirOut <- paste0(dirOutRoot, "metabric/combined/")
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T, showWarnings = F)
}
fOut <- paste0(dirOut, "results.txt")
write.table(df_metabric_results, file = fOut, row.names = F, sep = "\t")

batch_correction_method <- c("smode")

message(paste0("batch_correction: ", batch_correction_method))
dirOut <- paste0(dirOutRoot, "metabric/combined/", "rmbatch_", batch_correction_method, "/")
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T, showWarnings = F)
}

df_metabric <- df_metabric_results %>% filter(batch_correction == batch_correction_method)

# ===================================
# plot distributions of cell types
# ===================================
# --------------
# Across all samples
# --------------
dfPlot <- df_metabric
dfPlot <- reshape2::melt(dfPlot, measure.vars = c(cell_type_order), id.vars = c("SampleID"))

gg <- ggplot(dfPlot, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  xlab("cell-type") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(size = 15, face = "bold")) +
  guides(fill = FALSE)
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.boxplot.pdf")
if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 30, height = 18, units = "cm")
} else {
  ggsave(filename = fOut, plot = gg, width = 20, height = 18, units = "cm")
}

# --------------
# Stratify by molecular subtype
# --------------
dfPlot <- df_metabric
dfPlot <- reshape2::melt(dfPlot, measure.vars = c(cell_type_order), id.vars = c("SampleID", "SUBTYPE"))

gg <- ggplot(dfPlot, aes(x = variable, y = value, fill = SUBTYPE)) +
  geom_boxplot() +
  xlab("cell-type") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(size = 15, face = "bold"))
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.split_by_molecular_subtype.boxplot.pdf")
if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 50, height = 18, units = "cm")
} else {
  ggsave(filename = fOut, plot = gg, width = 35, height = 18, units = "cm")
}

# Analyse each cell-type separately, stratify by molecular subtype
# and calculate signifcance of difference between subtypes for each cell-type
cell_types <- cell_type_order
for(cell_type in cell_types)
{
  dfPlot <- df_metabric
  dfPlot <- reshape2::melt(dfPlot, measure.vars = c(cell_type_order), id.vars = c("SampleID", "SUBTYPE"))
  
  dfPlot <- dfPlot %>% filter(variable == cell_type)
  
  gg <- ggplot(dfPlot, aes(x = SUBTYPE, y = value, fill = SUBTYPE)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(size = 0.1, position = position_jitterdodge(jitter.width = 0.15)) +
    xlab(cell_type) +
    ylab("Proportion of cell type") +
    theme(axis.text = element_text(size = 15)) +
    theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_text(size = 15, face = "bold")) +
    geom_signif(comparisons = list(c("Bas", "Her2"),
                                   c("Bas", "LumA"),
                                   c("Bas", "LumB"),
                                   c("Bas", "Norm"),
                                   c("Her2", "LumA"),
                                   c("Her2", "LumB"),
                                   c("Her2", "Norm"),
                                   c("LumA", "LumB"),
                                   c("LumA", "Norm"),
                                   c("LumB", "Norm")),
                test = "wilcox.test", step_increase = 0.05,
                textsize = 3, size = 0.4, tip_length = 0.01)
  #              y_position = c(0.5, 0.45, 0.4), textsize = 3, size = 0.4, tip_length = 0.02)
  #gg
  
  fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.split_by_molecular_subtype.", cell_type,".boxplot.pdf")
  ggsave(filename = fOut, plot = gg, width = 20, height = 25, units = "cm")
}

#
dfPlot <- df_metabric
dfPlot <- reshape2::melt(dfPlot, measure.vars = c(cell_type_order), id.vars = c("SampleID", "SUBTYPE"))

gg <- ggplot(dfPlot, aes(x = variable, y = value, fill = SUBTYPE)) +
  geom_boxplot() +
  xlab("cell-type") +
  ylab("Proportion of cell type") +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(size = 15, face = "bold")) +
  facet_wrap(~SUBTYPE) +
  guides(fill = FALSE)
#gg

fOut <- paste0(dirOut, "celltype_cellfraction.comparison_of_proportions.facet_by_molecular_subtype.boxplot.pdf")
if(celltype_classification == "celltype_subset") {
  ggsave(filename = fOut, plot = gg, width = 40, height = 20, units = "cm")
} else {
  ggsave(filename = fOut, plot = gg, width = 30, height = 18, units = "cm")
}


# ==============================
# H-clustering
# ==============================
dfColourCodes <- data.frame(name = c("Basal", "Her2", "LumA", "LumB", "Normal", "normal",
                                     "BRCA_Basal", "BRCA_Her2", "BRCA_LumA", "BRCA_LumB", "BRCA_Normal",
                                     "Bas", "Norm", "NC",
                                     "LumB_Her2", "HER2_ER", "ER", "TNBC", "HER2", "metaplastic",
                                     "Her2_Amp", "Her2_Non_Amp", "TBC",
                                     "LuminalA", "LuminalB", "LuminalB_HER2",
                                     "Normal_breast",
                                     "LumA_SC", "LumB_SC", "Her2_SC", "Basal_SC",
                                     "cancer", "unassigned", "normal",
                                     "METABRIC_discovery", "METABRIC_validation"),
                            value = c("red", "pink", "darkblue", "cyan", "green", "green",
                                      "red", "pink", "darkblue", "cyan", "green",
                                      "red", "green", "grey50",
                                      "pink", "pink", "blue", "red", "purple", "yellow",
                                      "black", "lightgrey", "yellow",
                                      "darkblue", "cyan", "pink",
                                      "grey50",
                                      "darkblue", "cyan", "pink", "red",
                                      "black", "grey50", "green",
                                      "#7BD18C", "#FFB0BD")
)


dirOut <- paste0(dirOutRoot, "metabric/combined/", "rmbatch_", batch_correction_method, "/clustering/")
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T, showWarnings = F)
}

for(row_scale in c("none", "zscore"))
{
  dfPlot <- df_metabric_results %>% filter(batch_correction == batch_correction_method)
  
  rownames(dfPlot) <- dfPlot$SampleID
  matPlot <- dfPlot[, cell_type_order]
  matPlot <- t(matPlot)
  if(row_scale == "zscore") {
    matPlot <- scale_rows_zscore(matPlot)
  }
  
  cluster_method <- "average"
  cluster_distance <- "euclidean"

  meta_data_cols <- rev(c("SUBTYPE", "COHORT"))
  dfColAnno <- dfPlot[, meta_data_cols, drop = F]
  
  ht_clustered <- complexHeatmap_wrapper(matData = matPlot,
                                         dfColAnno = dfColAnno,
                                         dfAnnotationColourCodes = dfColourCodes,
                                         min_max = c(-2, 2),
                                         scale_rows = FALSE,
                                         name = "Cell-type proportion\n(z-score)",
                                         show_column_names = F,
                                         cluster_columns = T,
                                         clustering_distance_columns = cluster_distance,
                                         clustering_method_columns = cluster_method,
                                         cluster_rows = T,
                                         clustering_distance_rows = cluster_distance,
                                         clustering_method_rows = cluster_method,
                                         row_names_gp = gpar(fontsize = 10),
                                         row_names_side = "right"
  )
  dirPlot <- paste0(dirOut, "hclustering/")
  if(!dir.exists(dirPlot)) {
    dir.create(dirPlot, showWarnings = F, recursive = T)
  }
  fHeatmap <- paste0(dirPlot, "cell_fractions_heatmap.hclust.",
                     "cluster_distance_", cluster_distance, "_",
                     "cluster_method_", cluster_method, ".",
                     "rowscale_", row_scale,
                     ".pdf")
  fExt <- tools::file_ext(fHeatmap)
  
  height <- 15
  width <- 50
  if(fExt == "png") {
    png(file = fHeatmap, height = height, width = width, units = "cm", res = 600)
  } else if(fExt == "pdf") {
    pdf(file = fHeatmap, height = height*0.39, width = width*0.39)
  } else {
    stop(paste0("Unsupported file type for heatmap: ", fExt))
  }
  message(paste0("start: ", Sys.time()))
  ht_draw <- draw(ht_clustered)
  message(paste0("end: ", Sys.time()))
  dev.off()
  
  # ---------------------
  # consensus clustering
  # ---------------------
  if(run_consensus_clustering)
  {
    if(row_scale == "zscore") 
    {
      message("RUNNING CONSENSUS CLUSTERING")
      # -------------------
      dirPlot <- paste0(dirOut, "consensus/")
      if(!dir.exists(dirPlot)) {
        dir.create(dirPlot, showWarnings = F, recursive = T)
      }
      
      max_k <- 10
      nCores <- 4
      
      matSimilarity <- matPlot
      set.seed(123)
      cp <- cola::consensus_partition(matSimilarity, top_n = nrow(matSimilarity), mc.cores = nCores, max_k = max_k, partition_repeat = 100, sample_by = "column", scale_rows = F)
      best_k <- suggest_best_k(cp)
      message(paste0("BEST k = ", best_k))
      dir_cola_report <- dirPlot
      cola_report(cp, dir_cola_report))
      
      message("SAVING COLA RESULTS")
      file_cola_results <- paste0(dir_cola_report, "cola_results.rds")
      saveRDS(cp, file = file_cola_results)
      
      message("EXTRACTING CLUSTERS")
      message("SAMPLES")    
      df_sample_clusters <- cbind(get_classes(cp, k = best_k), get_membership(cp, k = best_k))
      df_sample_clusters[, "SAMPLE_ID"] <- rownames(df_sample_clusters)
      colnames(df_sample_clusters)[which(colnames(df_sample_clusters) == "class")] <- "cluster"
      
      # group by clusterID
      df_sample_clusters <- df_sample_clusters %>% arrange(cluster)
      
      df_sample_clusters <- dplyr::left_join(df_sample_clusters, df_metabric[, c("SampleID", "COHORT")], by = c("SAMPLE_ID" = "SampleID"))
      
      file_sample_clusters <- paste0(dir_cola_report, "sample_clusters", ".all", ".txt")
      write.table(x = df_sample_clusters, file = file_sample_clusters, sep = "\t", row.names = FALSE)
      
      # Also write tables for the separate cohorts
      for(c in unique(df_sample_clusters$COHORT)) 
      {
        df_cohort <- df_sample_clusters %>% filter(COHORT == c)
        file_sample_clusters <- paste0(dir_cola_report, "sample_clusters", ".", c, ".txt")
        write.table(df_cohort, file_sample_clusters, sep = "\t", row.names = F)      
      }
      
      # NOTE: skip this if celltype=major as it doesn't work due to low numbers of cell-types/rows
      message("CELL_TYPES")
      if(celltype_classification == "celltype_major")
      {
        row_split <- NULL
        row_k <- 1
        matPlot_ordered <- matPlot[, df_sample_clusters$SAMPLE_ID]
        
      } else {
        df_celltype_clusters <- get_signatures(cp, k = best_k)
        df_celltype_clusters[, "celltype"] <- rownames(matPlot)[df_celltype_clusters$which_row]
        colnames(df_celltype_clusters)[which(colnames(df_celltype_clusters) == "km")] <- "cluster"
        
        # group by clusterID
        df_celltype_clusters <- df_celltype_clusters %>% arrange(cluster)
        file_celltype_clusters <- paste0(dir_cola_report, "celltype_clusters.txt")
        write.table(x = df_celltype_clusters, file = file_celltype_clusters, sep = "\t", row.names = FALSE)
        
        row_split <- df_celltype_clusters$cluster
        row_k <- length(unique(df_celltype_clusters$cluster))
        matPlot_ordered <- matPlot[df_celltype_clusters$celltype, df_sample_clusters$SAMPLE_ID]
      }
      
      # consensus heatmap
      # columns
      dfColAnno_ordered <- dfColAnno[df_sample_clusters$SAMPLE_ID, , drop = F]
      dfColAnno_ordered$SAMPLE_ID <- rownames(dfColAnno_ordered)
      dfColAnno_ordered <- inner_join(dfColAnno_ordered, df_sample_clusters[, c("SAMPLE_ID", "cluster")], by = "SAMPLE_ID")
      rownames(dfColAnno_ordered) <- dfColAnno_ordered$SAMPLE_ID
      dfColAnno_ordered <- dfColAnno_ordered[df_sample_clusters$SAMPLE_ID, c("COHORT", "SUBTYPE", "cluster"), drop = F]
      # add colours  
      set.seed(123)      
      dfColourCodes <- rbind(dfColourCodes, data.frame(name = as.factor(unique(df_sample_clusters$cluster)), value = distinctColorPalette(best_k)))      
      
      set.seed(123)
      ht_consensus <- complexHeatmap_wrapper(matData = matPlot_ordered, 
                                             dfColAnno = dfColAnno_ordered, 
                                             dfAnnotationColourCodes = dfColourCodes, 
                                             min_max = c(-2, 2), 
                                             scale_rows = FALSE, 
                                             name = "Cell-type proportion\n(z-score)",
                                             show_column_names = F,
                                             cluster_columns = T,
                                             clustering_distance_columns = cluster_distance,
                                             clustering_method_columns = cluster_method,
                                             cluster_rows = T,
                                             clustering_distance_rows = cluster_distance,
                                             clustering_method_rows = cluster_method,                             
                                             row_names_gp = gpar(fontsize = 10), 
                                             row_names_side = "right",
                                             column_split=dfColAnno_ordered$cluster,
                                             row_split=row_split
      )
      
      fHeatmap <- paste0(dirPlot, "cell_fractions_heatmap.consensus_clustering.",
                         "column_k_", best_k, ".",
                         "row_k_", row_k,                   
                         "cluster_distance_", cluster_distance, "_",
                         "cluster_method_", cluster_method, ".",
                         "rowscale_", row_scale,
                         ".pdf")
      fExt <- tools::file_ext(fHeatmap)
      
      height <- 15
      width <- 50
      if(fExt == "png") {
        png(file = fHeatmap, height = height, width = width, units = "cm", res = 600)
      } else if(fExt == "pdf") {
        pdf(file = fHeatmap, height = height*0.39, width = width*0.39)  
      } else {
        stop(paste0("Unsupported file type for heatmap: ", fExt))
      }
      message(paste0("start: ", Sys.time()))
      set.seed(123)      
      ht_draw <- draw(ht_consensus)
      message(paste0("end: ", Sys.time()))
      dev.off()
    }
  }
  
}


# ---------------------------
message("R Session Info")
sessionInfo()

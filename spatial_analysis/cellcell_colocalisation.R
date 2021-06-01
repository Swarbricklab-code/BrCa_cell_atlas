# SPATIAL CELL-CELL CO-LOCALISATION ANALYSIS
# Sunny Z. Wu
# 
#
# Run in screen 
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_spatial
# R
# 
#
# 01: SETUP---------------------------------------

# setup
library(zeallot)
library(STutility)
library(ggplot2)
library(magrittr)
library(Seurat)
library(magick)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# for stats
library(rstatix)


# DIRECTORY
dir.create("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/cellcell_colocalisation/")
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/cellcell_colocalisation/")
dir.create("output")

# cbind.fill function 
require(plyr) # requires plyr for rbind.fill()
cbind.fill <- function(...) {                                                                                                                                                       
  transpoted <- lapply(list(...),t)                                                                                                                                                 
  transpoted_dataframe <- lapply(transpoted, as.data.frame)                                                                                                                         
  return (data.frame(t(rbind.fill(transpoted_dataframe))))                                                                                                                          
}

# correlation p-value function
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

# 02: LOAD DATA ---------------------------------------------------------------------

se.list <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/RDATA_Visium_brca_objects_stereoscope.Rdata")

# load stereoscope cluster IDs per annotation tier
temp_sample_id <- c(2,5,7:10)
temp_celltypes_CTP_minor <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_CTP_minor.Rdata")
temp_celltypes_minor <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_minor.Rdata")
temp_celltypes_subset <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_subset.Rdata")
temp_celltypes_major  <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_major.Rdata")

# 03: LOAD LUNDEBERG HER2+ DATASET --------------------------------------------

# load her2+ data per cell tier
for(tier in c("major", "minor","subset")){
  
  matrix_datalist <- readRDS(paste0("/share/ScratchGeneral/sunwu/projects/spatial/data/Andersson_etal_data/RDATA/RDATA_Lundeberg_etal_2021_stereosceop_tier_",tier,".Rdata"))
  
  m <- paste0("matrix_datalist_",tier)
  assign(m, matrix_datalist)

}


# 04: COMPUTE CORRELATION MATRIX FOR HER2+ DATASET ----------------------------------------------

for(tier in c("major", "minor","subset")){
  print(tier)
  
  matrix_datalist <- get(paste0("matrix_datalist_",tier))
  corMat.list <- NULL
  
  for(sample in names(matrix_datalist)){
    print(sample)
    corMat <- cor(matrix_datalist[[sample]])
    
    corMat.list[[sample]] <- corMat
  }
  
  n <- paste0("corMat.list_",tier)
  assign(n, corMat.list)
  
}

# 05: COMBINED CORRELATION MATRICES TNBC, ER & HER2 DATASETS  --------------------------------------------

# clusters of interest 
temp_clustersofinterest <- c("CAFs_MSC_iCAF_like", "CAFs_myCAF_like", "Epithelial_cancer","Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE")
temp_clustersofinterest <- c(temp_clustersofinterest,
                             grep("T_cells_",temp_celltypes_subset,value=T))
temp_clustersofinterest <- c(temp_clustersofinterest,
                             grep("B_cells_",temp_celltypes_subset,value=T))
temp_clustersofinterest <- temp_clustersofinterest[!temp_clustersofinterest %in% c("T_cells_c9_NK_cells_AREG", "T_cells_c10_NKT_cells_FCGR3A", "T_cells_c11_MKI67")]
temp_clustersofinterest <- c(temp_clustersofinterest, "T_cells_CD8+", "T_cells_CD4+")
temp_clustersofinterest <- c(temp_clustersofinterest,
                             grep("Endothelial_",temp_celltypes_minor,value=T))
temp_clustersofinterest <- c(temp_clustersofinterest,
                             grep("PVL_",temp_celltypes_minor,value=T))
temp_clustersofinterest <- c(temp_clustersofinterest, "PVL", "Endothelial", "T_cells", "B_cells")
temp_clustersofinterest <- c(temp_clustersofinterest, "Myeloid_c9_Macrophage_2_CXCL10")
temp_clustersofinterest <- temp_clustersofinterest[!temp_clustersofinterest %in% c("Endothelial_Lymphatic_LYVE1")]
temp_clustersofinterest <- temp_clustersofinterest[!temp_clustersofinterest %in% c("PVL_Cycling")]

# annotation data frame
temp_dataframe <- data.frame(row.names = temp_clustersofinterest)
temp_dataframe$lineage <- as.vector(temp_clustersofinterest)
for(major_lineage in c("B_cells", "CAFs", "T_cells", "Endothelial", "PVL", "Myeloid")){
  
  major_lineage_subsets <- grep(major_lineage, unique(temp_dataframe$lineage), value=T)
  
  for(subset in major_lineage_subsets){
    temp_dataframe$lineage[temp_dataframe$lineage == subset] <- major_lineage
  }
}

# combined with HER2+ lundeberg dataset
temp_sampleIDs <- NULL
for(sample in temp_sample_id){
  temp_sampleID <- paste0(unique(se.list[[sample]]@meta.data$patientid))
  print(temp_sampleID)
  
  # subset level
    temp_df_metagenes <- se.list[[sample]]@reductions$stereoscope_tnbcsubset@cell.embeddings
    temp_colnames <- gsub("/","_",temp_celltypes_subset)
    temp_colnames <- gsub("-","_",temp_colnames)
    
    # fix typo in cell cluster ID
    temp_colnames <- gsub("Differeniated","Differentiated",temp_colnames)
    # temp_colnames <- gsub("\\+","",temp_colnames)
    colnames(temp_df_metagenes) <- temp_colnames
  
  # minor level
  if(sample %in% 7:10){
    temp_vector <- get(paste0("temp_celltypes_CTP_minor"))
  } else {
    temp_vector <- get(paste0("temp_celltypes_minor"))
  }
  temp_vector <- gsub("-","_",temp_vector)
  
  temp_df_metagenes_minor <- se.list[[sample]]@reductions$stereoscope_tnbcminor@cell.embeddings
  colnames(temp_df_metagenes_minor) <- temp_vector
  
  # major level
  temp_df_metagenes_major <- se.list[[sample]]@reductions$stereoscope_tnbcmajor@cell.embeddings
  colnames(temp_df_metagenes_major) <- temp_celltypes_major
  
  temp_df_metagenes <- cbind(temp_df_metagenes, temp_df_metagenes_minor, temp_df_metagenes_major)
  
  temp_df_metagenes <- temp_df_metagenes[,colnames(temp_df_metagenes) %in% temp_clustersofinterest]
  
  temp_df_metagenes <- temp_df_metagenes[, !duplicated(colnames(temp_df_metagenes)),drop=F]
  
  corMat <- cor(temp_df_metagenes)
  colnames(corMat) <- gsub("\\+","",colnames(corMat))
  rownames(corMat) <- gsub("\\+","",rownames(corMat))
  
  temp_df_m <- reshape2::melt(corMat)
  
  rownames(temp_df_m) <- paste0(temp_df_m$Var1, "_vs_", temp_df_m$Var2)
  
  if(sample == 2){
    temp_dataframe_combined <- temp_df_m
    temp_dataframe_combined <- temp_dataframe_combined[,2,drop=F]
  }
  
  temp_df_m <- temp_df_m[,3,drop=F]
  temp_df_m <- temp_df_m[order(rownames(temp_df_m)),,drop=F]
  
  colnames(temp_df_m) <-  temp_sampleID
  
  if(sample == 2){
    temp_cormat_combined <- temp_df_m
  } else (
    temp_cormat_combined <- cbind(temp_cormat_combined, temp_df_m)
  )
  
  m <- paste0("temp_proportions_",temp_sampleID)
  assign(m, temp_df_metagenes)
  temp_sampleIDs <- c(temp_sampleIDs, temp_sampleID)
  
}

for(her2dataset in c(1:length(corMat.list_subset))){
  
  temp_sampleID <- paste0("Lundeberg_etal_",names(corMat.list_subset)[her2dataset]
  )
  print(temp_sampleID)
  temp_df_metagenes <- matrix_datalist_subset[[her2dataset]]
  temp_df_metagenes_minor <- matrix_datalist_minor[[her2dataset]]
  temp_df_metagenes_major <- matrix_datalist_major[[her2dataset]]
  temp_df_metagenes <- cbind(temp_df_metagenes,temp_df_metagenes_minor, temp_df_metagenes_major)
  temp_df_metagenes <- temp_df_metagenes[,colnames(temp_df_metagenes) %in% temp_clustersofinterest]
  temp_df_metagenes <- temp_df_metagenes[,! colnames(temp_df_metagenes) %in% c("B_cells_Memory.1", "B_cells_Naive.1", "Endothelial_RGS5.1", "Endothelial_ACKR1.1", "Endothelial_CXCL12.1")]
  corMat <- cor(temp_df_metagenes)
  colnames(corMat) <- gsub("\\+","",colnames(corMat))
  rownames(corMat) <- gsub("\\+","",rownames(corMat))
  
  temp_df_m <- reshape2::melt(corMat)
  rownames(temp_df_m) <- paste0(temp_df_m$Var1, "_vs_", temp_df_m$Var2)
  temp_df_m <- temp_df_m[,3,drop=F]
  temp_df_m <- temp_df_m[order(rownames(temp_df_m)),,drop=F]
  
  colnames(temp_df_m) <-  temp_sampleID
  
  temp_colnames <- c(colnames(temp_cormat_combined), colnames(temp_df_m))
  temp_cormat_combined <- cbind.fill(temp_cormat_combined, 
                                     temp_df_m)
  colnames(temp_cormat_combined) <-  temp_colnames
  
  m <- paste0("temp_proportions_",temp_sampleID)
  assign(m, temp_df_metagenes)
  temp_sampleIDs <- c(temp_sampleIDs, temp_sampleID)
  
}


# 06: PLOT CORRELATION HEATMAP -----------------------------------------------

# subtype annotation df
temp_dataframe_col <- data.frame(row.names = colnames(temp_cormat_combined),
                                 subtype = c("TNBC", "TNBC", "TNBC", "TNBC",
                                             "ER+", "ER+",
                                             "HER2+", "HER2+", "HER2+", "HER2+","HER2+","HER2+","HER2+"))
# clusters of interest to plot
temp_clusters_to_plot <- c("CAFs_MSC_iCAF_like_vs_CAFs_myCAF_like",
                           "CAFs_MSC_iCAF_like_vs_B_cells_Memory",
                           "CAFs_MSC_iCAF_like_vs_B_cells_Naive",
                           "CAFs_MSC_iCAF_like_vs_T_cells_CD4",
                           "CAFs_MSC_iCAF_like_vs_T_cells_CD8",
                           "CAFs_myCAF_like_vs_B_cells_Memory",
                           "CAFs_myCAF_like_vs_B_cells_Naive",
                           "CAFs_myCAF_like_vs_T_cells_CD4",
                           "CAFs_myCAF_like_vs_T_cells_CD8",
                           "PVL_vs_Endothelial",
                           "PVL_Immature_vs_PVL_Differentiated",
                           "PVL_Immature_vs_Endothelial_RGS5",
                           "PVL_Differentiated_vs_Endothelial_ACKR1",
                           "PVL_Differentiated_vs_Endothelial_CXCL12",
                           "Myeloid_c1_LAM1_FABP5_vs_Myeloid_c2_LAM2_APOE",
                           "Myeloid_c1_LAM1_FABP5_vs_T_cells_CD8",
                           "Myeloid_c1_LAM1_FABP5_vs_T_cells_CD4",
                           "Myeloid_c2_LAM2_APOE_vs_T_cells_CD8",
                           "Myeloid_c2_LAM2_APOE_vs_T_cells_CD4",
                           "Myeloid_c9_Macrophage_2_CXCL10_vs_T_cells_CD8",
                           "Myeloid_c9_Macrophage_2_CXCL10_vs_T_cells_CD4"
)

temp_cormat_combined_subset <- temp_cormat_combined[rownames(temp_cormat_combined) %in% temp_clusters_to_plot,,drop=F]

# cell type annotation dataframe
temp_dataframe <- data.frame(row.names = temp_clusters_to_plot)
temp_dataframe$lineage1 <- c(rep("CAFs",9), rep("PVL cells", 5), rep("Myeloid", 7))
temp_dataframe$lineage2 <- c("CAFs",
                             "B_cells",
                             "B_cells",
                             "T_cells",
                             "T_cells",
                             "B_cells",
                             "B_cells",
                             "T_cells",
                             "T_cells",
                             "Endothelial",
                             "PVL cells",
                             "Endothelial",
                             "Endothelial",
                             "Endothelial",
                             "Myeloid",
                             "T_cells",
                             "T_cells",
                             "T_cells",
                             "T_cells",
                             "T_cells",
                             "T_cells"
)
# DF for multiple testing filtering below
temp_df_significance_cortest <- temp_cormat_combined_subset

# fix up "/"s from rownames to not break plotting
temp_cormat_combined_subset <- temp_cormat_combined_subset[temp_clusters_to_plot,,drop=F]
rownames(temp_cormat_combined_subset) <- gsub("CAFs_","", rownames(temp_cormat_combined_subset))
rownames(temp_cormat_combined_subset) <- gsub("MSC_iCAF_like","MSC/iCAF-like", rownames(temp_cormat_combined_subset))
rownames(temp_cormat_combined_subset) <- gsub("myCAF_like","myCAF-like", rownames(temp_cormat_combined_subset))
rownames(temp_dataframe) <- rownames(temp_cormat_combined_subset) 

# colour pallete
temp_colours <- brewer.pal(6, "Set3")

# plot with pheatmap
pheatmap::pheatmap(temp_cormat_combined_subset,
                   filename = paste0("output/cellcell_colocalisation_figure.pdf"),
                   cluster_rows = F,
                   cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                   gaps_row = c(1,5,9, 14),
                   gaps_col = c(4,6),
                   breaks = seq(-0.5, 0.5, length.out = 101),
                   legend_breaks = seq(-0.5,0.5,0.1),
                   cellheight = 12.5,
                   cellwidth = 12.5,
                   annotation_col = temp_dataframe_col,
                   annotation_row = temp_dataframe,
                   annotation_colors = list(
                     lineage1 = c("CAFs" = temp_colours[1], 
                                  "PVL cells" = temp_colours[2],
                                  "Myeloid" = temp_colours[3]
                     ),
                     lineage2 = c("CAFs" = temp_colours[1], 
                                  "PVL cells" = temp_colours[2],
                                  "Myeloid" = temp_colours[3],
                                  "Endothelial" = temp_colours[4],
                                  "T_cells" = temp_colours[5],
                                  "B_cells" = temp_colours[6]
                     ),
                     subtype = c("TNBC" = "red",
                                 "ER+" = "blue",
                                 "HER2+" = "pink")
                   ),
                   na_col = "grey")


# 07: STATISTICAL SIGNIFICANCE AND MULTIPLE COMPARISONS TESTING ---------------

# clusters of interest to plot
temp_clusters_to_plot <- c("CAFs_MSC_iCAF_like", "CAFs_myCAF_like",
                           "B_cells_Memory",
                           "B_cells_Naive",
                           "T_cells_CD4+",
                           "T_cells_CD8+",
                           "PVL", "Endothelial",
                           "PVL_Immature","PVL_Differentiated",
                           "Endothelial_RGS5",
                           "Endothelial_ACKR1",
                           "Endothelial_CXCL12",
                           "Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE",
                           "Myeloid_c9_Macrophage_2_CXCL10"
)


# compute significance of Pearson correlations
pval_combined <- NULL
for(sample in temp_sampleIDs){
  
  temp_df_proportions <- get(paste0("temp_proportions_",sample))
  temp_df_proportions <- temp_df_proportions[,colnames(temp_df_proportions) %in% temp_clusters_to_plot,drop=F]
  corMat_test <- cor.test.p(temp_df_proportions)

  # remove duplicates prior to p-val adjustment
  temp_num <- unique(colnames(corMat_test))[1]
  corMat_test_m <- NULL
  for(celltype in unique(colnames(corMat_test))){
    
    corMat_test_subset <- corMat_test[,celltype,drop=F]
    corMat_test_subset <- reshape2::melt(corMat_test_subset)
    corMat_test_subset$p_val_adjust =
      p.adjust(corMat_test_subset$value,
               method = "BH")
    
    for(num in temp_num){
      corMat_test_subset <- corMat_test_subset[! corMat_test_subset$Var1 == num,,drop=F]
    }
    
    corMat_test_m <- rbind(corMat_test_m, corMat_test_subset)
    temp_num <- c(temp_num,celltype)
  }
  
  corMat_test_m$patientID <- sample
  pval_combined <- rbind(pval_combined, corMat_test_m)
}

# remove duplicate correlations
pval_combined$ID <- paste0(pval_combined$Var1,"_vs_",pval_combined$Var2)
for(num in temp_num){
  pval_combined <- pval_combined[! pval_combined$ID == paste0(num,"_",num),,drop=F]
}

# write to file
rownames(pval_combined) <- c(1: nrow(pval_combined))

pval_combined$significant <- NA
for(row in c(1:nrow(pval_combined))){
  if(pval_combined[row,"p_val_adjust"] < 0.05 ){
    pval_combined[row,"significant"] <- "Y"
  } else {
    pval_combined[row,"significant"] <- "N"
  }
}

write.csv(pval_combined, "pearson_corr_cellcell_colocalisation.csv")


# 08: PLOT CORRELATION HEATMAP WITHOUT NON-SIGNIFICANT -----------------------------------------------

# all non-significant value indicated by grey
temp_df_significance_cortest_filtered <- temp_df_significance_cortest
pval_combined$ID2 <- paste0(pval_combined$Var2,"_vs_",pval_combined$Var1)
pval_combined$ID <- gsub("\\+","",pval_combined$ID)
pval_combined$ID2 <- gsub("\\+","",pval_combined$ID2)

for(celltypes in unique(rownames(temp_df_significance_cortest_filtered))){
    print(celltypes)
    
    temp_pval_df <- NULL
    if(celltypes %in% pval_combined$ID){
      temp_pval_df <- pval_combined[pval_combined$ID == celltypes,,drop=F]
    }
    if(celltypes %in% pval_combined$ID2){
      temp_pval_df <- pval_combined[pval_combined$ID2 == celltypes,,drop=F]
    }
  
    for(sample in colnames(temp_df_significance_cortest_filtered)){
      # print(sample)
      temp_pval <- temp_pval_df[temp_pval_df$patientID == sample,"p_val_adjust",drop=T]
      
      if(temp_pval > 0.05){
        temp_df_significance_cortest_filtered[celltypes,sample] <- NA
      }
    }
    
    }

# order 
temp_clusters_to_plot <- c("CAFs_MSC_iCAF_like_vs_CAFs_myCAF_like",
                           "CAFs_MSC_iCAF_like_vs_B_cells_Memory",
                           "CAFs_MSC_iCAF_like_vs_B_cells_Naive",
                           "CAFs_MSC_iCAF_like_vs_T_cells_CD4",
                           "CAFs_MSC_iCAF_like_vs_T_cells_CD8",
                           "CAFs_myCAF_like_vs_B_cells_Memory",
                           "CAFs_myCAF_like_vs_B_cells_Naive",
                           "CAFs_myCAF_like_vs_T_cells_CD4",
                           "CAFs_myCAF_like_vs_T_cells_CD8",
                           "PVL_vs_Endothelial",
                           "PVL_Immature_vs_PVL_Differentiated",
                           "PVL_Immature_vs_Endothelial_RGS5",
                           "PVL_Differentiated_vs_Endothelial_ACKR1",
                           "PVL_Differentiated_vs_Endothelial_CXCL12",
                           "Myeloid_c1_LAM1_FABP5_vs_Myeloid_c2_LAM2_APOE",
                           "Myeloid_c1_LAM1_FABP5_vs_T_cells_CD8",
                           "Myeloid_c1_LAM1_FABP5_vs_T_cells_CD4",
                           "Myeloid_c2_LAM2_APOE_vs_T_cells_CD8",
                           "Myeloid_c2_LAM2_APOE_vs_T_cells_CD4",
                           "Myeloid_c9_Macrophage_2_CXCL10_vs_T_cells_CD8",
                           "Myeloid_c9_Macrophage_2_CXCL10_vs_T_cells_CD4"
)


# plot with pheatmap
# fix up "/"s from rownames to not break plotting
temp_df_significance_cortest_filtered <- temp_df_significance_cortest_filtered[temp_clusters_to_plot,,drop=F]
rownames(temp_df_significance_cortest_filtered) <- gsub("CAFs_","", rownames(temp_df_significance_cortest_filtered))
rownames(temp_df_significance_cortest_filtered) <- gsub("MSC_iCAF_like","MSC/iCAF-like", rownames(temp_df_significance_cortest_filtered))
rownames(temp_df_significance_cortest_filtered) <- gsub("myCAF_like","myCAF-like", rownames(temp_df_significance_cortest_filtered))
rownames(temp_dataframe) <- rownames(temp_df_significance_cortest_filtered) 

pheatmap::pheatmap(temp_df_significance_cortest_filtered,
                   filename = paste0("output/cellcell_colocalisation_figure_filtered_ns_pval_adjust.pdf"),
                   cluster_rows = F,
                   cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
                   gaps_row = c(1,5,9, 14),
                   gaps_col = c(4,6),
                   breaks = seq(-0.5, 0.5, length.out = 101),
                   legend_breaks = seq(-0.5,0.5,0.1),
                   cellheight = 12.5,
                   cellwidth = 12.5,
                   annotation_col = temp_dataframe_col,
                   annotation_row = temp_dataframe,
                   annotation_colors = list(
                     lineage1 = c("CAFs" = temp_colours[1], 
                                  "PVL cells" = temp_colours[2],
                                  "Myeloid" = temp_colours[3]
                     ),
                     lineage2 = c("CAFs" = temp_colours[1], 
                                  "PVL cells" = temp_colours[2],
                                  "Myeloid" = temp_colours[3],
                                  "Endothelial" = temp_colours[4],
                                  "T_cells" = temp_colours[5],
                                  "B_cells" = temp_colours[6]
                     ),
                     subtype = c("TNBC" = "red",
                                 "ER+" = "blue",
                                 "HER2+" = "pink")
                   ),
                   na_col = "grey")



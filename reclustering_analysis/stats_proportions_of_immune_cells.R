# PROPORTIONS OF IMMUNE CELLS
#
# Sunny Wu
# Run in screen 
# qrsh -pe smp 4 -l mem_requested=10G -P TumourProgression
# source activate r_seuratdev
# R
#
# 01: SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(rstatix)
library(reshape2)
library(cowplot)

# DIRECTORIES
dir.create("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/05_figures/03_immune_figures_v4_pval_adjusted/")
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/05_figures/03_immune_figures_v4_pval_adjusted/")

# 02: LOAD DATA ---------------------------------------------------------------

for(temp_cell_type in c("T_cells", "Myeloid")) {
  
  seurat_10X_integrated <-
    readRDS(paste0("/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/Jul2020_updated_objectfile/Rdata_final/02_subset_object_",
                   temp_cell_type,
                   ".Rdata"))
  
  n <- paste0("seurat_10X_integrated_",temp_cell_type)
  assign(n, seurat_10X_integrated)
}



print(unique(seurat_10X_integrated_T_cells@meta.data$celltype_subset))
print(unique(seurat_10X_integrated_Myeloid@meta.data$celltype_subset))

# temp_csv <- read.csv("/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/Jul2020_updated_objectfile/celltype_ids_newlabels.csv")

# 03: PROPORTIONS DF ------------------------------------------------------

# PLOT ALL CLUSTERS STATS
temp_colname <- "celltype_subset"
temp_cellprop_df_all_combined <- NULL
for(temp_cell_type in c("T_cells", "Myeloid")) {
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",temp_cell_type))
  Idents(seurat_10X_integrated) <- temp_colname
  seurat_10X_integrated@meta.data$subtype <- factor(seurat_10X_integrated@meta.data$subtype,
                                                    levels = rev(c("ER+", "HER2+", "TNBC")))
  
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  
  # estimate proportions
  temp_cluster_names <- 
    unique(seurat_10X_integrated@meta.data[,temp_colname])
  
  temp_Names <- 
    as.vector(unique(seurat_10X_integrated@meta.data$orig.ident))
  
  temp_cellprop_df_all <- NULL
  for(i in c(1:length(temp_Names))) {
    # print(i)
    print(temp_Names[i])
    temp_filtered_summary <- 
      subset(seurat_10X_integrated@meta.data,
             orig.ident == temp_Names[i])
    
    temp_molecular_subtype <- unique(temp_filtered_summary$subtype)
    
    temp_cellprop_df <-
      data.frame(unclass(table(temp_filtered_summary[,temp_colname])))
    colnames(temp_cellprop_df) <- "Freq"
    
    # append missing clusters
    temp_missing_clusters <- temp_cluster_names[! temp_cluster_names %in% unique(rownames(temp_cellprop_df))]
    if(length(temp_missing_clusters) > 0){
      print("no missing clusters")
      temp_missing_clusters <- data.frame(row.names = temp_missing_clusters,
                                          rep(0,times=length(temp_missing_clusters)))
      colnames(temp_missing_clusters) <- "Freq"
      temp_cellprop_df <- rbind(temp_cellprop_df, temp_missing_clusters)
    }
    
    # identify proportions
    temp_cellprop_df$proportion <- NULL
    temp_cellprop_df$proportion <- temp_cellprop_df[,1]/sum(temp_cellprop_df[,1])
    
    temp_cellprop_df <- 
      data.frame(sample = rep(temp_Names[i],
                              times = length(row.names(temp_cellprop_df))),
                 cell_type = row.names(temp_cellprop_df),
                 value = temp_cellprop_df[,1],
                 proportion = temp_cellprop_df[,2],
                 subtype = temp_molecular_subtype[[1]]
      )
    
    temp_cellprop_df_all <- 
      rbind(temp_cellprop_df_all,
            temp_cellprop_df)
    
  }
  
  
  temp_cellprop_df_all$subtype <- factor(temp_cellprop_df_all$subtype,
                                         levels = rev(c("ER+", "HER2+", "TNBC")))
  
  my_comparisons <- list(c("ER+", "TNBC"), c("ER+", "HER2+"),
                         c("TNBC", "HER2+")
  )
  
  temp_nrow <- length(temp_cluster_names)
  
  # plot graph with p vals
  temp_gene_boxplot <- ggboxplot(temp_cellprop_df_all, 
                                 x = "subtype",
                                 y="proportion",
                                 color = "subtype",
                                 # legend = "none",
                                 facet.by = "cell_type", 
                                 repel = F,
                                 nrow = temp_nrow,
                                 scales ="free_x",
                                 palette = rev(c("blue", "pink", "red"))) +
    stat_compare_means(method = "t.test",
                       paired = F,
                       label = "p.signif", 
                       ref.group = NULL,
                       hide.ns = T,
                       comparisons = my_comparisons,
                       bracket.size = 0.1) +  # Pairwise comparison against all
    # stat_compare_means(method = "anova", 
    #                    label.y = temp_max_label_y) +
    xlab(" ") + ylab("Proportion") +
    # rotate_x_text(angle = 45) + 
    coord_flip()
  
  temp_gene_boxplot <- ggpar(temp_gene_boxplot,
                             legend.title="Subtype")
  
  pdf(file = paste0("05_PROPORTIONS/01_proportions_combined_",temp_cell_type, ".pdf"), width = 6, height = 15)
  print(temp_gene_boxplot)
  dev.off()
  
  temp_cellprop_df_all_combined <- rbind(temp_cellprop_df_all_combined, temp_cellprop_df_all)
  
  # n <- paste0("temp_cellprop_df_all_",temp_cell_type)
  # assign(n, temp_cellprop_df_all)
  
}


# 04: STATISTICAL SIGNIFICANCE AND MULTIPLE TESTING ------------------------------------------------------------

# T-CELLS
# temp_cellprop_df_all_combined_Tcells <- temp_cellprop_df_all_combined[temp_cellprop_df_all_combined$cell_type %in% unique(grep("T_cells",temp_cellprop_df_all_combined$cell_type, value=T)),]

temp_cell_ids_of_interest <- c("T_cells_c6_IFIT1", "T_cells_c11_MKI67", "T_cells_c8_CD8+_LAG3")

temp_cellprop_df_all_sig <- temp_cellprop_df_all_combined[temp_cellprop_df_all_combined$cell_type %in% temp_cell_ids_of_interest,]

stat.test <- temp_cellprop_df_all_sig %>%
    group_by(cell_type) %>%
    t_test(proportion ~ subtype, p.adjust.method = "BH")
  stat.test <- as.data.frame(stat.test)
  stat.test$pvalmethod <- "BH"

write.csv(stat.test, "t_test_Tcells.csv")


# MYELOID
temp_cellprop_df_all_combined_Myeloid <- temp_cellprop_df_all_combined[temp_cellprop_df_all_combined$cell_type %in% unique(grep("Myeloid",temp_cellprop_df_all_combined$cell_type, value=T)),]

stat.test <- temp_cellprop_df_all_combined_Myeloid %>%
    group_by(cell_type) %>%
    t_test(proportion ~ subtype, p.adjust.method =  "BH")
  stat.test <- as.data.frame(stat.test)
  stat.test$pvalmethod <-  "BH"
write.csv(stat.test, "t_test_Myeloid.csv")


# 05: PLOT PROPORTIONS ----------------------------------------------------------

# PLOT BOXPLOTS
temp_colname <- "celltype_subset"
temp_cellprop_df_all_combined <- NULL
my_comparisons <- list(c("ER+", "TNBC"), c("ER+", "HER2+"),
                       c("TNBC", "HER2+")
)

## all clusters
for(temp_cell_type in c("T_cells", "Myeloid")) {
  temp_cellprop_df_all <-  temp_cellprop_df_all_combined[temp_cellprop_df_all_combined$cell_type %in% unique(grep(temp_cell_type,temp_cellprop_df_all_combined$cell_type, value=T)),]
  
  temp_nrow <- length(unique(temp_cellprop_df_all$cell_type))
  
  # plot graph with p vals
  temp_gene_boxplot <- ggboxplot(temp_cellprop_df_all, 
                                 x = "subtype",
                                 y="proportion",
                                 color = "subtype",
                                 # legend = "none",
                                 facet.by = "cell_type", 
                                 repel = F,
                                 nrow = temp_nrow,
                                 scales ="free_x",
                                 palette = rev(c("blue", "pink", "red"))) +
    stat_compare_means(method = "t.test",
                       paired = F,
                       label = "p.signif", 
                       ref.group = NULL,
                       hide.ns = T,
                       comparisons = my_comparisons,
                       bracket.size = 0.1) +  # Pairwise comparison against all
    # stat_compare_means(method = "anova", 
    #                    label.y = temp_max_label_y) +
    xlab(" ") + ylab("Proportion") +
    # rotate_x_text(angle = 45) + 
    coord_flip()
  
  temp_gene_boxplot <- ggpar(temp_gene_boxplot,
                             legend.title="Subtype")
  
  pdf(file = paste0("01_proportions_combined_",temp_cell_type, ".pdf"), width = 6, height = 15)
  print(temp_gene_boxplot)
  dev.off()
  
  temp_cellprop_df_all_combined <- rbind(temp_cellprop_df_all_combined, temp_cellprop_df_all)
  
  # n <- paste0("temp_cellprop_df_all_",temp_cell_type)
  # assign(n, temp_cellprop_df_all)
  
}


## clusters of interest
### T-cells
temp_cell_ids_of_interest <- c("T_cells_c6_IFIT1", "T_cells_c11_MKI67", "T_cells_c8_CD8+_LAG3")

temp_cellprop_df_all_sig <- temp_cellprop_df_all_combined[temp_cellprop_df_all_combined$cell_type %in% temp_cell_ids_of_interest,]

temp_gene_boxplot <- ggboxplot(temp_cellprop_df_all_sig, 
                               x = "subtype",
                               y="proportion",
                               color = "subtype",
                               legend = "none",
                               facet.by = "cell_type", 
                               repel = F,
                               nrow = 1,
                               scales ="free_x", 
                               palette = rev(c("blue", "pink", "red"))) +
  stat_compare_means(method = "t.test",
                     paired = F,
                     label = "p.signif", 
                     ref.group = NULL,
                     hide.ns = T,
                     comparisons = my_comparisons,
                     bracket.size = 0.1) +  # Pairwise comparison against all
  xlab(" ") + ylab("Proportion") +
  coord_flip()

temp_gene_boxplot <- ggpar(temp_gene_boxplot,
                           legend.title="Subtype")

pdf(file = paste0("02_proportions_figure_Tcells.pdf"), width = 7, height = 2)
print(temp_gene_boxplot)
dev.off()


### Myeloid cells
temp_cell_ids_of_interest <- c("Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE")
temp_cellprop_df_all_sig <- temp_cellprop_df_all_combined[temp_cellprop_df_all_combined$cell_type %in% temp_cell_ids_of_interest,]

temp_gene_boxplot <- ggboxplot(temp_cellprop_df_all_sig, 
                               x = "subtype",
                               y="proportion",
                               color = "subtype",
                               legend = "none",
                               facet.by = "cell_type", 
                               repel = F,
                               nrow = 2,
                               scales ="free_x", 
                               palette = rev(c("blue", "pink", "red"))) +
  stat_compare_means(method = "t.test",
                     paired = F,
                     label = "p.signif", 
                     ref.group = NULL,
                     hide.ns = T,
                     comparisons = my_comparisons,
                     bracket.size = 0.1) +  # Pairwise comparison against all
  xlab(" ") + ylab("Proportion") +
  coord_flip()

temp_gene_boxplot <- ggpar(temp_gene_boxplot,
                           legend.title="Subtype")

pdf(file = paste0("02_proportions_figure_Myeloid.pdf"), width = 3, height = 3)
print(temp_gene_boxplot)
dev.off()




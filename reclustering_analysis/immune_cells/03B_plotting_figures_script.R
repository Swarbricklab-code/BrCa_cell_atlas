# Celltype default reclustering script 
# 20191013
#
#
# 01 PARSE ARGUMENTS  ------------------------------------------------------------

# arguments from command line
temp_start_time <- Sys.time()
print(temp_start_time)
temp_args <-
  commandArgs(trailingOnly = T)

# 01 CELL TYPE
temp_cell_type <- 
  temp_args[1]

# 02 # CORES/WORKERS
temp_run_parallel_cores <- 
  temp_args[2]
temp_run_parallel_cores <- 
  as.numeric(temp_run_parallel_cores)
temp_run_parallel_cores <-
  temp_run_parallel_cores/2

# 04 MEM PER CORE (DIVIDIED BY 2 TO AVOID CRASHED IN JACKSTRAW)
temp_mem_per_core <- 
  temp_args[3]
temp_mem_per_core <-
  as.numeric(gsub("G", "", temp_mem_per_core))
temp_mem_per_core <-
  temp_mem_per_core/2

# 05 TEMP WD
temp_wd <- 
  temp_args[4]

# RUN PARALELL ?
temp_run_parallel <- T

if(temp_cell_type == "B_cells" | temp_cell_type == "Plasma_cells") {
  temp_PC_number <- 10
  temp_reduction_use <- "TSNE"
  temp_res <- 0.8
}
if(temp_cell_type == "T_cells") {
  temp_PC_number <- 30
  temp_reduction_use <- "UMAP"
  temp_res <- 1
}
if(temp_cell_type == "Myeloid_cells") {
  temp_PC_number <- 20
  temp_reduction_use <- "UMAP"
  temp_res <- 0.8
}

# temp_cell_type <- "Myeloid_cells"
# temp_run_parallel_cores <- 4
# temp_run_parallel_cores <-
#   temp_run_parallel_cores/2
# temp_mem_per_core <- 10
# temp_mem_per_core <-
#   temp_mem_per_core/2
# temp_run_parallel <- T
# temp_wd <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/run04/analysis_01_T_cells/"



# SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(future)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(grid)


# DIRECTORY ---------------------------------------------------------------

setwd(paste0(temp_wd))

dir.create("03_output_filtered_plotting")
setwd("03_output_filtered_plotting")

# LOAD DATA ---------------------------------------------------------------

for(type in c("integrated")) {
  
  seurat_10X_integrated <- readRDS(paste0("../02_output_filtered/RDATA_02_PROCESSED_FILTERED_RESCALED", type, ".Rdata"))
  
  # CLUSTER HIGHER RES FOR CERTAIN CELLTYPES
  if(temp_cell_type == "T_cells"){
    
      seurat_10X_integrated <-
        FindClusters(
          object = seurat_10X_integrated,
          graph.name = paste0("SUBSET", temp_PC_number), 
          resolution = temp_res
        )

  }

  # remove small clusters
  temp_colname <- paste0("SUBSET", temp_PC_number, "_res.",temp_res)
  Idents(seurat_10X_integrated) <- temp_colname
  
  # remove sampleIDs with < 30 cells
  print(table(seurat_10X_integrated@meta.data[,temp_colname]))
  temp_df <- as.data.frame(table(seurat_10X_integrated@meta.data[,temp_colname]))
  temp_df_filtered <- temp_df[temp_df$Freq > 5,]
  temp_df_filtered <- as.vector(temp_df_filtered$Var1)
  temp_removed <- as.vector(temp_df$Var1[!temp_df$Var1 %in% temp_df_filtered])
  
  if(length(temp_removed) > 0){
    print("trimming clusters")
    Idents(seurat_10X_integrated) <- temp_colname
    seurat_10X_integrated <- SubsetData(seurat_10X_integrated,
                              ident.use=temp_df_filtered)
    print(temp_removed)
  }
  
  # factorise clusterIDs
  print(unique(seurat_10X_integrated@meta.data[,temp_colname]))
  seurat_10X_integrated@meta.data[,temp_colname] <- 
    factor(seurat_10X_integrated@meta.data[,temp_colname],
           levels=str_sort(unique(seurat_10X_integrated@meta.data[,temp_colname]),numeric = T))
  print(unique(seurat_10X_integrated@meta.data[,temp_colname]))
  
  # assign
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
}

# DIMPLOTS ----------------------------------------------------------------

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 8,
      height = 8,
      res = 300,
      units = 'in'
    )
  }


for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  
  # DIMPLOT
  temp_png_function(paste0("01_FILTERED_UMAP_C_", type, ".png"))
  temp_dimplot <- DimPlot(
    object = seurat_10X_integrated,
    label.size = 8,
    pt.size = 0.5,
    reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
    label = T,
    group.by = paste0("SUBSET", temp_PC_number, "_res.",temp_res)
  )
  print(temp_dimplot)
  dev.off()
  
  # FEATUREPLOT RNA
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 10,
        height = 10,
        res = 300,
        units = 'in'
      )
    }
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  
  if(temp_cell_type == "B_cells" | temp_cell_type == "Plasma_cells") {
    temp_markers <-  c("CD19", "MS4A1", "BLK", "FCRL2", "TCL1A", "JCHAIN", "IGKC", "CD79A", "TNFRSF17")
  }
  if(temp_cell_type == "T_cells") {
    temp_markers <-  c("CD3D", "ITGAE", "CD4", "CD8A", "CD8B", "MKI67", "FOXP3", "NCR1", "XCL1", "KLRK1", "LAG3", "PDCD1", "CD274", "TOX", "CD86", "IL2RA")
  }
  if(temp_cell_type == "Myeloid_cells") {
    temp_markers <-  c("CD68","ITGAM","CD74", "LILRB4", "FCER1A", "CD163", "IFITM3", "CLEC9A", "S100A12", "RGCC", "PDCD1", "CD274", "CD80", "CD86",  "CD84")
    # 
    # temp_markers_2 <-  c("C5AR1","FCAR","HLA-DQA1","CD81","CD14","FLT3","CADM1","IL3RA","CD1C","IRF8","IRF4","THBD","CD5","FCGR3A")
    
    #cDC1
    temp_markers_1 <- c("THBD", "ICOSLG", "XCR1", "FLT3", "CD226", "LY75")
     # CD141 = THBD ICOSL = ICOSLG CD32=FCGR2A CD135 = FLT3 CD205 = LY75
    
    #cDC2
    temp_markers_2 <- c("CD5", #DC2 defining markers & top pre-DC marker vs pDC 
                        "CD1D", "CD1C", "CD2", "BTLA", "FCER1A", "C5AR1", "FCAR" 
                        ) # CD88 = C5AR1, CD89 = FCAR FCER1A
    #cDC3
    temp_markers_3 <- c("CD14", "CD163", "PVRL2", "CD36", "FCGR1A", "ITGAM", "MRC1" 
    ) # CD112 = NECTIN2, CD64 = FCGR1A, CD206 = MRC1 NECTIN2=PVRL2
    
    # mono_macro 
    temp_markers_4 <- c("CD14",  "FCGR1A", "ITGAM", "MRC1", "ITGAL", "ITGAX", "NCAM1", "FCGR3A"
    ) # CD11b = ITGAM, CD11C = ITGAX, CD56 = NCAM1, CD16 = FCGR3A CD11A = ITGAL
     # het on pDCs
    temp_markers_5 <- c("MME", "ENTPD1", "ITGA6") # heterogenous on pDC
    # cd49f = ITGA6
    
    # preDCs
    temp_markers_6 <- c("CD22", "CDH1", "TNFSF13B", "IL3RA", "CD244", "HAVCR2") # heterogenous on pDC
    # CD324=CDH1 CD268=BAFF=TNFSF13B CD123 CD123=IL3RA CD366=HAVCR2
    
    # pDCs
    temp_markers_7 <- c("IL3RA", "CLEC4C", "MME", "ENTPD1", "ITGA6") # heterogenous on pDC
    # CD324=CDH1 CD268=BAFF CD123 CD123=IL3RA CD366=HAVCR2
    
    # misc
    temp_markers_8 <- c("CD34", "ITGA7", "CD40", "CD86", "CD80", "HLA-DRA", "PECAM1")
    # CD10 = MME, CD39 = ENTPD1, CD49f = ITGA6 , CD303 = CLEC4C
  }  
  
  temp_png_function(paste0("02_FILTERED_featureplots_01_",type, ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0.25,
    reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
    min.cutoff = 0,
    order = T)
  print(temp_featureplot)
  dev.off()
  
  temp_png_function(paste0("03_FILTERED_vlnplot_01_",type, ".png"))
  tempvlnplot <- VlnPlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0,
    group.by = paste0("SUBSET", temp_PC_number, "_res.",temp_res)
    )
  print(tempvlnplot)
  dev.off()
  
  if(temp_cell_type == "Myeloid_cells") {
    for(num in c(1:8)) {
    
    temp_markers_subset <- get(paste0("temp_markers_",num))
    
    temp_png_function(paste0("04_FILTERED_featureplots_02_",num, ".png"))
    temp_featureplot <- FeaturePlot(
      object = seurat_10X_integrated,
      features = temp_markers_subset,
      pt.size = 0.25,
      reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
      min.cutoff = 0,
      order = T)
    print(temp_featureplot)
    dev.off()
    }
  }
  
  temp_png_function_2 <-
    function(x) {
      png(
        file = (x),
        width = 12,
        height = 6,
        res = 300,
        units = 'in'
      )
    }
  
  for(metadata in c("clinical_subtype")) {
    
    temp_png_function_2(paste0("05_FILTERED_",type,"_split_by_", metadata, ".png"))
    temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      pt.size = 0.5,
      reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
      split.by = paste0(metadata),
      label = T,
      label.size = 8,
      group.by = paste0("SUBSET", temp_PC_number, "_res.",temp_res) 
    )
    print(temp_dimplot)
    dev.off()
  }
  
  temp_png_function_2 <-
    function(x) {
      png(
        file = (x),
        width = 15,
        height = 15,
        res = 300,
        units = 'in'
      )
    }
  
  for(metadata in c("orig.ident")) {
    
    temp_png_function_2(paste0("06_FILTERED_",type,"_split_by_", metadata, ".png"))
    temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      pt.size = 0.5,
      reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
      split.by = paste0(metadata),
      label = T,
      label.size = 4,
      group.by = paste0("SUBSET", temp_PC_number, "_res.0.8") 
    )
    print(temp_dimplot)
    dev.off()
  }
  
  
  
  
}


# PLOTTING STROMAL PROPORTIONS WITH T-TEST  ---------------------------------------

library(reshape2)
library(ggpubr)
library(cowplot)

# if(temp_cell_type == "T_cells" | temp_cell_type == "B_cells"){
#     temp_split_lymphocytes_by_all_params <- c(F,T)
# } else temp_split_lymphocytes_by_all_params <- c(F)

temp_split_lymphocytes_by_all_params <- c(F)

for(type in c("integrated")) {

  for(temp_split_lymphocytes_by_all in temp_split_lymphocytes_by_all_params) {
    
    seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
    # generate dataframe of stromal proportions

      if(temp_split_lymphocytes_by_all){
        if(temp_cell_type == "T_cells" | temp_cell_type == "B_cells"){
          if(temp_cell_type == "T_cells"){
            temp_dir <- "B_cells"
          }
          if(temp_cell_type == "B_cells"){
            temp_dir <- "T_cells"
          }

          seurat_10X_integrated_lymphocytes <- readRDS(paste0("../../analysis_01_",temp_dir,"/02_output_filtered/RDATA_02_PROCESSED_FILTERED_RESCALED", type, ".Rdata"))
      
          seurat_10X_integrated@meta.data$SUBSET20_res.0.8_lymphocyte <-
            seurat_10X_integrated@meta.data[,paste0("SUBSET", temp_PC_number, "_res.",temp_res)]
      
          seurat_10X_integrated_lymphocytes@meta.data$SUBSET20_res.0.8_lymphocyte <-
            "other_cluster"
      
          seurat_10X_integrated <- merge(seurat_10X_integrated,
                                         y = seurat_10X_integrated_lymphocytes
                                         )
          temp_cluster_names <-
            unique(seurat_10X_integrated@meta.data$SUBSET20_res.0.8_lymphocyte)
          temp_cluster_names <-
            temp_cluster_names[! temp_cluster_names %in% "other_cluster"]
      
        }
        if(temp_cell_type == "Myeloid_cells"){
          temp_cluster_names <- 
            unique(seurat_10X_integrated@meta.data[,paste0("SUBSET", temp_PC_number, "_res.",temp_res)])
        }
      } else {
        temp_cluster_names <- 
          unique(seurat_10X_integrated@meta.data[,paste0("SUBSET", temp_PC_number, "_res.",temp_res)])
      }
      
      temp_Names <- 
        unique(seurat_10X_integrated@meta.data$orig.ident)
      
      temp_cluster_names <- factor(temp_cluster_names,
                                   levels=(0:length(temp_cluster_names)))
      temp_cluster_names <- sort(temp_cluster_names)
      
      
      for(i in c(1:length(temp_Names))) {
        print(i)
        print(temp_Names[i])
        temp_filtered_summary <- 
          subset(seurat_10X_integrated@meta.data,
                 orig.ident == temp_Names[i])
        
        temp_molecular_subtype <- unique(temp_filtered_summary$clinical_subtype)
    
        if(temp_split_lymphocytes_by_all){
          
          if(temp_cell_type == "T_cells" | temp_cell_type == "B_cells"){
          temp_cellprop_df <-
            data.frame(unclass(table(temp_filtered_summary$SUBSET20_res.0.8_lymphocyte)))
          colnames(temp_cellprop_df) <- "Freq"
          }
          if(temp_cell_type == "Myeloid_cells"){
            temp_cellprop_df <-
              data.frame(unclass(table(temp_filtered_summary[,paste0("SUBSET", temp_PC_number, "_res.",temp_res)])))
            colnames(temp_cellprop_df) <- "Freq"
          }
        } else {
          temp_cellprop_df <-
            data.frame(unclass(table(temp_filtered_summary[,paste0("SUBSET", temp_PC_number, "_res.",temp_res)])))
          colnames(temp_cellprop_df) <- "Freq"
        }
        
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
        
        # temp_cellprop_df <- temp_cellprop_df[!rownames(temp_cellprop_df) %in% "other_cluster",]
        
        n <- paste0("temp_cellprop_df_",
                    i)
        
        assign(n,
               (temp_cellprop_df <- 
                  data.frame(sample = rep(temp_Names[i],
                                          times = length(row.names(temp_cellprop_df))),
                             cell_type = row.names(temp_cellprop_df),
                             value = temp_cellprop_df[,1],
                             proportion = temp_cellprop_df[,2],
                             subtype = temp_molecular_subtype[[1]]
                  )))
        
      }
      
      temp_cellprop_df_all <-
        NULL
      
      for(i in c(1:length(temp_Names))) {
        
        temp_cell_pop_df <- get(paste0("temp_cellprop_df",
                                       "_",
                                       i))
        
        temp_cellprop_df_all <- 
          rbind(temp_cellprop_df_all,
                temp_cell_pop_df)
      }
      
      # temp_cellprop_df_all$cell_type <- factor(temp_cellprop_df_all$cell_type,
      #                                          levels = c(1,3,6,4,2,0,7,5))
      
      temp_cellprop_df_all$subtype <- factor(temp_cellprop_df_all$subtype,
                                             levels = c("ER+", "HER2+", "TNBC"))
      
      ## plot
      # temp_png_function <-
      #   function(x) {
      #     png(
      #       file = (x), 
      #       width = 4, 
      #       height = 5, 
      #       res = 300, 
      #       units = 'in'
      #     )
      #   }
      
      temp_labels <- 
        NULL
      
      for(cluster in temp_cluster_names){
        print(cluster)
        temp_cellprop_df_all_subset <- 
          temp_cellprop_df_all[temp_cellprop_df_all$cell_type %in% cluster,]
        
        temp_max_label_y <- 
          max(temp_cellprop_df_all_subset$proportion)+0.2
        
        my_comparisons <- list( c("ER+", "TNBC"), c("ER+", "HER2+"), c("TNBC", "HER2+") )
        
        # plot graph with p vals
        temp_gene_boxplot <- ggboxplot(temp_cellprop_df_all_subset, 
                                       x = "subtype",
                                       y="proportion",
                                       color = "subtype") +
          # stat_compare_means(comparisons = my_comparisons,
          #                    method = "t.test")
          # stat_compare_means(method = "anova") +        # Add global annova p-value
          # stat_compare_means(label = "p.signif", 
          #                    method = "t.test",
          #                    ref.group = ".all.", 
          #                    hide.ns = F) +  # Pairwise comparison against all
          stat_compare_means(method = "t.test",
                             ref.group = ".all.",
                             hide.ns = F,
                             comparisons = my_comparisons) +  # Pairwise comparison against all
          stat_compare_means(method = "anova", 
                             label.y = temp_max_label_y) +
          # scale_color_manual(values= c("#66c2a5", "#fc8d62", "#8da0cb","#a6d854")) +
          # ylab(expression(CD8+~T-cells~per~1~mm^{2})) +
          xlab(" ") + ylab(" ")
        
        
        # temp_png_function(paste0("stromal_proportions/proportions_cluster_",cluster,".png"))
        # print(temp_gene_boxplot)
        # dev.off()
        
        n <- paste0("temp_gene_boxplot_",cluster)
        assign(n, temp_gene_boxplot)
        
        
        m <- paste0("Cluster ", cluster)
        
        temp_labels <- append(temp_labels, m)
        
      }
      
      temp_pdf_function <-
        function(x) {
          pdf(
            file = (x),
            width = 14,
            height = 12,
            useDingbats=F
          )
        }
      
      
      temp_grid <- 
        plot_grid(plotlist=mget(paste0("temp_gene_boxplot_",temp_cluster_names)),
                  labels = temp_labels,
                  nrow = (round(length(temp_cluster_names)/4)),
                  ncol = 4)
      
      temp_pdf_function(paste0("07_proportions_cluster_combined_splitby_lymphocyte_",temp_split_lymphocytes_by_all,".pdf"))
      print(temp_grid)
      dev.off()
      
      
      }
}


# SEEING GAMMYS ANNOTATIONS FROM CID3838 OVERLAID -------------------------

seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))

# subset 3838 from integrated object
temp_metadata <- seurat_10X_integrated@meta.data
temp_metadata_3838 <- temp_metadata[temp_metadata$orig.ident %in% "CID3838",]

# annotations from Gammy
temp_csv <- 
  read.csv("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/CID3838_metadata.csv")
rownames(temp_csv) <- temp_csv$X

rownames(x = temp_csv) <- gsub("3838", "CID3838", row.names(temp_csv))

temp_csv <- temp_csv[colnames(temp_csv) %in% c("broad_lineage",
                                               "cell_type_ADT")]

print("number of non-matching barcodes")
temp_matching_barcodes <- 
  length(rownames(temp_csv))-as.numeric(length(setdiff(rownames(temp_csv), rownames(temp_metadata_3838))))
print(temp_matching_barcodes)

temp_csv <- temp_csv[rownames(temp_csv) %in% rownames(temp_metadata_3838),]


if(temp_matching_barcodes > 0){
  
  # append back to original metadata
  temp_metadata <- temp_metadata[!rownames(temp_metadata) %in% rownames(temp_csv),]
  temp_metadata_df <- data.frame(row.names = rownames(temp_metadata))
  temp_metadata_df$broad_lineage <- NA
  temp_metadata_df$cell_type_ADT <- NA
  
  temp_metadata <- rbind(temp_metadata_df,
                         temp_csv)
  
  # sort and append back to original data
  temp_metadata_sorted <- 
    temp_metadata[rownames(seurat_10X_integrated@meta.data),,drop=F]
  
  all.equal(rownames(temp_metadata_sorted),rownames(seurat_10X_integrated@meta.data))
  
  seurat_10X_integrated <- AddMetaData(seurat_10X_integrated,
                                       temp_metadata_sorted)
  
  Idents(seurat_10X_integrated) <- "orig.ident"
  seurat_10X_integrated_subset <- subset(seurat_10X_integrated,
                                      idents="CID3838")
}

if(temp_matching_barcodes > 0){
  
# plot
    temp_png_function <-
      function(x) {
        png(
          file = (x),
          width = 8,
          height = 8,
          res = 300,
          units = 'in'
        )
      }
    
    temp_num <- 8
    for(colname in c("cell_type_ADT", "broad_lineage")){
      temp_png_function(paste0("0",temp_num,"_3838annotations_", colname, ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated_subset,
        label.size = 8,
        pt.size = 1,
        reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
        label = F,
        group.by = colname
      )
      print(temp_dimplot)
      dev.off()
      temp_num <- temp_num+1
    }
    
    
    temp_png_function <-
      function(x) {
        png(
          file = (x),
          width = 14,
          height = 14,
          res = 300,
          units = 'in'
        )
      }
    temp_png_function(paste0("0",temp_num,"_3838annotations_splitby_", "orig.ident", ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated_subset,
        pt.size = 2,
        reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
        split.by = "cell_type_ADT",
        label = F,
        label.size = 4,
        group.by = paste0("SUBSET", temp_PC_number, "_res.",temp_res)
      )
      print(temp_dimplot)
      dev.off()
      
      saveRDS(seurat_10X_integrated,
              "RDATA_seurat_objected_filtered_CID3838_gammy_annotated.Rdata")
}





# CLUSTER AVG GE ---------------------------------------------------------

  print(cluster)

  DefaultAssay(seurat_10X_integrated) <- "RNA"
  temp_colname <- paste0("SUBSET", temp_PC_number, "_res.",temp_res)
  Idents(seurat_10X_integrated) <- temp_colname
  
  temp_cluster.allmarkers <- 
    read.csv(paste0("../02_output_filtered/07_MAST_FINDALLMARKERS_BULK_res.","0.8","_RNA_assay.csv"))
  
  temp_cluster.allmarkers <-
    (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(25,
                                                             avg_logFC))
  
  
  print(table(temp_cluster.allmarkers$cluster))
  print(length(unique(temp_cluster.allmarkers$gene)))
  
  temp_cluster.averages <- 
    AverageExpression(seurat_10X_integrated,
                      return.seurat = TRUE,
                      verbose = T,
                      assays="RNA",
                      features=unique(temp_cluster.allmarkers$gene))
  
  temp_data_frame <- 
    as.matrix(temp_cluster.averages@assays$RNA@data)
  
  temp_data_frame_hclust <- temp_data_frame
  
  
# plot
library("viridis")  
hmcol <- viridis(8)

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 8,
      height = 16,
      res = 300,
      units = 'in'
    )
  }

  # temp_png_function(paste0("11_pheatmap.png"))
  # pheatmap(temp_data_frame, 
  #                       color = rev(hmcol),
  #                       cluster_cols = T, 
  #                       cluster_rows = T, 
  #                       scale = "row", 
  #                       clustering_distance_cols = "correlation",
  #                       clustering_distance_rows = "correlation",
  #                       fontsize_row = 4,
  #                       show_rownames = T,
  #                       show_colnames = T, 
  #                       fontsize_col = 20,
  #                       annotation_legend = T,
  #                       cellheight= 2, 
  #                       cellwidth = 20,
  #                       gaps_col = NULL, 
  #                       annotation_names_col = T, 
  #                       angle_col = 45,
  #                       treeheight_row = 75, 
  #                       treeheight_col = 50,
  #                       legend = T,
  #                       border_color=FALSE)
  # dev.off()


{
  if(temp_cell_type == "B_cells" | temp_cell_type == "Plasma_cells") {
    sigGenes_v <- c("CD19", "MS4A1", "BLK", "FCRL2", "TCL1A", "JCHAIN", "IGKC", "CD79A", "TNFRSF17")
  }
  if(temp_cell_type == "T_cells"){
    sigGenes_v <- c("CD3D", "ITGAE", "CD4", "CD8A", "CD8B", "MKI67", "FOXP3", "NCR1", "XCL1", "KLRK1", "LAG3", "PDCD1", "CD274", "TOX", "CD86", "IL2RA")
  }
  if(temp_cell_type == "Myeloid_cells"){
    sigGenes_v <- unique(c("CD68","ITGAM","CD74", "LILRB4", "FCER1A", "CD163", "IFITM3", "CLEC9A", "S100A12", "RGCC", "PDCD1", "CD274", "CD80", "CD86",  "CD84", "THBD", "ICOSLG", "XCR1", "FLT3", "CD226", "LY75", "CD5",
                    "CD1D", "CD1C", "CD2", "BTLA", "FCER1A", "C5AR1", "FCAR", 
                    "CD14", "CD163", "PVRL2", "CD36", "FCGR1A", "ITGAM", "MRC1",
                    "CD14",  "FCGR1A", "ITGAM", "MRC1", "ITGAL", "ITGAX", "NCAM1", "FCGR3A",
                    "MME", "ENTPD1", "ITGA6",
                    "CD22", "CDH1", "TNFSF13B", "IL3RA", "CD244", "HAVCR2",
                    "IL3RA", "CLEC4C", "MME", "ENTPD1", "ITGA6",
                    "CD34", "ITGA7", "CD40", "CD86", "CD80", "HLA-DRA", "PECAM1"
                    )
                    )
  }
  
  rowMeta_df <- data.frame(Sig = rep("No", length(rownames(temp_data_frame_hclust))), 
                           stringsAsFactors = F,
                           row.names = rownames(temp_data_frame_hclust))
  for (gene_v in sigGenes_v) rowMeta_df[rownames(rowMeta_df) == gene_v, "Sig"] <- gene_v
  
}

heat <- pheatmap(temp_data_frame_hclust, 
                 color = rev(hmcol),
                 cluster_cols = T, 
                 cluster_rows = T, 
                 scale = "row", 
                 clustering_distance_cols = "correlation", 
                 clustering_distance_rows = "correlation",
                 # annotation_row = rowMeta_df,
                 fontsize_row = 7.5,
                 show_rownames = T,
                 show_colnames = T, 
                 fontsize_col = 10,
                 # annotation_legend = T,
                 cellheight= 1, 
                 cellwidth = 20,
                 gaps_col = NULL, 
                 # annotation_names_col = T, 
                 angle_col = 45,
                 # treeheight_row = 75, 
                 # treeheight_col = 50,
                 legend = F,
                 # border_color=FALSE,
                 main = paste0(temp_cell_type," DEGs"))

### Get order of genes after clustering
genesInHeatOrder_v <- heat$tree_row$labels[heat$tree_row$order]
whichSigInHeatOrder_v <- which(genesInHeatOrder_v %in% sigGenes_v)

print("exists in data frame")
for(gene in sigGenes_v){
  print(grep(paste0(gene),
             genesInHeatOrder_v,
             value=T))
}

whichSigInHeatOrderLabels_v <- genesInHeatOrder_v[whichSigInHeatOrder_v]


sigY <- 1 - (1/length(rownames(temp_data_frame_hclust)) * whichSigInHeatOrder_v)

### Change title
whichMainGrob_v <- which(heat$gtable$layout$name == "main")
heat$gtable$grobs[[whichMainGrob_v]] <- textGrob(label = paste0(temp_cell_type," DEGs"), 
                                                 gp = gpar(fontsize = 16))

### Remove rows
whichRowGrob_v <- which(heat$gtable$layout$name == "row_names")
heat$gtable$grobs[[whichRowGrob_v]] <- textGrob(label = whichSigInHeatOrderLabels_v,
                                                y = sigY,
                                                vjust = 1)

#
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x - unit(0.5, "npc"),
                           x1 = new.label$x + unit(0.25, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.6, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

temp_png_function(paste0("11_pheatmap.png"))
add.flag(heat,
         kept.labels = sigGenes_v,
         repel.degree = 1)
dev.off()
print(sigGenes_v)


# SAVE OBJECTS -------------------------------------------------

seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))

saveRDS(seurat_10X_integrated, 
          paste0("RDATA_01_PROCESSED_FILTERED", type, ".Rdata"))

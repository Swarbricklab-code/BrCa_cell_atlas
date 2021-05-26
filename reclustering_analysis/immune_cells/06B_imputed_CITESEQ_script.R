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
# temp_wd <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/run03_fixedmyeloidprolif/analysis_01_Myeloid_cells/"



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

dir.create("06_imputed_CITESeq")
setwd("06_imputed_CITESeq")

# LOAD DATA ---------------------------------------------------------------

temp_colname <- paste0("SUBSET", temp_PC_number, "_res.",temp_res)

for(type in c("integrated")) {
  
  seurat_10X_integrated <- readRDS(paste0("../03_output_filtered_plotting/RDATA_01_PROCESSED_FILTEREDintegrated.Rdata"))
  
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
}

# Path to imputed seurat objects
temp_path_to_objects <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/07_CITE_imputation_Nenad/run02_191121/output/"

for(norm in c("norm", "raw")){
  
  load(paste0(temp_path_to_objects,
              "imputed_",temp_cell_type,"_",
              norm,
              ".Rdata"))
  
  seurat_10X_integrated <- 
    get(paste0("seurat_10X_integrated_",type,"_filtered"))
  
  if(temp_cell_type == "T_cells"){
    cells <- subset(cells,
                    cells = rownames(seurat_10X_integrated@meta.data))
  }
  seurat_10X_integrated <- subset(seurat_10X_integrated,
                                  cells = rownames(cells@meta.data))
  
  if(temp_cell_type == "T_cells"){
    temp_df <- cells@assays$ADT
    temp_df <- temp_df[,colnames(seurat_10X_integrated)]
    print(all.equal(colnames(seurat_10X_integrated),colnames(temp_df)))  

    seurat_10X_integrated[["ADT"]] <- 
      CreateAssayObject(data = temp_df)
    
    
  } else {
    print(all.equal(colnames(seurat_10X_integrated),colnames(cells)))  
    seurat_10X_integrated@assays$ADT <- cells@assays$ADT
  }
  
  # factorise clusterIDs
  print(unique(seurat_10X_integrated@meta.data[,temp_colname]))
  seurat_10X_integrated@meta.data[,temp_colname] <- 
    factor(seurat_10X_integrated@meta.data[,temp_colname],
           levels=str_sort(unique(seurat_10X_integrated@meta.data[,temp_colname]),numeric = T))
  print(unique(seurat_10X_integrated@meta.data[,temp_colname]))

  m <- paste0("seurat_10X_CITEimputed_",norm)
  assign(m, seurat_10X_integrated)
  rm(cells)
}

# VLNPLOTS/FEATUREPLOTS OF CITESEQ ----------------------------------------------------------------

for(norm in c("norm", "raw")){
    
    dir.create(paste0("01_CITE_markers_",norm))
    print(norm)
    seurat_10X_CITEimputed <- 
      get(paste0("seurat_10X_CITEimputed_", norm))
    
    DefaultAssay(seurat_10X_CITEimputed) <- "ADT"
    
    temp_CITE_markers <- 
      row.names(seurat_10X_CITEimputed$ADT@data)
    
    temp_total_plots <- ceiling(length(temp_CITE_markers)/12)
    
    for(i in c(1:temp_total_plots)) {
      #vlnplot
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 14,
            height = 10,
            res = 300,
            units = 'in'
          )
        }
      
      print(i)
      temp_index_1 <- i+((i-1)*11)
      temp_index_2 <- i+((i-1)*11)+11
      temp_index <- c(temp_index_1:temp_index_2)          
      
      # print(temp_index)
      temp_marker_subset <- temp_CITE_markers[temp_index]
      temp_marker_subset <- temp_marker_subset[! temp_marker_subset %in% NA]
      print(temp_marker_subset)
      
      temp_png_function(paste0("01_CITE_markers_",norm,
                               "/01_FILTERED_imputed_CITE_",i,
                               ".png"))
      
      temp_vlnplot <- VlnPlot(object = seurat_10X_CITEimputed, 
                              features = temp_marker_subset, 
                              pt.size = 0, 
                              group.by = temp_colname, 
                              log = F)
      print(temp_vlnplot)
      dev.off()
      
      # featureplot
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 16,
            height = 12,
            res = 300,
            units = 'in'
          )
        }
      
      temp_png_function(paste0("01_CITE_markers_",norm,
                               "/02_FILTERED_imputed_CITE_",i,
                               ".png"))
      
      temp_featureplot <- FeaturePlot(
        object = seurat_10X_CITEimputed,
        features = temp_marker_subset,
        pt.size = 0.01,
        reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
        min.cutoff = 0,
        order = T)
      print(temp_featureplot)
      dev.off()
      
      
      
    }
  }


# CLUSTER AVERAGED HEATMAP -----------------------------------------------------------------

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 10,
      height = 18,
      res = 300,
      units = 'in'
    )
  }

for(norm in c("norm", "raw")){
  
  print(norm)
  seurat_10X_CITEimputed <- 
    get(paste0("seurat_10X_CITEimputed_", norm))
  
  Idents(seurat_10X_CITEimputed) <- temp_colname
  
  temp_CITE_markers <-
    row.names(seurat_10X_CITEimputed$ADT@data)
  
  temp_cluster.averages <- 
    AverageExpression(seurat_10X_CITEimputed,
                      return.seurat = TRUE,
                      verbose = T,
                      assays="ADT",
                      features=temp_CITE_markers)
  
  temp_data_frame <- 
    as.matrix(temp_cluster.averages@assays$ADT@data)
  
  temp_data_frame_hclust <- as.data.frame(temp_data_frame)

  temp_rowsums <- rowSums(temp_data_frame_hclust)
  temp_rowsums <- temp_rowsums[temp_rowsums == 0]
  temp_rowsums <- names(temp_rowsums)
  
  # repeat
  temp_cluster.averages <- 
    AverageExpression(seurat_10X_CITEimputed,
                      return.seurat = TRUE,
                      verbose = T,
                      assays="ADT",
                      features=temp_CITE_markers[!temp_CITE_markers %in% temp_rowsums])
  
  temp_data_frame <- 
    as.matrix(temp_cluster.averages@assays$ADT@data)
  
  temp_data_frame_hclust <- as.data.frame(temp_data_frame)
  
  # plot
  library("viridis")  
  hmcol <- viridis(24)
  
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
  
  temp_png_function(paste0("01_CITE_markers_",norm,
                           "/03_pheatmap.png"))
  pheatmap(temp_data_frame_hclust,
                        color = rev(hmcol),
                        cluster_cols = T,
                        cluster_rows = T,
                        scale = "row",
                        clustering_distance_cols = "correlation",
                        clustering_distance_rows = "correlation",
                        fontsize_row = 4,
                        show_rownames = T,
                        show_colnames = T,
                        fontsize_col = 20,
                        annotation_legend = T,
                        cellheight= 6,
                        cellwidth = 20,
                        gaps_col = NULL,
                        annotation_names_col = T,
                        angle_col = 45,
                        legend = T,
                        border_color=FALSE)
  dev.off()
}

# SAVE OBJECTS -------------------------------------------------


for(norm in c("norm", "raw")){
  
  print(norm)
  seurat_10X_CITEimputed <- 
    get(paste0("seurat_10X_CITEimputed_", norm))
  
  saveRDS(seurat_10X_CITEimputed, 
            paste0("RDATA_01_IMPUTEDCITESEQ_", norm, ".Rdata"))
}

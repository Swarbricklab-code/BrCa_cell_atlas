# Celltype default reclustering script 
# 20191013
# R v3.4.1
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

if(temp_cell_type == "B_cells" | temp_cell_type == "Plasma_cells") {
  temp_PC_number <- 10
}
if(temp_cell_type == "T_cells") {
  temp_PC_number <- 20
}
if(temp_cell_type == "Myeloid_cells") {
  temp_PC_number <- 20
}

# temp_cell_type <- "B_cells"
# temp_wd <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/run02/analysis_01_B_cells/"

# SETUP -------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(data.table)
library("viridis")  
library(pheatmap)


# DIRECTORY ---------------------------------------------------------------

setwd(paste0(temp_wd))

dir.create("05_GO")
setwd("05_GO")


# LOAD DEGS  ----------------------------------------------------

temp_cluster.allmarkers <- 
  read.csv(paste0("../02_output_filtered/07_MAST_FINDALLMARKERS_BULK_res.0.8_RNA_assay.csv"))

temp_colname <- paste0("SUBSET", temp_PC_number, "_res.0.8")

table(temp_cluster.allmarkers$cluster)

temp_clusterMarkers <- NULL
for (cluster in as.character(unique(temp_cluster.allmarkers$cluster))){
  
  temp_signatureEntrezIDs <- bitr(temp_cluster.allmarkers[temp_cluster.allmarkers$cluster == cluster,]$gene,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = "org.Hs.eg.db"
  )
  
  temp_cm <-
    data.frame(EntrezID = temp_signatureEntrezIDs$ENTREZID,
               clusterID = cluster)
  if (!is.null(temp_clusterMarkers))
  {
    temp_clusterMarkers <- rbind(temp_clusterMarkers, 
                                 temp_cm)
  } else {
    temp_clusterMarkers <- temp_cm
  }}


# RUN GO ------------------------------------------------------------------

for(ont in c("groupGO", "enrichKEGG", "enrichDO", "enrichGO")){
  print(ont)
  if(! ont %in% c("enrichKEGG","enrichDO")){
    temp_go <- compareCluster(
      EntrezID ~ clusterID,
      data = temp_clusterMarkers,
      fun = ont,
      OrgDb = "org.Hs.eg.db",
      readable = T
    )
  }
  
  if(ont == "enrichKEGG"){
    temp_go <- compareCluster(
      EntrezID ~ clusterID,
      data = temp_clusterMarkers,
      fun = ont,
      organism = "hsa"
    )
  }
  if(ont == "enrichDO"){
    temp_go <- compareCluster(
      EntrezID ~ clusterID,
      data = temp_clusterMarkers,
      fun = ont
    )
  }
  
  temp_go_df <- as.data.frame(temp_go)
  write.csv(temp_go_df,
            paste0("01_compareCluster_",temp_cell_type,"_",ont,".csv"))
  n <- paste0("temp_go_df_",temp_cell_type, "_",ont)
  assign(n,temp_go_df)
}

# GO ENRICHMENT PLOT ------------------------------------------------------

# filter df and plot
  for(ont in c("enrichKEGG", "enrichDO", "enrichGO", "groupGO")){
    print(ont)
    
    temp_go_df <- get(paste0("temp_go_df_",temp_cell_type, "_",ont))
    
    for(type in c("Description", "ID")) {
      print(type)
      
      temp_go_df_sorted <- 
        temp_go_df[ order( temp_go_df[, "clusterID"],
                           temp_go_df[, type] ),]
      
      if(! ont == "groupGO"){
        temp_go_df_sorted$logp_val <- (-(log10(temp_go_df_sorted$p.adjust)))
      } else {
        temp_go_df_sorted$logp_val <- temp_go_df_sorted$Count
      }
      
      temp_go_df_sorted <- temp_go_df_sorted[ , c("logp_val",
                                                  "clusterID",
                                                  type)]
      
      rownames(temp_go_df_sorted) <- NULL
      temp_go_combined_sorted_temp <- temp_go_df_sorted
      
      colnames(temp_go_df_sorted) <- c("logp_val", "clusterID", "type")
      
      # plot
      if(! ont == "groupGO"){
        temp_combind_summary_df_top <- 
          (temp_go_df_sorted %>% group_by(clusterID) %>% top_n(20,logp_val))
        temp_cluster_rows <- T
      } else {
        temp_combind_summary_df_top <- temp_go_df_sorted
        temp_cluster_rows <- F
      }
      
      print("   Total number of pathways")
      print(length(as.vector(unique(temp_combind_summary_df_top$type))))
      temp_combind_summary_df_top <- temp_go_df_sorted[temp_go_df_sorted$type %in% as.vector(unique(temp_combind_summary_df_top$type)),]
      
      # assign value to every pathway for every cluster (0s)
      temp_combind_summary_df_dcast <-
        dcast(data = temp_combind_summary_df_top,
              formula = type~clusterID,
              fun.aggregate = sum,
              value.var = "logp_val")
      
      #row names to variables
      rownames(temp_combind_summary_df_dcast) <- 
        temp_combind_summary_df_dcast$type
      
      temp_combind_summary_df_dcast <- 
        temp_combind_summary_df_dcast[, ! colnames(temp_combind_summary_df_dcast) %in% "type"]
      
      temp_data_frame_hclust <- as.matrix(temp_combind_summary_df_dcast)
      
      
      # pheatmap
      hmcol <- viridis(24)
      
      temp_png_function <-
        function(x) {
          png(
            file = (x), 
            width = 12, 
            height = 18, 
            res = 300, 
            units = 'in'
          )
        }
      
      temp_png_function(paste0("02_compareCluster_",temp_cell_type,"_",ont,"_",type,".png"))
      pheatmap(temp_data_frame_hclust, 
               color = rev(hmcol),
               cluster_cols = T, 
               cluster_rows = temp_cluster_rows, 
               scale = "row", 
               fontsize_row = 8,
               clustering_distance_cols = "correlation",
               clustering_distance_rows = "correlation",
               show_colnames = T, 
               fontsize_col = 20,
               annotation_legend = F,
               cellheight= 7.5, 
               cellwidth = 20,
               gaps_col = NULL, 
               annotation_names_col = T, 
               angle_col = 45,
               treeheight_row = 50, 
               legend = T,
               border_color=FALSE)
      dev.off()
      
      
    }
  }

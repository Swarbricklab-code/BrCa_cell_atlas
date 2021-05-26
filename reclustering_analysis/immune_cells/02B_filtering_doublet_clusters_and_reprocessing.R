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

# PCs
if(temp_cell_type == "B_cells" | temp_cell_type == "Plasma_cells") {
  temp_unfiltered_PC_number <- 10
  temp_PC_number <- 10
}
if(temp_cell_type == "T_cells") {
  temp_unfiltered_PC_number <- 30
  temp_PC_number <- 30
}
if(temp_cell_type == "Myeloid_cells") {
  temp_unfiltered_PC_number <- 30
  temp_PC_number <- 20
}

# SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(future)

# DIRECTORY ---------------------------------------------------------------

setwd(paste0(temp_wd))

dir.create("02_output_filtered")
setwd("02_output_filtered")

# LOAD DATA ---------------------------------------------------------------

for(type in c("integrated")) {
  
  seurat_10X_integrated <- readRDS(paste0("../01_output_unfiltered/RDATA_01_PROCESSED_", type, "_UNFILTERED.Rdata"))
  
  n <- paste0("seurat_10X_integrated_", type)
  assign(n, seurat_10X_integrated)
}




# CLUSTERS TO FILTER OUT --------------------------------------------------

for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_", type))
  
  Idents(seurat_10X_integrated) <- paste0("SUBSET", temp_unfiltered_PC_number, "_res.0.8")
  
  temp_ids <- unique(seurat_10X_integrated@meta.data[,paste0("SUBSET", temp_unfiltered_PC_number, "_res.0.8")])
  
  if(temp_cell_type == "B_cells") {
    # 6 t -cells
    # 8 poor quality
    temp_ids <- temp_ids[!temp_ids %in% c(6,8)]
  }
  if(temp_cell_type == "Plasma_cells") {
    # 0,2,7,10 EPCAM
    # 9,11 t -cells
    # 8 - collagens
    # 0, 10, 11, 2, 7, 8, 9 low UMIs
    temp_ids <- temp_ids[!temp_ids %in% c(0,2,7,8,9,10,11)]
  }
  if(temp_cell_type == "T_cells") {
    # 5,10 - poor quality & CAF and EPCAM markers & high mito
    # 14 myeloid proliferating cells
    # 4 JCHAIN and IGKs DRIVEN - plasma cells ?
    # 7,9 no DEGs high mito
    temp_ids <- temp_ids[!temp_ids %in% c(4,9,7,5,10,14)]
  }
  if(temp_cell_type == "Myeloid_cells") {
    # OLD v1 filtering notes
    # 22-26 small noisey clusters
    # 11 CAFs and T-cells
    # 9, 20  and 15 T-cell
    # 14 KRTs
    # 19,21 low quality
    # 7 low quality and T-cells

    # new filtering notes (with prol myeloid cells)
    # 11 T-cells
    # 9 KRTs
    # 3 all IG DEGs
    # 10, 14, 15,2,3,7 - poor quality
    # 18-27 small noisey clusters
    # 7, 15, 23:27 - no DEGs
    temp_ids <- temp_ids[!temp_ids %in% c(3,7,9,10,11,14,15,20:27)]
  }
  
  print("filtering out clusters")
  print(temp_ids)
}

# FILTER OUT LOW UMI CLUSTERS AND OTHER LINEAGE CLUSTERS ---------------------------------------

for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_", type))
  
  Idents(seurat_10X_integrated) <- paste0("SUBSET", temp_unfiltered_PC_number, "_res.0.8")
  
  print("raw")
  print(dim(seurat_10X_integrated))
  
  seurat_10X_integrated <- 
    SubsetData(seurat_10X_integrated,
               ident.use=temp_ids)  
  
  print("filtered")
  print(dim(seurat_10X_integrated))

  # remove sampleIDs with < 30 cells
  print(table(seurat_10X_integrated@meta.data$orig.ident))
  temp_df <- as.data.frame(table(seurat_10X_integrated@meta.data$orig.ident))
  temp_df_filtered <- temp_df[temp_df$Freq > 30,]
  temp_df_filtered <- as.vector(temp_df_filtered$Var1)
  temp_removed <- as.vector(temp_df$Var1[!temp_df$Var1 %in% temp_df_filtered])
  
  if(length(temp_removed) > 0){
    print("trimming samples")
    Idents(seurat_10X_integrated) <- "orig.ident"
    seurat_10X_integrated <- SubsetData(seurat_10X_integrated,
                                        ident.use=temp_df_filtered)
    print(temp_removed)
  }
  
  print(dim(seurat_10X_integrated))
  # print(table(seurat_10X_integrated@meta.data$orig.ident))
  
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
}


# REPROCESS FILTERED DATA -------------------------------------------------

# for(type in c("garnett", "integrated")) {
for(type in c("integrated")) {
  print(type)
  temp_seurat_10X_subset <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  
  temp_subset_cells_list <- SplitObject(temp_seurat_10X_subset, 
                                        split.by = "orig.ident")
  
  # reidentify variable features
  for (i in 1:length(temp_subset_cells_list)) {
    print(unique(temp_subset_cells_list[[i]]@meta.data$orig.ident))
    DefaultAssay(temp_subset_cells_list[[i]] ) <- "RNA"
    
    temp_subset_cells_list[[i]]@assays$integrated <- NULL
    
    temp_subset_cells_list[[i]] <- 
      NormalizeData(temp_subset_cells_list[[i]], 
                    verbose = T)
    
    temp_subset_cells_list[[i]] <- 
      FindVariableFeatures(temp_subset_cells_list[[i]], 
                           nfeatures = 5000,
                           verbose = F,
                           selection.method = "dispersion")
    print(length(temp_subset_cells_list[[i]]@assays$RNA@var.features))
  }
  
  
  if(temp_run_parallel){
    plan("multiprocess", 
         workers = temp_run_parallel_cores)
    options(future.globals.maxSize = temp_mem_per_core*2 * 1024^3)
  }
  
  # For integrating datasets that have low cell numbers (<200), you will 
  # encounter errors from the k.filter param in the FindIntegrationAnchors for 
  # nearest neighbours this breaks if you are trying to integrate a sample with 
  # <200 cells as there is not that many neighbours to search. Use the following
  # param to set this min to the min number of cells
  k.filter <- 
    min(200, min(sapply(temp_subset_cells_list, ncol)))
  print(k.filter)
  
  # if(k.filter < 30) {
  #   print(Sys.time())
  #   seurat_10X_anchors <- FindIntegrationAnchors(object.list = temp_subset_cells_list, 
  #                                                k.filter = k.filter, 
  #                                                dims = 1:(k.filter-1))
  #   print(Sys.time())
  #   
  # }
  
  if(k.filter >= 30) {
    
    print(Sys.time())
    seurat_10X_anchors <- FindIntegrationAnchors(object.list = temp_subset_cells_list, 
                                                 k.filter = k.filter,
                                                 anchor.features=5000, 
                                                 verbose = FALSE)
    print(Sys.time())
  }
  
  if(temp_run_parallel){
    plan("sequential")
  }
  
  temp_integration_start <- Sys.time()
  seurat_10X_integrated <- IntegrateData(anchorset = seurat_10X_anchors, 
                                         verbose = FALSE)
  temp_integration_finish <- Sys.time()
  
  # set assay back to integrated
  DefaultAssay(object = seurat_10X_integrated) <- 
    "integrated"
  
  # re-scale data
  seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated, 
                                     verbose = FALSE,
                                     vars.to.regress=c("nCount_RNA", "nFeature_RNA", "percent.mito"))
  
  # PCA 
  seurat_10X_integrated <- RunPCA(object = seurat_10X_integrated, 
                                  npcs = temp_PC_number, 
                                  verbose = T)
  
  # dimensional reduction
  # TSNE
  seurat_10X_integrated <-
    RunTSNE(object = seurat_10X_integrated,
            dims = 1:temp_PC_number,
            reduction.key = paste0("TSNESIG",temp_PC_number,"_"),
            reduction.name = paste0("TSNESIG",temp_PC_number)
    )
  # UMAP
  seurat_10X_integrated <- 
    RunUMAP(seurat_10X_integrated, 
            dims =  1:temp_PC_number,
            reduction.key = paste0("UMAPSIG",temp_PC_number,"_"),
            reduction.name = paste0("UMAPSIG",temp_PC_number),
            verbose = F
    )
  
  # reclust 
  seurat_10X_integrated <- 
    FindNeighbors(object = seurat_10X_integrated,
                  dims = 1:temp_PC_number,
                  graph.name = paste0("SUBSET", temp_PC_number)
    )
  
  # default clustering
  seurat_10X_integrated <-
    FindClusters(
      object = seurat_10X_integrated,
      graph.name = paste0("SUBSET", temp_PC_number), 
      resolution = 0.8
    )
  
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
}


# FILTERED DIMPLOTS AND FEATUREPLOTS -----------------------------------------------

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 14,
      height = 12,
      res = 300,
      units = 'in'
    )
  }

# for(type in c("garnett", "integrated")) {
for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  
  # DIMPLOT
  temp_png_function(paste0("04_FILTERED_UMAP_C_", type, ".png"))
  temp_dimplot <- DimPlot(
    object = seurat_10X_integrated,
    label.size = 4,
    pt.size = 0.5,
    reduction = paste0("TSNESIG",temp_PC_number),
    label = T,
    group.by = paste0("SUBSET", temp_PC_number, "_res.0.8")
  )
  print(temp_dimplot)
  dev.off()
  
  # FEATUREPLOT
  
  if(temp_cell_type == "B_cells"| temp_cell_type == "Plasma_cells") {
    temp_markers <-  c("CD19", "MS4A1", "BLK", "FCRL2", "TCL1A", "JCHAIN", "IGKC", "CD79A", "TNFRSF17")
  }
  if(temp_cell_type == "T_cells") {
    temp_markers <-  c("CD3D", "ITGAE", "CD4", "CD8A", "CD8B", "MKI67", "FOXP3", "NCR1", "XCL1", "KLRK1", "LAG3", "PDCD1", "CD274", "TOX", "CD86", "IL2RA")
  }
  if(temp_cell_type == "Myeloid_cells") {
    temp_markers <-  c("CD68","ITGAM","CD74", "LILRB4", "FCER1A", "CD163", "IFITM3", "CLEC9A", "S100A12", "RGCC", "PDCD1", "CD274", "CD80", "CD86",  "CD84")
    
    temp_markers_2 <-  c("C5AR1","FCAR","HLA-DQA1","CD81","CD14","FLT3","CADM1","IL3RA","CD1C","IRF8","IRF4","THBD","CD5","FCGR3A")
    
  }  
  
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  temp_png_function(paste0("05_FILTERED_featureplots_01_",type, ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0.25,
    reduction = paste0("TSNESIG",temp_PC_number),
    min.cutoff = 0,
    order = T)
  print(temp_featureplot)
  dev.off()
  
  if(temp_cell_type == "Myeloid_cells") {
    temp_png_function(paste0("05_FILTERED_featureplots_02_",type, ".png"))
    temp_featureplot <- FeaturePlot(
      object = seurat_10X_integrated,
      features = temp_markers_2,
      pt.size = 0.25,
      reduction = paste0("TSNESIG",temp_PC_number),
      min.cutoff = 0,
      order = T)
    print(temp_featureplot)
    dev.off()
    }
  
  temp_markers <- 
    c("PTPRC", "CD3D", "CD3E", "CD3G", "KRT8", "KRT18", "EPCAM", "KRT5", "KRT14", "KRT17","JCHAIN")
  temp_png_function(paste0("05_FILTERED_featureplots_02_nonlineage_markers_",type, ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0.25,
    reduction = paste0("TSNESIG",temp_PC_number),
    min.cutoff = 0,
    order = T)
  print(temp_featureplot)
  dev.off()
  
  for(metadata in c("orig.ident", "clinical_subtype")) {
    
    temp_png_function(paste0("06_FILTERED_",type,"_split_by_", metadata, ".png"))
    temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      label.size = 4,
      pt.size = 0.5,
      reduction = paste0("TSNESIG",temp_PC_number),
      split.by = paste0(metadata),
      label = F,
      group.by = paste0("SUBSET", temp_PC_number, "_res.0.8") 
    )
    print(temp_dimplot)
    dev.off()
  }
  
  
  
}


# RESCALE RNA -------------------------------------------------------------

# # rescale RNA assay
# for(type in c("integrated")) {
#   
#   print(type)
#   seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
#   
#   DefaultAssay(object = seurat_10X_integrated) <-
#     "RNA"
#   
#   # re-scale data
#   seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated,
#                                      verbose = FALSE,
#                                      vars.to.regress=c("nCount_RNA", "nFeature_RNA", "percent.mito"),
#                                      features = rownames(seurat_10X_integrated))
#   
#   n <- paste0("seurat_10X_integrated_",type,"_filtered")
#   assign(n, seurat_10X_integrated)
#   
# }

# SAVE RESCALED OBJECTS -------------------------------------------------

for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  
  saveRDS(seurat_10X_integrated, 
          paste0("RDATA_02_PROCESSED_FILTERED_RESCALED", type, ".Rdata"))
  
}

# DEG AT RES 0.8 ---------------------------------------------------------------------

# Integrated
for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  DefaultAssay(seurat_10X_integrated) <- "integrated"
  Idents(seurat_10X_integrated) <- paste0("SUBSET", temp_PC_number, "_res.0.8")
  
  temp_cluster.allmarkers <- FindAllMarkers(
    only.pos = T,
    object = seurat_10X_integrated,
    min.pct = 0.1, 
    logfc.threshold = 0.5,
    min.diff.pct = 0.1, 
    test.use = 'MAST', 
    print.bar = T
  ) 
  
  temp_cluster.allmarkers <- arrange(temp_cluster.allmarkers,
                                     (cluster),
                                     desc(avg_logFC))
  
  write.csv(temp_cluster.allmarkers,
            paste0("07_MAST_FINDALLMARKERS_BULK_res.0.8.csv"))
  
  temp_genes_for_heatmap <- 
    (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(20,
                                                             avg_logFC))
  
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 12,
        height = 16,
        res = 300,
        units = 'in'
      )
    }
  
  temp_png_function(paste0("07_MAST_FINDALLMARKERS_BULK_res.0.8_integrated.png"))
  print(DoHeatmap(
    object = seurat_10X_integrated,
    features = temp_genes_for_heatmap$gene,
    group.by = paste0("SUBSET", temp_PC_number, "_res.0.8"),
    raster=F)
  )
  dev.off()
  
}


# RNA
for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  Idents(seurat_10X_integrated) <- paste0("SUBSET", temp_PC_number, "_res.0.8")
  
  temp_cluster.allmarkers <- FindAllMarkers(
    only.pos = T,
    object = seurat_10X_integrated,
    min.pct = 0.1, 
    logfc.threshold = 0.25,
    min.diff.pct = 0.05, 
    test.use = 'MAST', 
    print.bar = T
  ) 
  
  temp_cluster.allmarkers <- arrange(temp_cluster.allmarkers,
                                     (cluster),
                                     desc(avg_logFC))
  
  write.csv(temp_cluster.allmarkers,
            paste0("07_MAST_FINDALLMARKERS_BULK_res.0.8_RNA_assay.csv"))
  
  temp_genes_for_heatmap <- 
    (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(20,
                                                             avg_logFC))
  
  seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated,
                                     verbose = FALSE,
                                     vars.to.regress=c("nCount_RNA", "nFeature_RNA", "percent.mito"),
                                     features = as.vector(temp_genes_for_heatmap$gene))
  
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 12,
        height = 16,
        res = 300,
        units = 'in'
      )
    }
  
  temp_png_function(paste0("07_MAST_FINDALLMARKERS_BULK_res.0.8_RNA.png"))
  print(DoHeatmap(
    object = seurat_10X_integrated,
    features = temp_genes_for_heatmap$gene,
    group.by = paste0("SUBSET", temp_PC_number, "_res.0.8"),
    raster=F)
  )
  dev.off()
  
}

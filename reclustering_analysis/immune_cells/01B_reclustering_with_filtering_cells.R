# Celltype default reclustering script 
# 20191013
# explore the differences between stromal subclasses from different subtypes of breast cancer and normal tissue
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
print(temp_cell_type)

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
}
if(temp_cell_type == "T_cells") {
  temp_PC_number <- 30
}
if(temp_cell_type == "Myeloid_cells") {
  temp_PC_number <- 30
}

# SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(future)

# LOAD LARGE ALIGNED DATASET ---------------------------------------------------------------

seurat_10X <-
  readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/02_integration/run02_new26samples/output/CCA_CCA23Jun2019/Output/Rdata/03_seurat_CCA_aligned_processed.Rdata")

# old object
# "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Jun2019/02_integration/output_run02_success/CCA_CCA23Jun2019/Output/Rdata/03_seurat_CCA_aligned_processed.Rdata"

# NEW CLINICAL METADATA ---------------------------------------------------

seurat_10X@meta.data$clinical_subtype <- seurat_10X@meta.data$subtype

# Change to metaplastic and normals and collapse HER2+
# seurat_10X@meta.data$clinical_subtype[seurat_10X@meta.data$orig.ident %in% "CID4523"] <- "TNBC"
# seurat_10X@meta.data$clinical_subtype[seurat_10X@meta.data$orig.ident %in% "CID4513"] <- "TNBC"

# seurat_10X@meta.data$clinical_subtype[seurat_10X@meta.data$orig.ident %in% "CID4530N"] <- "Normal"
# seurat_10X@meta.data$clinical_subtype[seurat_10X@meta.data$orig.ident %in% "CID4520N"] <- "Normal"
# seurat_10X@meta.data$clinical_subtype[seurat_10X@meta.data$orig.ident %in% "CID4378N"] <- "Normal"
# seurat_10X@meta.data$clinical_subtype[seurat_10X@meta.data$orig.ident %in% "CID4523N"] <- "Normal"


# DIRECTORY ---------------------------------------------------------------

setwd(paste0(temp_wd))

dir.create("01_output_unfiltered")
setwd("01_output_unfiltered")

# LOAD INTEGRATED CLUSTER DATA AND SUBSET CAFs AND ENDOTHELIAL CELLS --------------------------------------

for(type in c("integrated")) {
  
    Idents(object = seurat_10X) <-
      "int_cluster_major_PC_B_res.0.8"
  
    if(temp_cell_type == "B_cells") {
      temp_cell_type_ID <- c("B cells")
      cells_to_trim <- 30
    }
    if(temp_cell_type == "Plasma_cells") {
      temp_cell_type_ID <- c("Plasma cells")
      cells_to_trim <- 30
    }
    if(temp_cell_type == "T_cells") {
      temp_cell_type_ID <- "T cells"
      cells_to_trim <- 30

    }
    if(temp_cell_type == "Myeloid_cells") {
      temp_cell_type_ID <- "Myeloid cells"
      cells_to_trim <- 30
    }
    
    temp_seurat_10X_subset <- SubsetData(seurat_10X,
                                         ident.use=temp_cell_type_ID, 
                                     assay = "RNA")
    
    # Proliferating myeloid cells clustered with prolfierating T-cells in the above resolution seurat garnett calling
    # they do not in cluster 3 of UMAP_PC_A_res.0.4
    # include the differential cells here
    if(temp_cell_type == "Myeloid_cells") {
      Idents(object = seurat_10X) <-
        "int_PCA_A_res.0.4"
      
      temp_seurat_10X_subset_temp <- SubsetData(seurat_10X,
                                           ident.use=3, 
                                           assay = "RNA")
      
      temp_barcodes <- unique(c(colnames(temp_seurat_10X_subset),
                              colnames(temp_seurat_10X_subset_temp)))
      
      temp_seurat_10X_subset <- SubsetData(seurat_10X,
                                           cells=temp_barcodes, 
                                           assay = "RNA")
      rm(temp_seurat_10X_subset_temp)
    }
    

  print(dim(temp_seurat_10X_subset))
  print(table(temp_seurat_10X_subset@meta.data$garnett_cell))
  print(table(temp_seurat_10X_subset@meta.data$orig.ident))
  
  # remove sampleIDs with < 30 cells
  temp_df <- as.data.frame(table(temp_seurat_10X_subset@meta.data$orig.ident))
  temp_df_filtered <- temp_df[temp_df$Freq > cells_to_trim,]
  temp_df_filtered <- as.vector(temp_df_filtered$Var1)
  temp_removed <- as.vector(temp_df$Var1[!temp_df$Var1 %in% temp_df_filtered])
  
  if(length(temp_removed) > 0){
    print("trimming samples")
  Idents(temp_seurat_10X_subset) <- "orig.ident"
  temp_seurat_10X_subset <- SubsetData(temp_seurat_10X_subset,
                                       ident.use=temp_df_filtered)
  print(temp_removed)
  }
  
  # downsample cells to 42k cells number of cells works
  # n.cells <- 
  #   round(ncol(temp_seurat_10X_subset)*0.75)
  # temp_seurat_10X_subset <- SubsetData(temp_seurat_10X_subset, 
  #                                      cells = sample(x = colnames(temp_seurat_10X_subset),      
  #                                                         size = n.cells))
  # ncol(temp_seurat_10X_subset)
  

  n <- 
    paste0("temp_seurat_10X_subset_", type)
  assign(n, 
         temp_seurat_10X_subset)
  # rm(seurat_10X)
}
  
# RECLUSTERING ------------------------------------------------------------

for(type in c("integrated")) {
    print(type)
    temp_seurat_10X_subset <- get(paste0("temp_seurat_10X_subset_", type))
  
    temp_subset_cells_list <- SplitObject(temp_seurat_10X_subset,
                                          split.by = "orig.ident")
    
    # split object
    # Idents(temp_seurat_10X_subset) <- "orig.ident"
    # for(sampleID in unique(temp_seurat_10X_subset@meta.data$orig.ident)){
    #   temp_seurat_10X_subset_subset <- subset(temp_seurat_10X_subset,
    #                                           ident=sampleID)
    #   n <- paste0("temp_seurat_10X_subset_subset_",sampleID)
    #   assign(n, temp_seurat_10X_subset_subset)
    #   rm(temp_seurat_10X_subset_subset)
    # }
    # temp_subset_cells_list <-
    #   mget(ls(pattern = "temp_seurat_10X_subset_subset_CID"))

    # temp_subset_cells_list <- temp_subset_cells_list[temp_order]
    
    # reidentify variable features
    for (i in 1:length(temp_subset_cells_list)) {
      
      print(unique(temp_subset_cells_list[[i]]@meta.data$orig.ident))
      DefaultAssay(temp_subset_cells_list[[i]]) <- "RNA"
      
      # remove integrated assay to avoid this subsequent error during FindIntegrationANchors
      # Error in RowMergeMatrices(mat1 = mat1, mat2 = mat2, mat1_rownames = mat1.names,  : Need S4 class dgRMatrix for a sparse matrix
      # issue here but no solution https://github.com/satijalab/seurat/issues/1355
      
      temp_subset_cells_list[[i]]@assays$integrated <- NULL
      
      temp_subset_cells_list[[i]] <- 
        NormalizeData(temp_subset_cells_list[[i]], 
                      verbose = F)
      
      temp_subset_cells_list[[i]] <- 
        FindVariableFeatures(temp_subset_cells_list[[i]], 
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
    # 
    #   seurat_10X_anchors <- FindIntegrationAnchors(object.list = temp_subset_cells_list)
    # 
    #   print(Sys.time())
    # 
    # }
    
    if(k.filter >= 30) {

    print(Sys.time())
    seurat_10X_anchors <- FindIntegrationAnchors(object.list = temp_subset_cells_list, 
                                                 k.filter = k.filter,
                                                 verbose = F)
    print(Sys.time())
    }
    
    
    
    if(temp_run_parallel){
      plan("sequential")
    }
    
    temp_integration_start <- Sys.time()
    seurat_10X_integrated <- IntegrateData(anchorset = seurat_10X_anchors)
    temp_integration_finish <- Sys.time()
    
    # set assay back to integrated
    DefaultAssay(object = seurat_10X_integrated) <- 
      "integrated"
    
    # re-scale data
    if(temp_run_parallel){
      plan("multiprocess", 
           workers = temp_run_parallel_cores)
      options(future.globals.maxSize = temp_mem_per_core * 1024^3)
    }
    seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated, 
                                       verbose = FALSE,
                                       vars.to.regress=c("nCount_RNA", "nFeature_RNA", "percent.mito"))
    if(temp_run_parallel){
      plan("sequential")
    }
    
    # PCA 
    seurat_10X_integrated <- RunPCA(object = seurat_10X_integrated, 
                                    npcs = 100, 
                                    verbose = T)
    
    # dimensional reduction
    # UMAP
    seurat_10X_integrated <- 
      RunUMAP(seurat_10X_integrated, 
              dims =  1:temp_PC_number,
              reduction.key = paste0("UMAPSIG",temp_PC_number,"_"),
              reduction.name = paste0("UMAPSIG",temp_PC_number),
              verbose = F
      )
    
    # TSNE
    seurat_10X_integrated <-
      RunTSNE(object = seurat_10X_integrated,
              dims = 1:temp_PC_number,
              reduction.key = paste0("TSNESIG",temp_PC_number,"_"),
              reduction.name = paste0("TSNESIG",temp_PC_number)
      )
    
    # reclust 
    seurat_10X_integrated <- 
      FindNeighbors(object = seurat_10X_integrated,
                    dims = 1:temp_PC_number,
                    graph.name = paste0("SUBSET", temp_PC_number)
      )
    
    # default clustering at 0.8
    seurat_10X_integrated <-
      FindClusters(
        object = seurat_10X_integrated,
        graph.name = paste0("SUBSET", temp_PC_number),
        resolution = 0.8
      )
    
    n <- 
      paste0("seurat_10X_integrated_", type)
    assign(n, 
           seurat_10X_integrated)
}

# SAVE UNFILTERED OBJECTS -------------------------------------------------

for(type in c("integrated")) {
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_", type))
  
  saveRDS(seurat_10X_integrated, paste0("RDATA_01_PROCESSED_", type, "_UNFILTERED.Rdata"))
  
}

# DIMPLOTS AND FEATUREPLOTS -----------------------------------------------

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

for(type in c("integrated")) {
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_", type))

  # DIMPLOT
  temp_png_function(paste0("01_TSNE_C_", type, ".png"))
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
  if(temp_cell_type == "B_cells" | temp_cell_type == "Plasma_cells") {
    temp_markers <-  c("CD19", "MS4A1", "BLK", "FCRL2", "TCL1A", "JCHAIN", "IGKC", "CD79A", "TNFRSF17")
  }
  if(temp_cell_type == "T_cells") {
    temp_markers <-  c("CD3D", "ITGAE", "CD4", "CD8A", "CD8B", "MKI67", "FOXP3", "NCR1", "XCL1", "KLRK1", "LAG3", "PDCD1", "CD274", "TOX", "CD86", "IL2RA")
  }
  if(temp_cell_type == "Myeloid_cells") {
    temp_markers <-  c("CD68", "ITGAM", "CD74", "CD163", "CD84", "CSF3R", "LAG3", "PDCD1", "CD274", "CD80", "CD86")
  }
  
  temp_png_function(paste0("02_featureplots_01_markers_",type, ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0.25,
    reduction = paste0("TSNESIG",temp_PC_number),
    min.cutoff = 0,
    order = T)
  print(temp_featureplot)
  dev.off()
  
  # RNA
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  temp_markers <- 
    c("PTPRC", "CD3D", "CD3E", "CD3G", "KRT8", "KRT18", "EPCAM", "KRT5", "KRT14", "KRT17", "PDGFRB", "PDGFRA", "COL1A1", "MS4A1", "CD68", "CD14")
  temp_png_function(paste0("02_featureplots_02_nonlineage_markers_RNA",type, ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0.25,
    reduction = paste0("TSNESIG",temp_PC_number),
    min.cutoff = 0,
    order = T)
  print(temp_featureplot)
  dev.off()
  
  
  
}




# METRICS
temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 12,
      height = 12,
      res = 300,
      units = 'in'
    )
  }

temp_png_function_2 <-
  function(x) {
    png(
      file = (x),
      width = 15,
      height = 5,
      res = 300,
      units = 'in'
    )
  }

for(type in c("integrated")) {
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_", type))
  

  temp_png_function(paste0("03_featureplots_metrics_",type, ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mito"),
    pt.size = 1,
    reduction = paste0("TSNESIG",temp_PC_number),
    min.cutoff = 0, 
    order = T)
  print(temp_featureplot)
  dev.off()
  
  
  temp_png_function_2(paste0("03_vlnplots_metrics_",type, ".png"))
  temp_vlnplot <- VlnPlot(seurat_10X_integrated,
                          features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
                          pt.size = 0,
                          ncol = 3, 
                          group.by = paste0("SUBSET", temp_PC_number, "_res.0.8"), 
                          log = T)
  print(temp_vlnplot)
  dev.off()

}



# DEG AT RES 0.8 ---------------------------------------------------------------------

# Integrated
for(type in c("integrated")) {
  
  print(type)
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_", type))
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
            paste0("04_MAST_FINDALLMARKERS_BULK_res.0.8.csv"))
  
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
  
  temp_png_function(paste0("04_MAST_FINDALLMARKERS_BULK_res.0.8_integrated.png"))
  print(DoHeatmap(
    object = seurat_10X_integrated,
    features = temp_genes_for_heatmap$gene,
    group.by = paste0("SUBSET", temp_PC_number, "_res.0.8"),
    raster=F)
  )
  dev.off()
  
}






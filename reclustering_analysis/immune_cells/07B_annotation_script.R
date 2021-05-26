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
# temp_wd <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/run04/analysis_01_Myeloid_cells/"



# SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(future)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(grid)
library(ggplot2)

# DIRECTORY ---------------------------------------------------------------

setwd(paste0(temp_wd))

dir.create("07_annotated_objects")
setwd("07_annotated_objects")

# LOAD DATA ---------------------------------------------------------------

temp_colname <- paste0("SUBSET", temp_PC_number, "_res.",temp_res)

for(type in c("integrated")) {
  
  seurat_10X_integrated <- readRDS(paste0("../06_imputed_CITESeq/RDATA_01_IMPUTEDCITESEQ_raw.Rdata"))
  
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
}


# ANNOTATE CLUSTERS -------------------------------------------------------


for(type in c("integrated")) {
  
  seurat_10X_integrated <- readRDS(paste0("../06_imputed_CITESeq/RDATA_01_IMPUTEDCITESEQ_raw.Rdata"))
  Idents(seurat_10X_integrated) <- temp_colname
  
  if(temp_cell_type == "B_cells"){
    seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
                                          `0` = "Memory B-cells",
                                          `1` = "Memory B-cells",
                                          `2` = "Memory B-cells",
                                          `3` = "Naive B-cells",
                                          `4` = "Memory B-cells",
                                          `5` = "Naive B-cells",
                                          `6` = "Memory B-cells",
                                          `7` = "Memory B-cells",
                                          `8` = "Memory B-cells",
                                          `9` = "Memory B-cells",
                                          `10` = "Memory B-cells")
    temp_major_label <- "B-cells"
    seurat_10X_integrated@meta.data$celltype_major <- temp_major_label
    seurat_10X_integrated@meta.data$celltype_minor <- seurat_10X_integrated@active.ident
    seurat_10X_integrated@meta.data$celltype_subset <- seurat_10X_integrated@active.ident
  }
  if(temp_cell_type == "Plasma_cells"){
    # seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
    #                                       `0` = "Plasmablasts",
    #                                       `1` = "Plasmablasts",
    #                                       `2` = "Plasmablasts",
    #                                       `3` = "Plasmablasts",
    #                                       `4` = "Plasmablasts",
    #                                       `5` = "Plasmablasts",
    #                                       `6` = "Plasmablasts",
    #                                       `7` = "Plasmablasts",
    #                                       `8` = "Plasmablasts")
    temp_major_label <- "Plasmablasts"
    seurat_10X_integrated@meta.data$celltype_major <- temp_major_label
    seurat_10X_integrated@meta.data$celltype_minor <- temp_major_label
    seurat_10X_integrated@meta.data$celltype_subset <- temp_major_label
    
  }
  if(temp_cell_type == "T_cells"){
    seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
                                          `0` = "CD8+ T-cells 1",
                                          `1` = "CD4+ T-cells 1",
                                          `2` = "CD4+ T-cells 2",
                                          `3` = "T-Regs",
                                          `4` = "CD4+ T-cells 2", # collapsed with c2 per simons suggestion
                                          `5` = "CD8+ T-cells 5",
                                          `6` = "Tfh cells",
                                          `7` = "Unassigned 1",
                                          `8` = "Exhausted CD8+ T-cells",
                                          `9` = "NK cells",
                                          `10` = "Cycling T-cells",
                                          `11` = "Unassigned 2", # to remove per simons suggestion
                                          `12` = "NKT",
                                          `13` = "CD8+ T-cells 4",
                                          `14` = "CD8+ T-cells 2",
                                          `15` = "Unassigned 6", # removed c15, which was unique to mainly 1 patient and had weird DEG per simons suggestion
                                          `16` = "Unassigned 3",
                                          `17` = "Unassigned 4",
                                          `18` = "Unassigned 5")
    seurat_10X_integrated@meta.data$celltype_subset <- seurat_10X_integrated@active.ident
    
    Idents(seurat_10X_integrated) <- temp_colname
    seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
                                          `0` = "CD8+ T-cells",
                                          `1` = "CD4+ T-cells",
                                          `2` = "CD4+ T-cells",
                                          `3` = "CD4+ T-cells",
                                          `4` = "CD4+ T-cells",
                                          `5` = "CD8+ T-cells",
                                          `6` = "CD4+ T-cells",
                                          `7` = "Unassigned",
                                          `8` = "CD8+ T-cells",
                                          `9` = "NK cells",
                                          `10` = "Cycling",
                                          `11` = "Unassigned",
                                          `12` = "NKT",
                                          `13` = "CD8+ T-cells",
                                          `14` = "CD8+ T-cells",
                                          `15` = "CD8+ T-cells",
                                          `16` = "Unassigned",
                                          `17` = "Unassigned",
                                          `18` = "Unassigned")
    seurat_10X_integrated@meta.data$celltype_minor <- seurat_10X_integrated@active.ident
    
    temp_major_label <- "T-cells"
    seurat_10X_integrated@meta.data$celltype_major <- temp_major_label
    
  }
  if(temp_cell_type == "Myeloid_cells"){
    seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
                                          `0` = "Macrophage 1",
                                          `1` = "Monocyte 1",
                                          `2` = "LAM 2",
                                          `3` = "LAM 1",
                                          `4` = "LAM 1",
                                          `5` = "Monocyte 2",
                                          `6` = "Macrophage 2",
                                          `7` = "Cycling Myeloid",
                                          `8` = "cDC2",
                                          `9` = "pDC",
                                          `10` = "cDC1",
                                          `11` = "LAMP3+ DC",
                                          `12` = "Unassigned",
                                          `13` = "Macrophage 3")
    seurat_10X_integrated@meta.data$celltype_subset <- seurat_10X_integrated@active.ident
    
    Idents(seurat_10X_integrated) <- temp_colname
    seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
                                          `0` = "Macrophage",
                                          `1` = "Monocyte",
                                          `2` = "Macrophage",
                                          `3` = "Macrophage",
                                          `4` = "Macrophage",
                                          `5` = "Monocyte",
                                          `6` = "Macrophage",
                                          `7` = "Cycling",
                                          `8` = "DCs",
                                          `9` = "DCs",
                                          `10` = "DCs",
                                          `11` = "DCs",
                                          `12` = "Unassigned",
                                          `13` = "Macrophage")
    seurat_10X_integrated@meta.data$celltype_minor <- seurat_10X_integrated@active.ident
    
    
    temp_major_label <- "Myeloid"
    seurat_10X_integrated@meta.data$celltype_major <- temp_major_label
    
  }
  
  
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
  
}


# RECLUSTER SPECIFIC CELLS TYPES ------------------------------------------------

if(temp_cell_type == "Myeloid_cells"){
  
    for(type in c("integrated")) {
  
    temp_PC_number <- 20
    
    seurat_10X_integrated <- get( paste0("seurat_10X_integrated_",type,"_filtered"))
    Idents(seurat_10X_integrated) <- "celltype_subset"
    seurat_10X_integrated_subset <- subset(seurat_10X_integrated,
                                           idents="Monocyte 2")
    
    DefaultAssay(seurat_10X_integrated_subset) <- "integrated"
    
    seurat_10X_integrated_subset <- ScaleData(object = seurat_10X_integrated_subset, 
                                       verbose = T)
    
    # PCA 
    seurat_10X_integrated_subset <- RunPCA(object = seurat_10X_integrated_subset, 
                                    npcs = temp_PC_number, 
                                    verbose = T)
    
    seurat_10X_integrated_subset <- 
      RunUMAP(seurat_10X_integrated_subset, 
              dims =  1:temp_PC_number,
              reduction.key = paste0("UMAPSIG",temp_PC_number,"_"),
              reduction.name = paste0("UMAPSIG",temp_PC_number),
              verbose = F
      )
    
    # reclust 
    seurat_10X_integrated_subset <- 
      FindNeighbors(object = seurat_10X_integrated_subset,
                    dims = 1:temp_PC_number,
                    graph.name = paste0("SUBSET", temp_PC_number)
      )
    
    # default clustering
    seurat_10X_integrated_subset <-
      FindClusters(
        object = seurat_10X_integrated_subset,
        graph.name = paste0("SUBSET", temp_PC_number), 
        resolution = 0.4
      )
    
    
    
    # plot
    dir.create("monocyte_3_cluster")

    # DIMPLOT
    for(reduction in c("UMAP")){
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 6,
            height = 6,
            res = 300,
            units = 'in'
          )
        }
      
      temp_png_function(paste0("monocyte_3_cluster/01_",reduction,"_dimplot_integrated.png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated_subset,
        label.size = 4,
        pt.size = 0.5,
        reduction = paste0(reduction,"SIG",temp_PC_number),
        label = T
      )
      print(temp_dimplot)
      dev.off()
      
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
      
      # FEATUREPLOT
      DefaultAssay(seurat_10X_integrated_subset) <- "RNA"
      temp_png_function(paste0("monocyte_3_cluster/02_",reduction,"_featureplot.png"))
      temp_featureplot <- FeaturePlot(
        object = seurat_10X_integrated_subset,
        features = c("FCGR3A", "CDKN1C",	"CD79B",	"TCF7L2",	"LINC01272",	"RP11-362F19.1",	"PTP4A3",	"LST1",	"PECAM1",	"FCGR3A",	"C11orf21",	"IFITM2",	"CFD",	"PPM1N",	"CD52",	"TESC",	"CLEC12A"), 
        pt.size = 0.25,
        reduction = paste0(reduction,"SIG",temp_PC_number),
        min.cutoff = 0,
        order = T)
      print(temp_featureplot)
      dev.off()
      
      
    }
    
    seurat_10X_integrated_subset <- RenameIdents(seurat_10X_integrated_subset,
                                          `0` = "Monocyte 2",
                                          `1` = "Monocyte 2",
                                          `2` = "Monocyte 2",
                                          `3` = "Monocyte 3",
                                          `4` = "Monocyte 2")
    seurat_10X_integrated_subset@meta.data$celltype_subset <- seurat_10X_integrated_subset@active.ident
    
    
    # plot annotated subset object
    for(reduction in c("UMAP")){
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 6,
            height = 6,
            res = 300,
            units = 'in'
          )
        }
      
      temp_png_function(paste0("monocyte_3_cluster/03_",reduction,"_dimplot_integrated.png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated_subset,
        label.size = 4,
        pt.size = 0.5,
        reduction = paste0(reduction,"SIG",temp_PC_number),
        label = T
      )
      print(temp_dimplot)
      dev.off()
    }
    
    # REANNOTATE
    temp_df <- 
      seurat_10X_integrated@meta.data
    temp_df <- 
      temp_df[,colnames(temp_df) %in% "celltype_subset",drop=F]
    
    temp_df_two <- 
      seurat_10X_integrated_subset@meta.data
    temp_df_two <- 
      temp_df_two[,colnames(temp_df_two) %in% "celltype_subset",drop=F]
    
    temp_df <- 
      temp_df[!rownames(temp_df) %in% rownames(temp_df_two),,drop=F]
    
    temp_df <- rbind(temp_df,
                     temp_df_two)  
    temp_df <- 
      temp_df[rownames(seurat_10X_integrated@meta.data),,drop=F]
      
    print(all.equal(rownames(temp_df),
                    rownames(seurat_10X_integrated@meta.data)))
    
    seurat_10X_integrated@meta.data$celltype_subset <-
      temp_df$celltype_subset
    
    n <- paste0("seurat_10X_integrated_",type,"_filtered")
    assign(n, seurat_10X_integrated)
    
  }
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
  for(metadata in c("celltype_minor", "celltype_subset", temp_colname)){
  temp_png_function(paste0("01_ANNOTATED_", metadata, ".png"))
  temp_dimplot <- DimPlot(
    object = seurat_10X_integrated,
    label.size = 4,
    pt.size = 0.5,
    repel = T,
    reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
    label = T,
    group.by = metadata
  )
  print(temp_dimplot)
  dev.off()

  
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 12,
        height = 6,
        res = 300,
        units = 'in'
      )
    }
  temp_png_function(paste0("02_ANNOTATED_splitbysubtype_", metadata, ".png"))
  temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      pt.size = 0.5,
      reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
      split.by = "subtype",
      label = T,
      label.size = 4,
      group.by = metadata
    )
    print(temp_dimplot)
    dev.off()
    }
}



# DEG ---------------------------------------------------------------------

if(! temp_cell_type == "Plasma_cells" ){
  
  dir.create("03_DGE")
  seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
  
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  Idents(seurat_10X_integrated) <- "celltype_subset"
  
  # remove unassigned
  temp_ids <- unique(seurat_10X_integrated@meta.data$celltype_subset)
  temp_ids <- as.vector(temp_ids[!temp_ids %in% grep("Unassigned", temp_ids,value = T)])
  seurat_10X_integrated <- subset(seurat_10X_integrated,
                                  idents=temp_ids)
  
  temp_cluster.allmarkers <- FindAllMarkers(
    only.pos = T,
    object = seurat_10X_integrated,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    min.diff.pct = 0.1,
    test.use = 'MAST',
    print.bar = F, 
    verbose = F
  )
  temp_cluster.allmarkers <- arrange(temp_cluster.allmarkers,
                                     (cluster),
                                     desc(avg_logFC))
  
  write.csv(temp_cluster.allmarkers,
            paste0("03_DGE/01_MAST_FINDALLMARKERS_RNA_annotated_celltype_unassignedremoved.csv"))

  temp_genes_for_heatmap <- 
    (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(10,
                                                             avg_logFC))
  
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 16,
        height = 18,
        res = 300,
        units = 'in'
      )
    }
  
  seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated,
                                     verbose = FALSE,
                                     vars.to.regress=c("nCount_RNA", "nFeature_RNA", "percent.mito"),
                                     features = as.vector(temp_genes_for_heatmap$gene))
  
  temp_png_function(paste0("03_DGE/02_MAST_FINDALLMARKERS_RNA_annotated_celltype.png"))
  print(DoHeatmap(
    object = seurat_10X_integrated,
    features = temp_genes_for_heatmap$gene,
    group.by = "celltype_subset",
    raster=F,
    size=8,
    angle = 90)
  )
  dev.off()
  
}


# VLNPLOTS/FEATUREPLOTS OF CITESEQ ----------------------------------------------------------------

    dir.create(paste0("04_CITE_markers_raw"))
    dir.create(paste0("04_CITE_markers_raw/all_markers"))

    seurat_10X_CITEimputed <- 
      get(paste0("seurat_10X_integrated_",type,"_filtered"))
    
    DefaultAssay(seurat_10X_CITEimputed) <- "ADT"
    
    # all markers    
    temp_CITE_markers <- 
      row.names(seurat_10X_CITEimputed$ADT@data)
    
    temp_length <- ceiling(length(temp_CITE_markers)/12)
    
    for(i in c(1:temp_length)) {
      #vlnplot
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 16,
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
      
      temp_png_function(paste0("04_CITE_markers_raw/all_markers/01_imputed_CITE_",i,
                               ".png"))
      
      temp_vlnplot <- VlnPlot(object = seurat_10X_CITEimputed, 
                              features = temp_marker_subset, 
                              pt.size = 0, 
                              group.by = "celltype_subset", 
                              log = F)
      print(temp_vlnplot)
      dev.off()
      
      # featureplot
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

      temp_png_function(paste0("04_CITE_markers_raw/all_markers/02_imputed_CITEfeatureplot_",i,
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
    
    # markers of interest
    # temp_marker_subset <- c()
    # 
    # temp_png_function(paste0("04_CITE_markers_raw/02_markers_of_interest_vlnplot.png"))
    # temp_vlnplot <- VlnPlot(object = seurat_10X_CITEimputed, 
    #                         features = temp_marker_subset, 
    #                         pt.size = 0, 
    #                         group.by = "celltype_subset", 
    #                         log = F)
    # print(temp_vlnplot)
    # dev.off()
    # 
    # temp_png_function(paste0("04_CITE_markers_raw/03_markers_of_interest_featureplot.png"))
    # temp_featureplot <- FeaturePlot(
    #   object = seurat_10X_CITEimputed,
    #   features = temp_marker_subset,
    #   pt.size = 0.01,
    #   reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
    #   min.cutoff = 0,
    #   order = T)
    # print(temp_featureplot)
    # dev.off()
    
    
# PROPORTIONS ----------------------------------------------------------

dir.create("05_PROPORTIONS")

library(reshape2)
library(ggpubr)
library(cowplot)

    temp_colname <- "celltype_subset"

    temp_subset <- 
      get(paste0("seurat_10X_integrated_",type,"_filtered"))
    Idents(temp_subset) <- temp_colname
    
    temp_subset@meta.data$subtype <- factor(temp_subset@meta.data$subtype,
                                            levels = rev(c("ER+", "HER2+", "TNBC")))
    
    DefaultAssay(temp_subset) <- "RNA"
    
    # remove sampleIDs with less than 30 cells from stats
    # print(table(temp_subset@meta.data$orig.ident))
    # temp_df <- as.data.frame(table(temp_subset@meta.data$orig.ident))
    # temp_df_filtered <- temp_df[temp_df$Freq > 30,]
    # temp_df_filtered <- as.vector(temp_df_filtered$Var1)
    # temp_removed <- as.vector(temp_df$Var1[!temp_df$Var1 %in% temp_df_filtered])
    # 
    # if(length(temp_removed) > 0){
    #   print("trimming samples")
    #   Idents(temp_subset) <- "orig.ident"
    #   temp_subset <- SubsetData(temp_subset,
    #                             ident.use=temp_df_filtered)
    #   print(temp_removed)
    # }
    
    # estimate proportions
    temp_cluster_names <- 
      unique(temp_subset@meta.data[,temp_colname])
    
    temp_Names <- 
      as.vector(unique(temp_subset@meta.data$orig.ident))
    
    # temp_cluster_names <- factor(temp_cluster_names,
    #                              levels=(0:length(temp_cluster_names)))
    # temp_cluster_names <- sort(temp_cluster_names)
    
    for(i in c(1:length(temp_Names))) {
      # print(i)
      # print(temp_Names[i])
      temp_filtered_summary <- 
        subset(temp_subset@meta.data,
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
    
    
    temp_cellprop_df_all$subtype <- factor(temp_cellprop_df_all$subtype,
                                           levels = rev(c("ER+", "HER2+", "TNBC")))
    # temp_cellprop_df_all$cell_type <- factor(temp_cellprop_df_all$cell_type,
    #                                          levels=(0:length(temp_cluster_names)))
    
    
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
                                   scales ="free_y") +
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
    
    
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 6,
            height = 15,
            res = 300,
            units = 'in'
          )
        }
          
    temp_png_function(paste0("05_PROPORTIONS/01_proportions_combined.png"))
    print(temp_gene_boxplot)
    dev.off()
    
    
    # temp_pdf_function(paste0("07_PROPORTIONS/02_proportions_combined_",cluster,".pdf"))
    # print(temp_gene_boxplot)
    # dev.off()
    


# PROPORTIONS PER PATIENT PER CLUSTER  ------------------------------------

    temp_colname <- "celltype_subset"
    Idents(seurat_10X_integrated) <- temp_colname
    
    temp_df <- seurat_10X_integrated@meta.data
    temp_df_combined <- NULL
    for(cluster in unique(seurat_10X_integrated@meta.data$SUBSET30_res.1)){
      temp_df_subset <- temp_df[temp_df$SUBSET30_res.1 == cluster,,drop=F]
      temp_table <- as.data.frame(table(temp_df_subset$orig.ident))
      temp_table$proportion <- temp_table[,2]/sum(temp_table[,2])
      temp_table$cluster <- cluster
      temp_table$celltype_subset <- as.vector(unique(temp_df_subset$celltype_subset))   
      temp_table <- temp_table[order(-temp_table$Freq),,drop=F]
      
      temp_df_combined <- rbind(temp_df_combined,
                                temp_table)
      

      rm(temp_table)
      
    }
    
    write.csv(temp_df_combined,
              "05_PROPORTIONS/02_cluster_proportions_per_patient.csv")
    
      
    # plot
    temp_png_function <-
      function(x) {
        png(
          file = (x), 
          width = 6, 
          height = 10, 
          res = 300, 
          units = 'in'
        )
      }
    
    temp_df_combined$cluster <- factor(temp_df_combined$cluster,
                                       levels=rev(str_sort(unique(temp_df_combined$cluster),numeric=T)))
    
    temp_ggplot <- ggplot(temp_df_combined, 
                          aes(x=cluster,
                              fill=Var1,
                              y=proportion)) +
      geom_bar(stat="identity", 
               position="fill") +
      theme(axis.text.x=element_blank(),
            legend.title=element_text(size=10), 
            legend.text=element_text(size=7)) +
      scale_y_continuous(position = "right") + 
      xlab("Cluster") +
      ylab("Proportion per patient") +
      guides(fill=guide_legend(title="Cell Type")) +
      coord_flip()
    
    
    temp_png_function(paste0("05_PROPORTIONS/03_cluster_proportions_per_patient.png"))
    print(temp_ggplot)
    dev.off()
    
    
    temp_ggplot <- ggplot(temp_df_combined, 
                          aes(x=celltype_subset,
                              fill=Var1,
                              y=proportion)) +
      geom_bar(stat="identity", 
               position="fill") +
      theme(axis.text.x=element_blank(),
            legend.title=element_text(size=10), 
            legend.text=element_text(size=7)) +
      scale_y_continuous(position = "right") + 
      xlab("Cell type") +
      ylab("Proportion per patient") +
      guides(fill=guide_legend(title="Cell Type")) +
      coord_flip()
    
    
    temp_png_function(paste0("05_PROPORTIONS/04_celltype_proportions_per_patient.png"))
    print(temp_ggplot)
    dev.off()
    
    
    
    
        
# SAVE OBJECTS -------------------------------------------------

seurat_10X_integrated <- get(paste0("seurat_10X_integrated_",type,"_filtered"))
saveRDS(seurat_10X_integrated, 
        paste0("RDATA_01_PROCESSED_FILTERED_IMPUTEDCITE_RAW.Rdata"))



# RECLUSTERING SEURAT OBJECTS FOR STROMAL CELLS
# Sunny Z. Wu
# 
#
# Run in screen 
# qrsh -pe smp 24 -l mem_requested=10G -P TumourProgression
# source activate r_spatial
# R
#
# 01: SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(future)
library(data.table)

# DIRECTORIES
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Jun2019/04_reclustering_analysis/analysis01_mesenchymal_subtype_normals/")

dir.create("output_run02")
setwd("output_run02")


# 02: LOAD DATASET ---------------------------------------------------------------

# load UNFILTERED integrated seurat object
seurat_10X <- readRDS("/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/Jun2019_final_primary_set/03_integration/output/CCA/Output/Rdata/03_seurat_CCA_aligned_processed.Rdata")

# 03: SUBSET AND SPLIT OBJECT --------------------------------------------------------

Idents(seurat_10X) <- paste0("garnett_ext")
seurat_10X <- subset(seurat_10X,
                      idents=c("CAFs", "Endothelial"))

temp_sample_list <- SplitObject(seurat_10X, 
                                split.by = "orig.ident")

# set k.filter to the smallest sample
temp_k_filter <- print(min(sapply(temp_sample_list, ncol)))

for (i in 1:length(temp_sample_list)) {
  
  print(unique(temp_sample_list[[i]]@meta.data$orig.ident))
  DefaultAssay(temp_sample_list[[i]]) <- "RNA"
  
  temp_sample_list[[i]]@assays$integrated <- NULL
  
  # remove integrated assay to avoid this subsequent error during FindIntegrationANchors
  # Error in RowMergeMatrices(mat1 = mat1, mat2 = mat2, mat1_rownames = mat1.names,  : Need S4 class dgRMatrix for a sparse matrix
  # issue here but no solution https://github.com/satijalab/seurat/issues/1355
  
  temp_sample_list[[i]] <- 
    NormalizeData(temp_sample_list[[i]], 
                  verbose = F)
  
  temp_sample_list[[i]] <- 
    FindVariableFeatures(temp_sample_list[[i]], 
                         verbose = F,
                         do.plot = F,
                         nfeatures = 5000,
                         selection.method="dispersion")
  
  # temp_sample_list[[i]] <- ScaleData(
  #   object = temp_sample_list[[i]],
  #   display.progress = F
  # )
  
  print(length(temp_sample_list[[i]]@assays$RNA@var.features))
  print(length(temp_sample_list[[i]]@assays$RNA@scale.data))
  
}


# 04: RECLUSTERING (IN PARALLEL) ------------------------------------------

# future params
plan("multiprocess", 
     workers = 12)
options(future.globals.maxSize = 10 * 1024^3)

print(Sys.time())
seurat_10X_anchors <- FindIntegrationAnchors(object.list = temp_sample_list, 
                                             verbose = F, 
                                             k.filter = temp_k_filter,
                                             anchor.features=5000)
print(Sys.time())

print(Sys.time())
seurat_10X_integrated <- IntegrateData(anchorset = seurat_10X_anchors, 
                                       verbose = F)
print(Sys.time())


DefaultAssay(object = seurat_10X_integrated) <- 
  "integrated"

# re-scale data
seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated, 
                                   verbose = FALSE)

# PCA 
plan("sequential")
seurat_10X_integrated <- RunPCA(object = seurat_10X_integrated, 
                                npcs = 100, 
                                verbose = T)

# dimensional reduction
# UMAP
seurat_10X_integrated <- 
  RunUMAP(seurat_10X_integrated, 
          dims =  1:30,
          reduction.key = paste0("UMAPCOMBINEDFILTERED",30,"_"),
          reduction.name = paste0("UMAPCOMBINEDFILTERED",30),
          verbose = F
  )

# TSNE
seurat_10X_integrated <-
  RunTSNE(object = seurat_10X_integrated,
          dims = 1:30,
          reduction.key = paste0("TSNECOMBINEDFILTERED",30,"_"),
          reduction.name = paste0("TSNECOMBINEDFILTERED",30)
  )

# find neighbours step 
seurat_10X_integrated <- 
  FindNeighbors(object = seurat_10X_integrated,
                dims = 1:30,
                graph.name = paste0("COMBINEDFILTERED", 30)
  )

# clustering 
seurat_10X_integrated <-
  FindClusters(
    object = seurat_10X_integrated,
    graph.name = paste0("COMBINEDFILTERED", 30),
    resolution = c(0.4,0.6,0.8,1.2)
  )




# 05: PLOT MARKERS --------------------------------------------------------

temp_markers <-
  c("COL1A1","PDGFRA", "FAP", "PDPN", "CXCL12", "PDGFRB", "ACTA2",
    "CD36", "MCAM", "MYH11", "RGS5",
    "PECAM1", "VWF", "CD34", "ACKR1", "DLL4", "LYVE1", "MKI67")

for(reduction in c("TSNE", "UMAP")){
  print(reduction)
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
  
  temp_colname <- paste0("COMBINEDFILTERED", 30, "_res.0.4")
  temp_reduction <- paste0(reduction,"COMBINEDFILTERED",30)
  
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  temp_png_function(paste0("01_ALL_STROMAL_CELLS/",reduction,"_01_featureplot.png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = temp_markers,
    pt.size = 0.1,
    reduction = temp_reduction,
    order = T)
  
  print(temp_featureplot)
  dev.off()
  
  # dimplot tSNE
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 10, 
        height = 10, 
        res = 450, 
        units = 'in'
      )
    }
  
  temp_dimplot <- DimPlot(
    object = seurat_10X_integrated,
    label.size = 6,
    pt.size = 0.1,
    label = T,
    reduction = temp_reduction,
    group.by = temp_colname 
  )
  
  # get legend
  library(ggpubr)
  leg <- get_legend(temp_dimplot)
  leg <- as_ggplot(leg)
  
  
  # remove legend
  temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
    xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))
  
  temp_png_function(paste0("01_ALL_STROMAL_CELLS/",reduction,"_02_dimplot.png"))
  print(temp_dimplot)
  dev.off()
  
  # print legend
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x), 
        width = 4, 
        height = 3, 
        useDingbats = F
      )
    }
  temp_pdf_function("01_ALL_STROMAL_CELLS/03_legend.pdf")
  print(leg)
  dev.off()
}    


# 06: ANNOTATE MAJOR LINEAGES ---------------------------------------------

Idents(seurat_10X_integrated) <- paste0("COMBINEDFILTERED", 30, "_res.0.4")

seurat_10X_integrated <- RenameIdents(seurat_10X_integrated,
                                               `0` = "Endothelial",
                                               `1` = "CAFs",
                                               `2` = "PVL",
                                               `3` = "CAFs",
                                               `4` = "Endothelial",
                                               `5` = "PVL",
                                               `6` = "Endothelial",
                                               `7` = "Lymphatic_endothelial",
                                               `8` = "Cycling")

# 07: RECLUSTER EACH MAJOR CELL TYPE ---------------------------------------------------------------------

temp_colname <- paste0("COMBINEDFILTERED", 30, "_res.0.4")
Idents(seurat_10X_integrated) <- temp_colname

for(cluster in c("CAFs", "PVL", "Endothelial")){
  
  temp_subset <- subset(seurat_10X_integrated,
                        idents=cluster)
  
  DefaultAssay(object = temp_subset) <- 
    "integrated"
  
  # re-scale data
  temp_subset <- ScaleData(object = temp_subset, 
                           verbose = FALSE)
  
  # PCA 
  plan("sequential")
  temp_subset <- RunPCA(object = temp_subset, 
                        npcs = 20, 
                        verbose = T)
  
  # dimensional reduction
  # UMAP
  temp_subset <- 
    RunUMAP(temp_subset, 
            dims =  1:20,
            reduction.key = paste0("UMAPFILTEREDSUBSET",20,"_"),
            reduction.name = paste0("UMAPFILTEREDSUBSET",20),
            verbose = F
    )
  
  # TSNE
  temp_subset <-
    RunTSNE(object = temp_subset,
            dims = 1:20,
            reduction.key = paste0("TSNEFILTEREDSUBSET",20,"_"),
            reduction.name = paste0("TSNEFILTEREDSUBSET",20)
    )
  
  # reclust 
  temp_subset <- 
    FindNeighbors(object = temp_subset,
                  dims = 1:20,
                  graph.name = paste0("FILTEREDSUBSET", 20)
    )
  
  # default clustering at 0.8
  temp_subset <-
    FindClusters(
      object = temp_subset,
      graph.name = paste0("FILTEREDSUBSET", 20),
      resolution = c(0.4,0.6,0.8,1.2)
    )
  
  assign(paste0("temp_subset_",cluster),
         temp_subset)
}



# SAVE RDS ----------------------------------------------------------------

dir.create("RData")
for(cluster in c("CAFs", "PVL", "Endo")){

  temp_subset <- get(paste0("temp_subset_",cluster))
  saveRDS(temp_subset,
          paste0("RData/RData_",cluster,".Rdata"))
}

saveRDS(seurat_10X_integrated,
        paste0("RData/RData_stromal_combined.Rdata"))


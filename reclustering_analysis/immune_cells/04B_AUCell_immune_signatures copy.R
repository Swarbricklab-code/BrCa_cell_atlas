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

# temp_cell_type <- "B_cells"
# temp_wd <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/run02/analysis_01_B_cells/"



# SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(Seurat)
library(GSEABase)
library(AUCell)
library(reshape2)
library(NMF)
library(Matrix)
library(cowplot)

# DIRECTORY ---------------------------------------------------------------

setwd(paste0(temp_wd))

dir.create("04_output_AUCell")
setwd("04_output_AUCell")

# LOAD DATA ---------------------------------------------------------------

for(type in c("integrated")) {
  
  seurat_10X_integrated <- readRDS(paste0("../03_output_filtered_plotting/RDATA_01_PROCESSED_FILTEREDintegrated.Rdata"))
  
  n <- paste0("seurat_10X_integrated_",type,"_filtered")
  assign(n, seurat_10X_integrated)
}

# read genesets
temp_genesets <- 
  "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/04_reclustering_analysis/01_immune_reclustering/combined_immune_signatures.csv"

  temp_xcell_genesets <- 
    read.csv(temp_genesets)
  
  temp_GeneSetCollection_list <- NULL
  
  for(i in c(1:ncol(temp_xcell_genesets))) {
    n <- paste0("temp_set_name",
                i)
    assign(n,
           colnames(temp_xcell_genesets[i]))
    
    temp_set_name <- get(paste0("temp_set_name", 
                                i))
    
    temp_set <- na.omit(temp_xcell_genesets[i])
    
    colnames(temp_set) <- "gene_set"
    
    temp_set <- GeneSet(unique(as.vector(temp_set$gene_set)),
                        setName = temp_set_name)
    temp_GeneSetCollection_list <- append(temp_GeneSetCollection_list,
                                          temp_set)
    
  }
  
  temp_gene_set_collection <- 
    GeneSetCollection(temp_GeneSetCollection_list)
  
  rm(list = ls(pattern = "temp_set_name"))


# EXPORT MATRIX FOR AUCELL ------------------------------------------------------------------

  temp_run_by_sampleID <- T
  
  if(temp_run_by_sampleID){
    temp_sampleID <- 
      unique(seurat_10X_integrated_integrated_filtered@meta.data$orig.ident)
    
    Idents(seurat_10X_integrated_integrated_filtered) <- "orig.ident"
    
    for(sampleID in temp_sampleID) {
      print(sampleID)
      
      temp_seurat_subset <-
        subset(seurat_10X_integrated_integrated_filtered, idents=sampleID)
      
      temp_exprMatrix <- GetAssayData(object = temp_seurat_subset, 
                                      assay = "RNA", 
                                      slot = "data")
      
      temp_exprMatrix <- 
        Matrix(temp_exprMatrix,
               sparse = T)
      
      n <- paste0("temp_exprMatrix_",sampleID)
      assign(n, temp_exprMatrix)
    }
    rm(temp_seurat_subset,temp_exprMatrix)
  }
  

# AUCELL ------------------------------------------------------------------

  if(temp_run_by_sampleID){
    
    for(sampleID in temp_sampleID) {
      print(sampleID)
      temp_exprMatrix <- get(paste0("temp_exprMatrix_",sampleID))
      
      dim(temp_exprMatrix)
      
      print("Buildrankings")
      temp_cells_rankings <- 
        AUCell_buildRankings(temp_exprMatrix, 
                             nCores = 1, 
                             plotStats = F)
      print("subset genes")
      # subset gene sets
      temp_subsetgeneSets <- 
        subsetGeneSets(temp_gene_set_collection, 
                       rownames(temp_exprMatrix)) 
      
      print("CalcAUC")
      # calculate area under the curve
      temp_cells_AUC <- 
        AUCell_calcAUC(geneSets = temp_subsetgeneSets, 
                       rankings = temp_cells_rankings, 
                       aucMaxRank = ceiling(0.05 * nrow(temp_cells_rankings)), 
                       nCores = 1, 
                       verbose = T)
      
      # save(temp_cells_AUC, 
      #      file=paste0("output/",sampleID,"_01_cells_AUC_aucMaxRank.RData"))
      
      print("Transposing matrix")
      #transpose matrix for seurat metadata assignment
      temp_cells_AUC_matrix <- 
        t(as.data.frame(getAUC(temp_cells_AUC)))
      
      # saveRDS(temp_cells_AUC_matrix, 
      #         paste0("output/",sampleID,"_02_cells_AUC_matrix_transposed.Rdata"))
      
      n <- paste0("temp_cells_AUC_matrix_",sampleID)
      assign(n,temp_cells_AUC_matrix)
      print(paste0(sampleID, " done"))
    }
    rm(temp_cells_AUC_matrix,temp_cells_AUC,temp_subsetgeneSets)
  }
  
# COMBINED AND ADD TO SEURAT OBJECT ---------------------------------------

  if(temp_run_by_sampleID){
    
    temp_cells_AUC_matrix <- 
      NULL
    
    for(sampleID in temp_sampleID) {
      print(sampleID)
      temp_cells_AUC_matrix_subset <- 
        get(paste0("temp_cells_AUC_matrix_",sampleID))
      
      print(dim(temp_cells_AUC_matrix_subset))
      
      temp_cells_AUC_matrix <-
        rbind(temp_cells_AUC_matrix,
              temp_cells_AUC_matrix_subset)
      
      print(dim(temp_cells_AUC_matrix))
      rm(temp_cells_AUC_matrix_subset)
      
    }
  }
  
  temp_cells_AUC_matrix_sorted <-
    temp_cells_AUC_matrix[rownames(seurat_10X_integrated_integrated_filtered@meta.data),,drop=FALSE]
  
  temp_cells_AUC_matrix_sorted <- 
    as.data.frame.matrix(temp_cells_AUC_matrix_sorted)
  
    temp_first_geneset <- 
      (ncol(seurat_10X_integrated_integrated_filtered@meta.data) + 1)
    
    seurat_10X_integrated_integrated_filtered <- AddMetaData(seurat_10X_integrated_integrated_filtered, 
                              metadata = temp_cells_AUC_matrix_sorted)
    
    saveRDS(seurat_10X_integrated_integrated_filtered,
            "RDATA_04_PROCESSED_FILTERED_RESCALED_AUCell.Rdata")


# PLOTTING ----------------------------------------------------------------
      
      dir.create("01_featureplots_vlnplots")
      
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
    
    temp_png_function_2 <- 
      function(x) {
        png(
          file = (x), 
          width = 6, 
          height = 5, 
          res = 300, 
          units = 'in'
        )
      }
    
    
    
      temp_last_geneset <-
        (ncol(seurat_10X_integrated_integrated_filtered@meta.data))
      
      temp_gene_set_names <-
        colnames(seurat_10X_integrated_integrated_filtered@meta.data[temp_first_geneset:temp_last_geneset])
      
      for(i in c(1:length(temp_gene_set_names))) {
        
        temp_gene_set_name <- 
          (temp_gene_set_names[i])
        print(temp_gene_set_name)
        
        temp_png_function(paste0("01_featureplots_vlnplots/",
                                 i, "_featureplot",
                                 temp_gene_set_name,
                                 ".png"))
        temp_featureplot <- FeaturePlot(
          object = seurat_10X_integrated_integrated_filtered,
          features = temp_gene_set_name,
          order = T,
          pt.size = 0.5,
          reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
          min.cutoff = "q50")
        print(temp_featureplot)
        dev.off()
        
        temp_png_function_2(paste0("01_featureplots_vlnplots/",
                                 i, "_vlnplot",
                                 temp_gene_set_name,
                                 ".png"))
        temp_vlnplot <- VlnPlot(
          object = seurat_10X_integrated_integrated_filtered,
          features = temp_gene_set_name,
          pt.size = 0,
          group.by = paste0("SUBSET", temp_PC_number, "_res.",temp_res),
          log = F)
        print(temp_vlnplot)
        dev.off()
      }
      
      

# FEATUREPLOT SPLIT BY PATIENT --------------------------------------------

      dir.create("02_featureplots_by_patient")
      
      temp_png_function <- 
        function(x) {
          png(
            file = (x), 
            width = 8, 
            height = 40, 
            res = 300, 
            units = 'in'
          )
        }
    
      
      for(i in c(1:length(temp_gene_set_names))) {
        
        temp_gene_set_name <- 
          (temp_gene_set_names[i])
        print(temp_gene_set_name)
        Idents(seurat_10X_integrated_integrated_filtered) <- "orig.ident"
        
        temp_median_visual_cutoff <- median(seurat_10X_integrated_integrated_filtered@meta.data[,temp_gene_set_name])
        
        # for(sampleID in unique(seurat_10X_integrated_integrated_filtered@meta.data$orig.ident)){
        #   print(sampleID)
        #   temp_subset <- subset(seurat_10X_integrated_integrated_filtered,
        #                         idents=sampleID)
          
          temp_featureplot <- FeaturePlot(
            object = seurat_10X_integrated_integrated_filtered,
            features = temp_gene_set_name,
            order = T,
            pt.size = 1,
            reduction = paste0(temp_reduction_use,"SIG",temp_PC_number),
            min.cutoff = "q50", 
            by.col = F,
            split.by = "orig.ident", 
            coord.fixed = T
            )
        #   temp_featureplot <- 
        #     temp_featureplot + labs(title = sampleID)
        # 
        # assign(paste0("temp_featureplot_",sampleID),
        #        temp_featureplot)
        # }
        
        # temp_grid <-
        #   plot_grid(plotlist=mget(paste0("temp_featureplot_",unique(seurat_10X_integrated_integrated_filtered@meta.data$orig.ident))),
        #             ncol = 5)
        # temp_png_function("temp.png")
        temp_png_function(paste0("02_featureplots_by_patient/",
                                 i, "_featureplot_by_patient",
                                 temp_gene_set_name,
                                 ".png"))
        print(temp_featureplot)
        dev.off()
              
        
        }
      
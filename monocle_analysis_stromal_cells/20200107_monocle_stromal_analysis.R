# Monocle analysis on stromal cells
# SUNNY WU
#
# 01: SETUP -------------------------------------------------------------------

library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(grid)

# 02: DIRECTORIES -------------------------------------------------------------

dir.create("MONOCLE_ANALYSIS")
setwd("MONOCLE_ANALYSIS")

# 03: LOAD DATA ------------------------------------------------------------

for(cluster in c("CAFs", "PVL", "Endo")){
  temp_subset <- readRDS(paste0("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/05_figures/05_STROMAL/RData/RData_",cluster,".Rdata"))
  n <-  paste0("temp_subset_",cluster)
  assign(n,temp_subset)
}


res <- 0.8
npcs <- 20



# 04: MONOCLE 2 PROCESSING ---------------------------------------------------------------

dir.create("monocle_new_ordering")

temp_colname <- paste0("FILTEREDSUBSET", npcs, "_res.",res)
Idents(temp_subset) <- temp_colname

library(monocle)
for(cluster in c("CAFs", "PVL", "Endo")){
  
  if(cluster == "CAFs"){
    temp_ngenes_to_use <- 1000
  }
  if(cluster == "PVL"){
    temp_ngenes_to_use <- 500
  }
  if(cluster == "Endo"){
    temp_ngenes_to_use <- 1000
  }
  print(cluster)
  dir.create(paste0("monocle_new_ordering/",cluster))
  
  temp_subset <- get(paste0("temp_subset_",cluster))
  
  #Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(temp_subset@assays$integrated@data), 'sparseMatrix')
  
  pd <- new('AnnotatedDataFrame', data = temp_subset@meta.data)
  
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  #Construct monocle cds
  HSMM <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily = uninormal())
  
  #Run ordering algorithm
   print(" running DEG")
   diff_test_res <-
     differentialGeneTest(HSMM,
                          fullModelFormulaStr = paste0("~",temp_colname),
                          cores = 4)
  
   diff_test_res_filtered <- 
     diff_test_res[diff_test_res$qval < 0.001,,drop=F]
   
   diff_test_res_filtered <-
     diff_test_res_filtered[order(diff_test_res_filtered$qval),]
  
   write.csv(diff_test_res_filtered,
             paste0("monocle_new_ordering/",
                    cluster,
                    "/ordering_genes_",temp_ngenes_to_use,"_used.csv"))
   
   temp_ordering_genes <-
     row.names(diff_test_res)[order(diff_test_res$qval)][1:temp_ngenes_to_use]
   

  # var_genes <- temp_subset[["integrated"]]@var.features[1:2000]
  ordering_genes <- temp_ordering_genes
  
  HSMM <- setOrderingFilter(HSMM, 
                            ordering_genes)
  print(dim(exprs(HSMM)))
  
  ## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
  HSMM <- reduceDimension(HSMM,
                          norm_method="none", 
                          reduction_method="DDRTree")
  
  HSMM <- orderCells(HSMM)
  
  n <- paste0("HSMM_",cluster)
  assign(n, HSMM)
}

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

# plot 
temp_diff_genes_to_plot <- c("ACTA2", "TAGLN", "CXCL12", "PECAM1", 
                             "ACKR1", "DLL4", "COL1A1", "ALDH1A1", 
                             "ZEB1", "RGS5", "CD36", "MYH11")

for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  HSMM <- get(paste0("HSMM_",cluster))
  # First decide what you want to color your cells by
  # print(head(pData(HSMM)))
  
  temp_diff_genes_to_plot_filtered <- 
    temp_diff_genes_to_plot[temp_diff_genes_to_plot %in% rownames(HSMM)]
  
  ## order cells change colors and theta to match your plot
  for(type in c(temp_colname, "Pseudotime", "orig.ident", "State", "subtype")) {
    print(type)
    temp_ggplot <- plot_cell_trajectory(HSMM, 
                                        color_by = type,
                                        show_branch_points = T,
                                        show_tree = TRUE,
                                        cell_size = 0.5) 
    temp_png_function(paste0("monocle_new_ordering/",cluster,"/01_plot_cell_trajectory_",type,".png"))
    print(temp_ggplot)
    dev.off()
    if(! type == "Pseudotime"){
      HSMM_new <- HSMM
      colnames(pData(HSMM_new))[which(names(pData(HSMM_new)) == type)] <- "newcolname"
      
      temp_png_function(paste0("monocle_new_ordering/",cluster,"/01_plot_cell_trajectory_",type,"_facet.png"))
      temp_ggplot <- plot_cell_trajectory(HSMM_new, 
                                          color_by = "newcolname",
                                          show_branch_points = T,
                                          show_tree = TRUE,
                                          cell_size = 0.5) +
        facet_wrap(~newcolname, nrow = 2)
      print(temp_ggplot)
      dev.off()
    }
  }
  
  temp_png_function(paste0("monocle_new_ordering/",cluster, "/02_plot_genes_in_pseudotime.png"))
  temp_plot <-   
    plot_genes_in_pseudotime(HSMM[temp_diff_genes_to_plot_filtered,], 
                             color_by = "State",
                             ncol=3)
  print(temp_plot)
  dev.off()
  
  temp_png_function(paste0("monocle_new_ordering/",cluster, "/03_plot_genes_in_branched_pseudotime.png"))
  temp_plot <-   
    plot_genes_branched_pseudotime(HSMM[temp_diff_genes_to_plot_filtered,], 
                                   color_by = "State",
                                   ncol=3)
  print(temp_plot)
  dev.off()
  
}

# 05: ASSIGN ROOT STATE  ----------------------------------------------

for(cluster in c("CAFs", "PVL", "Endo")){
  if(cluster == "CAFs"){
    rootstate <- 4
  }
  if(cluster == "PVL"){
    rootstate <- 1
  }
  if(cluster == "Endo"){
    rootstate <- 3
  }
  print(cluster)
  HSMM <- get(paste0("HSMM_",cluster))
  
  # set root state
  HSMM <-
    orderCells(HSMM,
               root_state = rootstate)
  
  temp_png_function(paste0("monocle_new_ordering/",cluster,"/04_plot_cell_trajectory_Pseudotime_newROOTSTATE.png"))
  temp_ggplot <- plot_cell_trajectory(HSMM, 
                                      color_by = "Pseudotime",
                                      show_branch_points = T,
                                      show_tree = TRUE,
                                      cell_size = 0.5) 
  print(temp_ggplot)
  dev.off()
  
  n <- paste0("HSMM_", cluster)
  assign(n,HSMM)
}

# replot
temp_diff_genes_to_plot <- c("ACTA2", "TAGLN", "CXCL12", "PECAM1", 
                             "ACKR1", "DLL4", "COL1A1", "ALDH1A1", 
                             "ZEB1", "RGS5", "CD36", "MYH11")

for(cluster in c("CAFs", "PVL", "Endo")){
  
  HSMM <- get(paste0("HSMM_",cluster))
  
  # First decide what you want to color your cells by
  # print(head(pData(HSMM)))
  
  temp_diff_genes_to_plot_filtered <- 
    temp_diff_genes_to_plot[temp_diff_genes_to_plot %in% rownames(HSMM)]
  
  ## order cells change colors and theta to match your plot
  for(type in c(temp_colname, "Pseudotime", "orig.ident", "State")) {
    temp_ggplot <- plot_cell_trajectory(HSMM, 
                                        color_by = type,
                                        show_branch_points = T,
                                        show_tree = TRUE,
                                        cell_size = 0.5) 
    temp_png_function(paste0("monocle_new_ordering/",cluster,"/05_plot_cell_trajectory_newROOTSTATE_",type,".png"))
    print(temp_ggplot)
    dev.off()
  }
  
  temp_png_function(paste0("monocle_new_ordering/",cluster, "/06_plot_genes_in_pseudotime_newROOTSTATE.png"))
  temp_plot <-   
    plot_genes_in_pseudotime(HSMM[temp_diff_genes_to_plot_filtered,], 
                             color_by = "State",
                             ncol=3)
  print(temp_plot)
  dev.off()
  
  temp_png_function(paste0("monocle_new_ordering/",cluster, "/07_plot_genes_in_branched_pseudotime_newROOTSTATE.png"))
  temp_plot <-   
    plot_genes_branched_pseudotime(HSMM[temp_diff_genes_to_plot_filtered,], 
                                   color_by = "State",
                                   ncol=3)
  print(temp_plot)
  dev.off()
  
}




# 06: MONOCLE DEG PSEUDOTIME --------------------------------------------------

for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  HSMM <- get(paste0("HSMM_",cluster))
  
  to_be_tested <- row.names(fData(HSMM))
  cds_subset <- HSMM[to_be_tested,]
  
  temp_diff_test_res <- differentialGeneTest(HSMM,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)")
  
  temp_diff_test_res <- arrange(temp_diff_test_res,
                                pval)
  
  n <- paste0("temp_diff_test_res_", cluster)
  assign(n,temp_diff_test_res)
  
  write.csv(temp_diff_test_res, 
            paste0("monocle_new_ordering/",cluster, "/08_differentialGeneTest_pseudotime.csv"))
  
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

# temp_colname <- paste0("FILTEREDSUBSET", npcs, "_res.",res)
temp_colname <- "State"
for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  
  HSMM <- 
    get(paste0("HSMM_",cluster))
  temp_diff_test_res <- 
    get(paste0("temp_diff_test_res_", cluster))
  
  temp_diff_test_res <- arrange(temp_diff_test_res,
                                pval)
  
  temp_genes_to_plot <- temp_diff_test_res$gene_short_name[1:20]
  
  temp_genes_to_plot <- row.names(subset(fData(HSMM),
                                         gene_short_name %in% temp_genes_to_plot))
  
  cds_subset <- HSMM[temp_genes_to_plot,]
  
  temp_png_function(paste0("monocle_new_ordering/",cluster, "/09_plot_genes_in_pseudotime.png"))
  temp_plot <-   
    plot_genes_in_pseudotime(cds_subset, 
                             color_by = temp_colname,
                             ncol= 4)
  print(temp_plot)
  dev.off()
  
  temp_genes_to_plot <- temp_diff_test_res$gene_short_name[21:40]
  
  temp_genes_to_plot <- row.names(subset(fData(HSMM),
                                         gene_short_name %in% temp_genes_to_plot))
  
  cds_subset <- HSMM[temp_genes_to_plot,]
  
  temp_png_function(paste0("monocle_new_ordering/",cluster, "/10_plot_genes_in_pseudotime.png"))
  temp_plot <-   
    plot_genes_in_pseudotime(cds_subset, 
                             color_by = temp_colname,
                             ncol= 4)
  print(temp_plot)
  dev.off()
  
}

# 07: MONOCLE EXPLORE BRANCHES ------------------------------------------------

for(cluster in c("CAFs", "PVL", "Endo")){
  
  dir.create(paste0("monocle_new_ordering/",cluster,"/11_branches/"))
  print(cluster)
  
  HSMM <- 
    get(paste0("HSMM_",cluster))
  
  temp_num_branches <- 
    length(HSMM@auxOrderingData$DDRTree$branch_points)
  
  print(paste0("number of branches = ",temp_num_branches))
  
  for(branch in c(1:temp_num_branches)){
    print(paste0("analyzing branch ",branch))
    dir.create(paste0("monocle_new_ordering/",cluster,"/11_branches/branch_",branch))
    
    # BEAM
    # downsample cells to speed up analysis
    # sampleCells <- detectGenes(HSMM)
    # sampleCells_subset <- sampleCells[fData(sampleCells)$num_cells_expressed > 100, ]
    
    # run beam
    print("Running BEAM")
    temp_start <- Sys.time()
    BEAM_res <- BEAM(HSMM, 
                     branch_point = branch, 
                     cores = 4)
    print(paste0("BEAM took ",(Sys.time()-temp_start), " secs"))
    
    BEAM_res <- BEAM_res[order(BEAM_res$qval),]
    write.csv(BEAM_res,
              paste0("monocle_new_ordering/",cluster,"/11_branches/branch_",branch,"/01_BEAM_res.csv"))
    
    BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
    BEAM_res_genes <- as.vector(BEAM_res$gene_short_name[1:50])
    
    # plot
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
    
    temp_png_function(paste0("monocle_new_ordering/",cluster, "/11_branches/branch_",branch,"/02_plot_genes_branched_heatmap.png"))
    temp_ggplot <- plot_genes_branched_heatmap(HSMM[BEAM_res_genes,],
                                               branch_point = branch,
                                               num_clusters = 2,
                                               cores = 4,
                                               use_gene_short_name = T,
                                               show_rownames = T)
    print(temp_ggplot)
    dev.off()
    
    # plot in psuedotime
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
    
    temp_png_function(paste0("monocle_new_ordering/",cluster, "/11_branches/branch_",branch,"/03_plot_genes_branched_pseudotime.png"))
    temp_plot <-
      plot_genes_branched_pseudotime(HSMM[BEAM_res_genes,],
                                     branch_point = branch,
                                     color_by = "State", 
                                     ncol = 5)
    print(temp_plot)
    dev.off()
    
  }
}



# 08: STATE ANNOTATION IN SEURAT OBJECTS ----------------------------------------------------

for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  
  HSMM <- get(paste0("HSMM_",cluster))
  temp_subset <- get(paste0("temp_subset_",cluster))
  
  temp_metadata_df <- pData(HSMM)
  print(all.equal(rownames(temp_metadata_df), rownames(temp_subset@meta.data)))
  temp_subset@meta.data$State <- temp_metadata_df$State
  
  print(table(temp_subset@meta.data$State))
  
  if(cluster == "CAFs"){
    temp_subset@meta.data$State <- factor(temp_subset@meta.data$State,
                                          levels = c(4,3,2,5,1))
    Idents(temp_subset) <- "State"
    temp_subset <- RenameIdents(temp_subset,
                                `1` = "Myofibroblast-like CAFs 2",
                                `2` = "Transitioning CAFs",
                                `3` = "Inflammatory-CAFs 2",
                                `4` = "Inflammatory-CAFs 1",
                                `5` = "Myofibroblast-like CAFs 1")
    temp_subset@meta.data$celltype_subset <- temp_subset@active.ident
    
    Idents(temp_subset) <- "State"
    temp_subset <- RenameIdents(temp_subset,
                                `1` = "Myofibroblast-like CAFs",
                                `2` = "Myofibroblast-like CAFs", # transitioning cells called as differenaited in minor call
                                `3` = "Inflammatory-CAFs",
                                `4` = "Inflammatory-CAFs",
                                `5` = "Myofibroblast-like CAFs")
    temp_subset@meta.data$celltype_minor <- temp_subset@active.ident
    
    temp_major_label <- "CAFs"
    temp_subset@meta.data$celltype_major <- temp_major_label
    
  }
  if(cluster == "PVL"){
    temp_subset@meta.data$State <- factor(temp_subset@meta.data$State,
                                          levels = c(5,6,4,3,2,7,1))
    Idents(temp_subset) <- "State"
    temp_subset <- RenameIdents(temp_subset,
                                `1` = "Differentiated PVL 3",
                                `2` = "Transitioning PVL 2",
                                `3` = "Differentiated PVL 1",
                                `4` = "Transitioning PVL 1",
                                `5` = "Immature PVL 1",
                                `6` = "Immature PVL 2",
                                `7` = "Differentiated PVL 2")
    temp_subset@meta.data$celltype_subset <- temp_subset@active.ident
    
    Idents(temp_subset) <- "State"
    temp_subset <- RenameIdents(temp_subset,
                                `1` = "Differentiated PVL",
                                `2` = "Differentiated PVL",
                                `3` = "Differentiated PVL",
                                `4` = "Differentiated PVL",
                                `5` = "Immature PVL",
                                `6` = "Immature PVL",
                                `7` = "Differentiated PVL")
    temp_subset@meta.data$celltype_minor <- temp_subset@active.ident
    
    temp_major_label <- "PVL"
    temp_subset@meta.data$celltype_major <- temp_major_label
    
  }
  if(cluster == "Endo"){
    temp_subset@meta.data$State <- factor(temp_subset@meta.data$State,
                                          levels = c(3,2,1))
    Idents(temp_subset) <- "State"
    temp_subset <- RenameIdents(temp_subset,
                                `1` = "Tip-like endothelial cells",
                                `2` = "Pericyte-like endothelial cells",
                                `3` = "Stalk-like endothelial cells")
    
    temp_major_label <- "Endothelial"
    temp_subset@meta.data$celltype_major <- temp_major_label
    temp_subset@meta.data$celltype_minor <- temp_subset@active.ident
    temp_subset@meta.data$celltype_subset <- temp_subset@active.ident
    
  }
  
  print(table(temp_subset@meta.data$celltype_major))
  print(table(temp_subset@meta.data$celltype_minor))
  print(table(temp_subset@meta.data$celltype_subset))
  
  n <- paste0("temp_subset_",cluster)
  assign(n,
         temp_subset)
}

# 09: SAVE MONOCLE ------------------------------------------------------------

dir.create("monocle_new_ordering/Rdata/")

for(cluster in c("CAFs", "PVL", "Endo")){
  HSMM <- get(paste0("HSMM_",cluster))
  
  saveRDS(HSMM,
          paste0("monocle_new_ordering/Rdata/monocle_object_",cluster,".Rdata"))
  
  # readRDS(paste0("monocle_new_ordering/Rdata/monocle_object_",cluster,".Rdata"))
  # n <- paste0("HSMM_",cluster)
  # assign(n, HSMM)
}








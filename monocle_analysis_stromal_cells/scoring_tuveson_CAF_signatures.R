# STROMAL CELLS - SCORE TUVESON CAF SIGNATURES
# SUNNY Z WU
#
#
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_seurat_dev
# R
#
# 01: SETUP -------------------------------------------------------------------

library(GSEABase)
library(AUCell)
library(reshape2)
library(NMF)
library(Matrix)
library(cowplot)

# directories
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/stromal_AUCEll/")
dir.create("02_AUCELL_SIGNATURES")


# 02: READ SIGNATURES -----------------------------------------------------

# # read genesets
temp_genesets <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/05_figures/05_STROMAL_FIGURES_v2/gene_signatures.csv"

temp_xcell_genesets <-
  read.csv(temp_genesets)

temp_xcell_genesets <- temp_xcell_genesets[,2:6]

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

# 03: LOAD ANNOTATED OBJECTS ---------------------------------------------

# load all mesenchymal cells 
seurat_10X_integrated <- readRDS("../05_STROMAL_FIGURES_v3/Rdata/Stromal_ANNOTATED_object.Rdata")


# 04: RUN AUCELL ------------------------------

temp_exprMatrix <- GetAssayData(object = seurat_10X_integrated,
                                assay = "RNA",
                                slot = "data")

temp_exprMatrix <-
  Matrix(temp_exprMatrix,
         sparse = T)

# AUCELL
dim(temp_exprMatrix)

print("Buildrankings")
temp_cells_rankings <-
  AUCell_buildRankings(temp_exprMatrix,
                       nCores = 1,
                       plotStats = F)

# subset gene sets
temp_subsetgeneSets <-
  subsetGeneSets(temp_gene_set_collection,
                 rownames(temp_exprMatrix))

# calculate area under the curve
temp_cells_AUC <-
  AUCell_calcAUC(geneSets = temp_subsetgeneSets,
                 rankings = temp_cells_rankings,
                 aucMaxRank = ceiling(0.05 * nrow(temp_cells_rankings)),
                 nCores = 1,
                 verbose = T)

#transpose matrix for seurat metadata assignment
temp_cells_AUC_matrix <-
  t(as.data.frame(getAUC(temp_cells_AUC)))

temp_cells_AUC_matrix_sorted <-
  temp_cells_AUC_matrix[rownames(seurat_10X_integrated@meta.data),,drop=FALSE]

temp_cells_AUC_matrix_sorted <-
  as.data.frame.matrix(temp_cells_AUC_matrix_sorted)

print(all.equal(rownames(temp_cells_AUC_matrix_sorted),
                rownames(seurat_10X_integrated@meta.data)))

temp_first_geneset <-
  (ncol(seurat_10X_integrated@meta.data) + 1)

seurat_10X_integrated_AUCell <- AddMetaData(seurat_10X_integrated,
                                              metadata = temp_cells_AUC_matrix_sorted)


temp_last_geneset <-
  (ncol(seurat_10X_integrated_AUCell@meta.data))



# 05: FEATUEREPLOTS AND VLNPLOTS ------------------------------------------------------------

temp_gene_set_names <-
  c("iCAF_human_PDAC_Elyada_et_al_2019",
    "myCAF_human_PDAC_Elyada_et_al_2019",
    "iCAF_mouse_PDAC_Elyada_et_al_2019",
    "myCAF_mouse_PDAC_Elyada_et_al_2019",
    "apCAF_mouse_PDAC_Elyada_et_al_2019")    

colnames(seurat_10X_integrated_AUCell@meta.data)[87:91] <- 
  temp_gene_set_names

# PLOT
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
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 16,
      height = 12,
      useDingbats=F
    )
  }

temp_pdf_function(paste0("02_AUCELL_SIGNATURES/01_featureplot.pdf"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated_AUCell,
    features = temp_gene_set_names,
    order = T,
    pt.size = 0.01,
    reduction = temp_reduction)
  print(temp_featureplot)
  dev.off()

# vlnplot
temp_png_function <- 
  function(x) {
    png(
      file = (x), 
      width = 6, 
      height = 18, 
      res = 300, 
      units = 'in'
    )
  }
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 6,
      height = 18,
      useDingbats=F
    )
  }

temp_pdf_function(paste0("02_AUCELL_SIGNATURES/02_vlnplot.pdf"))
temp_vlnplot <- VlnPlot(object = seurat_10X_integrated_AUCell,
                        features = temp_gene_set_names,
                        pt.size = 0.01,
                        group.by = "celltype_subset",
                        ncol = 1)
print(temp_vlnplot)
dev.off()

temp_pdf_function(paste0("02_AUCELL_SIGNATURES/03_vlnplot_nodots.pdf"))
temp_vlnplot <- VlnPlot(object = seurat_10X_integrated_AUCell,
                        features = temp_gene_set_names,
                        pt.size = 0,
                        group.by = "celltype_subset", 
                        ncol = 1)
print(temp_vlnplot)
dev.off()


# 06: HEATMAPS ----------------------------------------------------------------

temp_df <- seurat_10X_integrated_AUCell@meta.data[,temp_gene_set_names]
temp_df$cluster <- seurat_10X_integrated_AUCell@meta.data$celltype_subset
temp_df_m <- melt(temp_df)
temp_df_m_agg <-
  aggregate(.~cluster+variable,
            temp_df_m,
            mean)

temp_df_m_agg_dcast <-
  dcast(data = temp_df_m_agg,
        formula = variable~cluster,
        fun.aggregate = sum,
        value.var = "value")

rownames(temp_df_m_agg_dcast) <- 
  temp_df_m_agg_dcast$variable

temp_df_m_agg_dcast <- 
  temp_df_m_agg_dcast[, ! colnames(temp_df_m_agg_dcast) %in% "variable"]

temp_df_m_agg_dcast <- as.matrix(temp_df_m_agg_dcast)

library("viridis")  
hmcol <- inferno(24)

rownames(temp_df_m_agg_dcast) <- c("iCAF_human_PDAC_Elyada_et_al_2019",
                                   "myCAF_human_PDAC_Elyada_et_al_2019",
                                   "iCAF_mouse_PDAC_Elyada_et_al_2019",
                                   "myCAF_mouse_PDAC_Elyada_et_al_2019",
                                   "apCAF_mouse_PDAC_Elyada_et_al_2019")

pheatmap(temp_df_m_agg_dcast, filename = paste0("02_AUCELL_SIGNATURES/04_heatmap.pdf"),
         color = rev(hmcol),
         cluster_cols = T, 
         cluster_rows = T, 
         scale = "row", 
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         fontsize_row = 10,
         show_rownames = T,
         show_colnames = T, 
         fontsize_col = 10,
         cellheight= 25, 
         cellwidth = 15,
         gaps_col = NULL, 
         angle_col = 45,
         treeheight_col = 50,
         legend = T
)

# SAVE RDS ----------------------------------------------------------------

saveRDS(seurat_10X_integrated_AUCell,
        "../05_STROMAL_FIGURES_v3/Rdata/Stromal_ANNOTATED_object_AUCell_Elyada_et_al.Rdata")






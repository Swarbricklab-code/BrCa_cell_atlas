# SPATIAL GENE MODULE ANALYSIS
# Sunny Z. Wu
# 
#
# Run in screen 
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_spatial
# R
# 
#
# 01: SETUP---------------------------------------

# setup
library(zeallot)
library(STutility)
library(ggplot2)
library(magrittr)
library(Seurat)
library(magick)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# For scoring metagenes
library(AUCell)
library(GSEABase)
library(Matrix)

# for stats
library(rstatix)


# DIRECTORY
dir.create("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/brca_gene_modules/")
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/brca_gene_modules/")
dir.create("output")


# 02: LOAD DATA ---------------------------------------------------------------

se.list <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/RDATA_Visium_brca_objects_stereoscope.Rdata")

# load gene module gene sets
temp_metagenes <- read.delim("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/11_visium_brca_metagenes/gene_clusters.txt")

# load stereoscope cluster IDs per annotation tier
temp_sample_id <- c(2,5,7:10)
temp_celltypes_CTP_minor <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_CTP_minor.Rdata")
temp_celltypes_minor <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_minor.Rdata")
temp_celltypes_subset <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_subset.Rdata")
temp_celltypes_major  <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019//14_spatial_CAF_Tcell_signalling/RDATA/celltypes_major.Rdata")


# 03: MAKE GENE SET COLLECTION -----------------------------------------

# make genesetcollection for AUCell
temp_GeneSetCollection_list <- NULL
for(cluster in unique(temp_metagenes$cluster)){
  print(cluster)
  temp_metagenes_cluster <- temp_metagenes[temp_metagenes$cluster == cluster,,drop=F]
  temp_metagenes_cluster <- as.vector(temp_metagenes_cluster$gene)
  # qc length of genesets, should be 200
  print(length(temp_metagenes_cluster))
  
  
  temp_set_name <- paste0("metagene_cluster_",cluster)
  
  # make gene set
  temp_set <- GeneSet(unique(as.vector(temp_metagenes_cluster)),
                      setName = temp_set_name)
  
  temp_GeneSetCollection_list <- append(temp_GeneSetCollection_list,
                                        temp_set)
  
  # n <- paste0("temp_metagenes_", cluster)
  # assign(n, temp_metagenes_cluster)
}

temp_gene_set_collection <- 
  GeneSetCollection(temp_GeneSetCollection_list)

# 04: EXPORT MATRICES FOR AUCELL ----------------------------------------------

for(sample in c(1:length(se.list))){
  
  temp_seurat_object <- se.list[[sample]]
  temp_matrix <- GetAssayData(temp_seurat_object,
                              assay="SCT")
  temp_matrix <- 
    Matrix(temp_matrix,
           sparse = T)
  
  n <- paste0("temp_matrix_",sample)
  assign(n, temp_matrix)
}

# 05: RUN AUCELL --------------------------------------------------------------

for(sample in c(1:length(se.list))){
  print(sample)
  temp_matrix <- get(paste0("temp_matrix_",sample))
  
  temp_cells_rankings <- 
    AUCell_buildRankings(temp_matrix, 
                         nCores = 1, 
                         plotStats = F)
  
  # subset gene sets
  temp_subsetgeneSets <- 
    subsetGeneSets(temp_gene_set_collection, 
                   rownames(temp_matrix)) 
  
  # calculate area under the curve
  temp_cells_AUC <- 
    AUCell_calcAUC(geneSets = temp_subsetgeneSets, 
                   rankings = temp_cells_rankings, 
                   aucMaxRank = ceiling(0.05 * nrow(temp_cells_rankings)))
  
  #transpose matrix for seurat metadata assignment
  temp_cells_AUC_matrix <- 
    t(as.data.frame(getAUC(temp_cells_AUC)))
  
  n <- paste0("temp_cells_AUC_matrix_",sample)
  assign(n, temp_cells_AUC_matrix)
  
}


# 06: APPEND TO SEURAT OBJECT --------------------------------------------------

for(sample in c(1:length(se.list))){
  print(sample)
  temp_seurat_object <- se.list[[sample]]
  temp_cells_AUC_matrix <- get(paste0("temp_cells_AUC_matrix_",sample))
  
  temp_cells_AUC_matrix_sorted <-
    temp_cells_AUC_matrix[rownames(temp_seurat_object@meta.data),,drop=FALSE]
  
  temp_cells_AUC_matrix_sorted <- 
    as.data.frame.matrix(temp_cells_AUC_matrix_sorted)
  
  if(all.equal(rownames(temp_cells_AUC_matrix_sorted),
               rownames(temp_seurat_object@meta.data))){
    print("rownames aligned")
    
    temp_seurat_object <- AddMetaData(temp_seurat_object, 
                                      metadata = temp_cells_AUC_matrix_sorted)
    
    colnames(temp_cells_AUC_matrix_sorted) <- gsub("metagene_cluster",
                                                   "metagenecluster",
                                                   colnames(temp_cells_AUC_matrix_sorted))
    
    # Set sc.data as cell.embeddings
    cell.embeddings <- temp_cells_AUC_matrix_sorted
    
    reduction.data <- CreateDimReducObject (
      embeddings = as.matrix(cell.embeddings),
      loadings = matrix(),
      assay = "SCT",
      key = "metagenecluster_"
    )
    
    # Store reduction object in reductions slot
    temp_seurat_object[["metageneclusters"]] <- reduction.data
    
    se.list[[sample]] <- temp_seurat_object
    
    # n <- paste0("temp_seurat_object_",sample)
    # assign(n, temp_seurat_object)
  } else {
    print("rownames not aligned")
  }
  
}







# SAVE RDS ----------------------------------------------------------------

saveRDS(se.list, "RDATA_Visium_brca_objects_stereoscope_AUCELL_gene_modules.Rdata")

# 07: FILTER CANCER ONLY SPOTS ---------------------------

# epithelial cancer histogram
temp_df_m_combined_cancer <- NULL
for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list[[sample]]@meta.data$rep))
  
  for(tier in c("minor")){
    
    if(tier == "minor" && sample %in% 7:10){
      temp_vector <- get(paste0("temp_celltypes_CTP_",tier))
    } else {
      temp_vector <- get(paste0("temp_celltypes_",tier))
    }
    
    if(tier == "minor"){
      temp_df <- se.list[[sample]]@reductions$stereoscope_tnbcminor@cell.embeddings
      temp_cellwidth <- 0.25
      temp_cellheight <- 10
      temp_pdf_function <-
        function(x) {
          pdf(file = (x),
              width = 12, 
              height = 6)
        }
    }
    colnames(temp_df) <- temp_vector
    temp_df <- as.data.frame(temp_df)
    temp_df$seurat_clusters <- as.vector(se.list[[sample]]@meta.data$seurat_clusters)
    temp_df$Classification <- as.vector(se.list[[sample]]@meta.data$Classification)
    
    temp_df$barcode <- rownames(temp_df)
    temp_df_m <- reshape2::melt(temp_df)
    temp_df_m$seurat_clusters <- factor(temp_df_m$seurat_clusters,
                                        levels=stringr::str_sort(unique(temp_df_m$seurat_clusters), numeric = T))
    temp_df_m <- temp_df_m[temp_df_m$variable == "Epithelial_cancer",,drop=F]
    # rank
    temp_df_m <- temp_df_m[order(-temp_df_m$value),,drop=T]
    temp_df_m$barcode_number <- c(1:length(temp_df_m$barcode))
    temp_df_m$sample <- temp_sampleID
    temp_df_m$patient <- unique(se.list[[sample]]@meta.data$patientid)
    temp_df_m$rounded <- round(temp_df_m$value, digits = 2)
    
    print(paste0(temp_sampleID," 0 epithelial = ", 
                 (nrow(temp_df_m[temp_df_m$rounded == 0,,drop=F])/nrow(temp_df_m))*100,
                 "%"
    )
    )
    
    temp_df_m_combined_cancer <- rbind(temp_df_m_combined_cancer,
                                       temp_df_m)
    
    
  }
  
  temp_df_m_plot <- temp_df_m
  temp_df_m_plot <- temp_df_m_plot[!temp_df_m_plot$Classification == "Artefact",,drop=F]
  temp_df_m_plot <- temp_df_m_plot[!temp_df_m_plot$Classification == "NA",,drop=F]
  temp_df_m_plot <- temp_df_m_plot[!temp_df_m_plot$Classification == "",,drop=F]
  
  temp_ggplot <- ggplot(temp_df_m_plot,
                        aes(x = barcode_number,
                            y = value,
                            fill = Classification)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = c(0,0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(breaks = seq(0,5000,500)) +
    # theme(axis.text.x = element_blank()) +
    ylab("Epithelial Cancer Deconvolution") +
    xlab("Spots (ranked)") +
    facet_wrap(. ~ Classification,
               scales="free_x",
               ncol = 2) + # free x removes white space
    geom_hline(yintercept=0.1, 
               color = "red", 
               size=0.5)
  
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
  
  temp_png_function(paste0("output/02_epithelialcancer_decon_",temp_sampleID,"_facet_by_path.png"))
  print(temp_ggplot)
  dev.off()
}

# combined plot all samples
temp_df_m_combined_cancer <- temp_df_m_combined_cancer[!temp_df_m_combined_cancer$Classification %in% c("Necrosis", "Uncertain", "Artefact","NA",""),,drop=F]

temp_ggplot <- ggplot(temp_df_m_combined_cancer,
                      aes(x = barcode_number,
                          y = value,
                          fill = Classification)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = c(0,0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = seq(0,5000,1000)) +
  # theme(axis.text.x = element_blank()) +
  ylab("Epithelial Cancer Deconvolution") +
  xlab("Spots (ranked)") +
  facet_wrap(. ~ patient,
             scales="free_x",
             ncol = 3) + # free x removes white space
  geom_hline(yintercept=0.1, 
             color = "red", 
             size=0.5)

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 10,
      height = 6,
      res = 300,
      units = 'in'
    )
  }

temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 10,
      height = 6,
      useDingbats = F
    )
  }

temp_pdf_function(paste0("output/03_epithelialcancer_histogram.pdf"))
print(temp_ggplot)
dev.off()

# filter epithelial cancer by 10%
temp_df_m_combined_cancer <- NULL
for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list[[sample]]@meta.data$rep))
  
  for(tier in c("minor")){
    
    if(tier == "minor" && sample %in% 7:10){
      temp_vector <- get(paste0("temp_celltypes_CTP_",tier))
    } else {
      temp_vector <- get(paste0("temp_celltypes_",tier))
    }
    
    if(tier == "minor"){
      temp_df <- se.list[[sample]]@reductions$stereoscope_tnbcminor@cell.embeddings
      temp_cellwidth <- 0.25
      temp_cellheight <- 10
      temp_pdf_function <-
        function(x) {
          pdf(file = (x),
              width = 12, 
              height = 6)
        }
    }
    colnames(temp_df) <- temp_vector
    temp_df <- as.data.frame(temp_df)
    temp_df$seurat_clusters <- as.vector(se.list[[sample]]@meta.data$seurat_clusters)
    temp_df$Classification <- as.vector(se.list[[sample]]@meta.data$Classification)
    
    temp_df$barcode <- rownames(temp_df)
    temp_df_m <- reshape2::melt(temp_df)
    temp_df_m$seurat_clusters <- factor(temp_df_m$seurat_clusters,
                                        levels=stringr::str_sort(unique(temp_df_m$seurat_clusters), numeric = T))
    temp_df_m <- temp_df_m[temp_df_m$variable == "Epithelial_cancer",,drop=F]
    # rank
    temp_df_m <- temp_df_m[order(-temp_df_m$value),,drop=T]
    temp_df_m$barcode_number <- c(1:length(temp_df_m$barcode))
    temp_df_m$sample <- temp_sampleID
    temp_df_m$rounded <- round(temp_df_m$value, digits = 2)
    
    print(paste0(temp_sampleID," 0 epithelial = ", 
                 (nrow(temp_df_m[temp_df_m$rounded == 0,,drop=F])/nrow(temp_df_m))*100,
                 "%"
    )
    )
    
    temp_df_m_combined_cancer <- rbind(temp_df_m_combined_cancer,
                                       temp_df_m)
    
    
  }
  
  temp_df_m_plot <- temp_df_m
  temp_df_m_plot <- temp_df_m_plot[!temp_df_m_plot$Classification == "Artefact",,drop=F]
  temp_df_m_plot <- temp_df_m_plot[!temp_df_m_plot$Classification == "NA",,drop=F]
  temp_df_m_plot <- temp_df_m_plot[!temp_df_m_plot$Classification == "",,drop=F]
  
  temp_ggplot <- ggplot(temp_df_m_plot,
                        aes(x = barcode_number,
                            y = value,
                            fill = Classification)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = c(0,0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(breaks = seq(0,5000,500)) +
    # theme(axis.text.x = element_blank()) +
    ylab("Epithelial Cancer Deconvolution") +
    xlab("Spots (ranked)") +
    facet_wrap(. ~ Classification,
               scales="free_x",
               ncol = 2) # free x removes white space
  
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
  
  temp_png_function(paste0("output/02_epithelialcancer_decon_",temp_sampleID,"_facet_by_path.png"))
  print(temp_ggplot)
  dev.off()
}


temp_cutoff <- 0.1
# filter epithelial cancer
temp_df_m_combined_filtered <- NULL
for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list[[sample]]@meta.data$rep))
  temp_df_m_combined_subset <- temp_df_m_combined_cancer[temp_df_m_combined_cancer$sample == temp_sampleID,,drop=F]
  temp_num <- nrow(temp_df_m_combined_subset)
  temp_df_m_combined_subset <- temp_df_m_combined_subset[temp_df_m_combined_subset$rounded > temp_cutoff,,drop=F]
  print(paste0(temp_sampleID,
               " filtered ",
               (temp_num-nrow(temp_df_m_combined_subset))))
  temp_df_m_combined_filtered <- rbind(temp_df_m_combined_filtered, temp_df_m_combined_subset)
}


# filter all normal regions
# sort(unique(temp_df_m_combined_filtered$Classification))
temp_df_m_combined_filtered <- 
  temp_df_m_combined_filtered[!temp_df_m_combined_filtered$Classification %in% 
                                c("Normal Gland", "Normal + stroma + lymphocytes","Necrosis", 
                                  "Stroma + adipose tissue",
                                  "Artefact", "", "NA"),,drop=F]

# filter objects
se.list_filtered <- se.list
for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list_filtered[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list_filtered[[sample]]@meta.data$rep))
  
  temp_df_m_combined_subset <- temp_df_m_combined_filtered[temp_df_m_combined_filtered$sample == temp_sampleID,,drop=F]
  
  temp_seurat_object <- se.list_filtered[[sample]]
  
  temp_seurat_object <- subset(temp_seurat_object,
                               cells=as.vector(temp_df_m_combined_subset$barcode))
  
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 12,
        height = 10,
        res = 300,
        units = 'in'
      )
    }
  
  
  temp_png_function(paste0("output/04_DimOverlay_filtered_by_epithelialdecon_",temp_sampleID,".png"))
  print(DimOverlay(temp_seurat_object,
                   reduction = "metageneclusters",
                   dims = 1:7,
                   pt.size = 0.25,
                   pt.alpha = 0.5,
                   ncols.samples = 1,
                   sample.label = F,
                   ncols.dims = 3,
                   cols = c("navyblue", "cyan", "yellow", "red", "dark red"),
                   center.zero = FALSE,
                   value.scale = "samplewise")
  )
  dev.off()
  
  se.list_filtered[[sample]] <- temp_seurat_object
  rm(temp_seurat_object)
}


# SAVE RDS ----------------------------------------------------------------

saveRDS(se.list_filtered, "RDATA_Visium_brca_objects_stereoscope_AUCELL_gene_modules_cancer_filtered_10%.Rdata")

# 08: STATS - MAKE DF ALL 6 SAMPLES --------------------

# combined df
temp_df_combined <- NULL
for(sample in c(2,5,7:10)){
  print(sample)
  temp_sampleID <- paste0(unique(se.list_filtered[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list_filtered[[sample]]@meta.data$rep))
  
  temp_colnames <- grep("meta",colnames(se.list_filtered[[sample]]@meta.data),value=T)
  temp_df_metagenes <- se.list_filtered[[sample]]@meta.data[, colnames(se.list_filtered[[sample]]@meta.data) %in% temp_colnames,drop=F]
  temp_df_metagenes_m <- reshape2::melt(temp_df_metagenes)
  
  temp_df_metagenes_m$sample <- temp_sampleID
  temp_df_metagenes_m$patient <- unique(se.list_filtered[[sample]]@meta.data$patientid)
  temp_df_metagenes_m$subtype <- unique(se.list_filtered[[sample]]@meta.data$subtype)
  
  temp_df_combined <- rbind(temp_df_combined, temp_df_metagenes_m)
  
  rm(temp_df_metagenes_m)
  
}

# plot
temp_df_combined$subtype <- factor(temp_df_combined$subtype,
                                   levels=c("TNBC", "ER"))
temp_df_combined$patient <- factor(temp_df_combined$patient,
                                   levels=c("CID44971", "CID4465", 
                                            "1160920F",
                                            "1142243F",
                                            "CID4535", "CID4290"))

temp_df_combined$variable <- gsub("metagene_cluster_","GM", temp_df_combined$variable)


# 09: STATS - T-TEST MULTIPLE TESTING -------------------------------------

# stats pariwise by sample
stat.test <- temp_df_combined %>%
  group_by(variable) %>%
  t_test(value ~ patient, p.adjust.method = "BH")
stat.test <- as.data.frame(stat.test)
stat.test$pvalmethod <- "BH"
write.csv(stat.test, "t_test_by_sample.csv")

# stats pariwise by subtype
# no multiple testing correction required for one comparison of TNBC vs ER
stat.test <- temp_df_combined %>%
  group_by(variable) %>%
  t_test(value ~ subtype, p.adjust.method = "none")
stat.test <- as.data.frame(stat.test)

# adjusted p-values
stat.test$p_val_adjust =
  p.adjust(stat.test$p,
           method = "BH")

stat.test$pvalmethod <- "BH"
write.csv(stat.test, "t_test_by_clinical_subtype.csv")

# 10: STATS - PLOTS ---------------------------------------------------------------

# comparisons
my_comparisons <- list(c("1142243F", "CID4535"),
                       c("1142243F", "CID4290"),
                       c("1160920F", "CID4535"),
                       c("1160920F", "CID4290"),
                       c("CID44971", "CID4535"),
                       c("CID44971", "CID4290"),
                       c("CID4465", "CID4535"),
                       c("CID4465", "CID4290")
)

# boxplot comparemeans
temp_gene_boxplot <- ggboxplot(temp_df_combined, 
                               x = "patient",
                               y= "value",
                               fill = "subtype",
                               # legend = "none",
                               facet.by = "variable",
                               repel = F,
                               nrow=1,
                               size=0.05, 
                               scales ="free_x",
                               outlier.size = 0.001
                               # ylim=c(0,100000)
) + stat_compare_means(method = "t.test",
                       paired = F,
                       label = "p.signif", 
                       hide.ns = T,tip.length = 0.01,
                       comparisons = my_comparisons,
                       bracket.size = 0.001,
                       # label.y = 0.45,
                       # symnum.args = symnum.args
) +  # Pairwise comparison against all
  xlab("Sample ID") + ylab("Signature Score (AUCell)") + theme(axis.text.x = element_text(angle=270, 
                                                                                          size=8),
                                                               axis.title.y = element_text(size=15)) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_y_continuous(position = "left") + 
  guides(fill=guide_legend(title="Clinical Subtype"))

pdf(paste0("output/05_boxplot_facet_by_celltype_all_samples_for_figure_cancer_filtered_10percent.pdf"), 
    width = 10, 
    height = 5)
print(temp_gene_boxplot)
dev.off()

# comparemeans grouped by subtype
my_comparisons <- list(c("TNBC", "ER")
)

# boxplot comparemeans
library(ggpubr)
temp_gene_boxplot <- ggboxplot(temp_df_combined, 
                               x = "subtype",
                               y= "value",
                               fill = "subtype",
                               # legend = "none",
                               facet.by = "variable",
                               repel = F,
                               nrow=1,
                               size=0.05, 
                               scales ="free_x",
                               outlier.size = 0.001
                               # ylim=c(0,100000)
) + stat_compare_means(method = "t.test",
                       paired = F,
                       label = "p.signif", 
                       hide.ns = T,
                       comparisons = my_comparisons,
                       tip.length = 0.01,
                       bracket.size = 0.001,
                       # label.y = temp_label_y,
                       # symnum.args = symnum.args
) +  # Pairwise comparison against all
  xlab("Subtype") + ylab("Signature Score (AUCell)") + theme(axis.text.x = element_text(angle=270, 
                                                                                        size=8),
                                                             axis.title.y = element_text(size=15)) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_y_continuous(position = "left") + 
  guides(fill=guide_legend(title="Clinical Subtype"))

pdf(paste0("output/06_boxplot_facet_by_celltype_all_samples_for_figure_cancer_filtered_10percent_by_subtype.pdf"), 
    width = 10, 
    height = 5)
print(temp_gene_boxplot)
dev.off()

# 11: METAGENE CLUSTERS CORRELATION PER SAMPLE ----------------------------------------

# correlation p-value function
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

temp_pval_combined <- NULL
for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list_filtered[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list_filtered[[sample]]@meta.data$rep))
  
  temp_colnames <- grep("meta",colnames(se.list_filtered[[sample]]@meta.data),value=T)
  temp_df_metagenes <- se.list_filtered[[sample]]@meta.data[, colnames(se.list_filtered[[sample]]@meta.data) %in% temp_colnames,drop=F]
  
  corMat <- cor(temp_df_metagenes)
  corMat_test <- cor.test.p(temp_df_metagenes)
  
  temp_breaks <-  seq(-1, 1, length.out = 101)
  temp_legendbreak <- seq(-1,1,0.2)
  
  rownames(corMat) <- gsub("metagene_cluster_","GM",rownames(corMat))
  colnames(corMat) <- gsub("metagene_cluster_","GM",colnames(corMat))
  
  rownames(corMat_test) <- gsub("metagene_cluster_","GM",rownames(corMat_test))
  colnames(corMat_test) <- gsub("metagene_cluster_","GM",colnames(corMat_test))
  
  
  pheatmap::pheatmap(corMat,
                     filename = paste0("output/07_corMat_sample_filtered_",temp_sampleID,".pdf"),
                     main = temp_sampleID,
                     cluster_rows = T,cellwidth = 18,cellheight = 18,
                     cluster_cols = T,
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                     breaks = temp_breaks,treeheight_col = 25, treeheight_row = 25,
                     legend_breaks = temp_legendbreak,
                     na_col = "grey")
  
  # filter original matrix for significance
  # all non-significant value indicated by grey
  for(GM1 in unique(rownames(corMat))){
    for(GM2 in unique(colnames(corMat))){
      temp_pval <- corMat_test[GM1,GM2]
      if(temp_pval > 0.05){
        print(paste0("correlation ", GM1, GM2))
        print("not significant")
        corMat[GM1,GM2] <- NA
      }
    }
  }
  
  pheatmap::pheatmap(corMat,
                     filename = paste0("output/08_corMat_sample_filtered_ns_",temp_sampleID,".pdf"),
                     main = temp_sampleID,
                     cluster_rows = T,cellwidth = 18,cellheight = 18,
                     cluster_cols = T,
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                     breaks = temp_breaks,treeheight_col = 25, treeheight_row = 25,
                     legend_breaks = temp_legendbreak,
                     na_col = "grey")
  
  # remove duplicates prior to p-val adjustment
  temp_GM <- c(1)
  corMat_test_m <- NULL
  for(GM in unique(colnames(corMat_test))){
    
    corMat_test_subset <- corMat_test[,GM,drop=F]
    corMat_test_subset <- reshape2::melt(corMat_test_subset)
    
    for(num in temp_GM){
      corMat_test_subset <- corMat_test_subset[! corMat_test_subset$Var1 == paste0("GM",num),,drop=F]
    }
    
    corMat_test_m <- rbind(corMat_test_m, corMat_test_subset)
    temp_GM <- c(temp_GM,length(temp_GM)+1)
    }
  
  colnames(corMat_test_m) <- c("GM1", "GM2", "p_val")
  corMat_test_m$GM_ID1 <- paste0(corMat_test_m$GM2,"_",corMat_test_m$GM1)
  corMat_test_m$sample <- temp_sampleID
  
  # adjusted p-values
  corMat_test_m$p_val_adjust =
    p.adjust(corMat_test_m$p_val,
             method = "BH")
  
  temp_pval_combined <- rbind(temp_pval_combined, corMat_test_m)
}
rownames(temp_pval_combined) <- c(1: nrow(temp_pval_combined))

temp_pval_combined$significant <- NA
for(row in c(1:nrow(temp_pval_combined))){
  if(temp_pval_combined[row,"p_val"] < 0.05 ){
    temp_pval_combined[row,"significant"] <- "Y"
  } else {
    temp_pval_combined[row,"significant"] <- "N"
  }
}

write.csv(temp_pval_combined, "pearson_corr_gene_modules.csv")


# 12: METAGENE CLUSTERS: PLOT MG3vsMG4 ACROSS ALL SAMPLES ----------------------------------------------

for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list_filtered[[sample]]@meta.data$patientid))
  
  temp_colnames <- grep("meta",colnames(se.list_filtered[[sample]]@meta.data),value=T)
  temp_df_metagenes <- se.list_filtered[[sample]]@meta.data[, colnames(se.list_filtered[[sample]]@meta.data) %in% temp_colnames,drop=F]
  
  corMat <- cor(temp_df_metagenes)
  temp_df_m <- reshape2::melt(corMat)
  
  rownames(temp_df_m) <- paste0(temp_df_m$Var1, "_", temp_df_m$Var2)
  
  if(sample == 2){
    temp_dataframe_combined <- temp_df_m
    temp_dataframe_combined <- temp_dataframe_combined[,2,drop=F]
  }
  
  temp_df_m <- temp_df_m[,3,drop=F]
  temp_df_m <- temp_df_m[order(rownames(temp_df_m)),,drop=F]
  
  colnames(temp_df_m) <-  temp_sampleID
  
  # temp_df_m$sample <- sample_id
  # temp_df_m$subtype <- subtype
  
  # temp_df_m_df <- data.frame(value = temp_df_m$value)
  # rownames(temp_df_m_df) <- paste0("sample_",sample,
  #                                  "_", temp_df_m$Var1,
  #                                  "_", temp_df_m$Var2)
  
  if(sample == 2){
    temp_cormat_combined <- temp_df_m
  } else (
    temp_cormat_combined <- cbind(temp_cormat_combined, temp_df_m)
  )
}

# plot by major cell lineage
temp_dataframe_col <- data.frame(row.names = colnames(temp_cormat_combined),
                                 sample = colnames(temp_cormat_combined),
                                 subtype = c("TNBC", "TNBC", "TNBC", "TNBC",
                                             "ER+", "ER+"))

# plot
temp_clusternames <- grep("meta",colnames(se.list_filtered[[sample]]@meta.data),value=T)

for(subset_lineage in temp_clusternames){
  print(subset_lineage)
  
  temp_cormat_combined_subset <- temp_cormat_combined[rownames(temp_cormat_combined) %in% grep(paste0("^",subset_lineage), rownames(temp_cormat_combined),value=T),,drop=F]
  
  temp_cormat_combined_subset <- temp_cormat_combined_subset[!rownames(temp_cormat_combined_subset) %in% paste0(subset_lineage,"_",subset_lineage),,drop=F]
  
  pheatmap::pheatmap(temp_cormat_combined_subset,
                     filename = paste0("output/09_pheatmap_",subset_lineage,".pdf"),
                     main = subset_lineage,
                     cluster_rows = T,
                     cluster_cols = T,
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                     breaks = seq(-1, 1, length.out = 101),
                     legend_breaks = seq(-1,1,0.2),
                     annotation_col = temp_dataframe_col,
                     fontsize_row = 6, 
                     na_col = "black")
}

# just MG3 vs MG4
for(subset_lineage in "metagene_cluster_3"){
  print(subset_lineage)
  
  temp_cormat_combined_subset <- temp_cormat_combined[rownames(temp_cormat_combined) %in% grep(paste0("^",subset_lineage), rownames(temp_cormat_combined),value=T),,drop=F]
  
  temp_cormat_combined_subset <- temp_cormat_combined_subset[!rownames(temp_cormat_combined_subset) %in% paste0(subset_lineage,"_",subset_lineage),,drop=F]
  
  temp_cormat_combined_subset <- temp_cormat_combined_subset[rownames(temp_cormat_combined_subset) %in% paste0(subset_lineage,"_","metagene_cluster_4"),,drop=F]
  
  temp_dataframe_col_subset <- temp_dataframe_col[,2,drop=F]
  
  # rownames(temp_cormat_combined_subset) <- gsub("metagene_cluster_","MG",rownames(temp_cormat_combined_subset))
  rownames(temp_cormat_combined_subset) <- "MG3_vs_MG4"
  
  pdf(file = paste0("output/10_MG3_vs_MG4pheatmap_",subset_lineage,".pdf"), height = 3, width = 7
  )
  pheatmap::pheatmap(temp_cormat_combined_subset,
                     main = subset_lineage,
                     cluster_rows = F,
                     cluster_cols = T,
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                     breaks = seq(-1, 1, length.out = 101),
                     legend_breaks = seq(-1,1,0.5),
                     annotation_col = temp_dataframe_col_subset,
                     annotation_colors = list(
                       subtype = c("ER+" = "blue", 
                                   TNBC = "red"
                       )),
                     cellwidth = 50, 
                     cellheight = 25,
                     show_rownames = T,
                     show_colnames = T)
  dev.off()
  
}


# 13: METAGENE CLUSTERS: PLOT ALL GMs ACROSS ALL SAMPLES ----------------------------------------------

for(sample in c(2,5,7:10)){
  temp_sampleID <- paste0(unique(se.list_filtered[[sample]]@meta.data$patientid))
  
  temp_colnames <- grep("meta",colnames(se.list_filtered[[sample]]@meta.data),value=T)
  temp_df_metagenes <- se.list_filtered[[sample]]@meta.data[, colnames(se.list_filtered[[sample]]@meta.data) %in% temp_colnames,drop=F]
  
  corMat <- cor(temp_df_metagenes)
  temp_df_m <- reshape2::melt(corMat)
  
  rownames(temp_df_m) <- paste0(temp_df_m$Var1, "_", temp_df_m$Var2)
  
  if(!sample == 2){
    temp_df_m <- temp_df_m[,3,drop=F]
  }
  
  temp_df_m <- temp_df_m[order(rownames(temp_df_m)),,drop=F]
  
  if(sample == 2){
    colnames(temp_df_m) <-  c("Var1", "Var2", temp_sampleID)
  } else {
    colnames(temp_df_m) <-  temp_sampleID
  }
  
  # temp_df_m$sample <- sample_id
  # temp_df_m$subtype <- subtype
  
  # temp_df_m_df <- data.frame(value = temp_df_m$value)
  # rownames(temp_df_m_df) <- paste0("sample_",sample,
  #                                  "_", temp_df_m$Var1,
  #                                  "_", temp_df_m$Var2)
  
  if(sample == 2){
    temp_cormat_combined <- temp_df_m
  } else (
    temp_cormat_combined <- cbind(temp_cormat_combined, temp_df_m)
  )
}

temp_dataframe_col <- data.frame(row.names = colnames(temp_cormat_combined)[3:8],
                                 sample = colnames(temp_cormat_combined)[3:8],
                                 subtype = c("TNBC", "TNBC", "TNBC", "TNBC",
                                             "ER+", "ER+"))
temp_cormat_combined$reverse <- paste0(temp_cormat_combined$Var2, "_", temp_cormat_combined$Var1)
temp_clusternames <- grep("meta",colnames(se.list_filtered[[sample]]@meta.data),value=T)
temp_cormat_combined_trimmed <- NULL
temp_dataframe_row <- temp_cormat_combined[,colnames(temp_cormat_combined) %in% c("Var1", "Var2"),drop=F]
temp_dataframe_row$Var1 <- gsub("metagene_cluster_","GM",temp_dataframe_row$Var1)
temp_dataframe_row$Var2 <- gsub("metagene_cluster_","GM",temp_dataframe_row$Var2)

for(subset_lineage in temp_clusternames){
  print(subset_lineage)
  
  temp_cormat_combined_subset <- temp_cormat_combined[rownames(temp_cormat_combined) %in% grep(paste0("^",subset_lineage), rownames(temp_cormat_combined),value=T),,drop=F]
  
  temp_cormat_combined_subset <- temp_cormat_combined_subset[!rownames(temp_cormat_combined_subset) %in% paste0(subset_lineage,"_",subset_lineage),,drop=F]
  
  #remove reverse
  if(! subset_lineage %in% "metagene_cluster_1"){
    temp_cormat_combined_subset <- temp_cormat_combined_subset[!rownames(temp_cormat_combined_subset) %in% c(rownames(temp_cormat_combined_trimmed), temp_cormat_combined_trimmed$reverse),,drop=F]
  }
  
  temp_cormat_combined_trimmed <- rbind(temp_cormat_combined_trimmed, temp_cormat_combined_subset)
}


temp_dataframe_row <- temp_dataframe_row[rownames(temp_dataframe_row) %in% rownames(temp_cormat_combined_trimmed),,drop=F]
temp_cormat_combined_trimmed <- temp_cormat_combined_trimmed[,!colnames(temp_cormat_combined_trimmed) %in% c("Var1", "Var2","reverse"),drop=F]

temp_dataframe_col_subset <- temp_dataframe_col[,2,drop=F]

rownames(temp_cormat_combined_trimmed) <- gsub("metagene_cluster_","GM",rownames(temp_cormat_combined_trimmed))
rownames(temp_dataframe_row) <- gsub("metagene_cluster_","GM",rownames(temp_dataframe_row))
pheatmap::pheatmap(temp_cormat_combined_trimmed,
                   filename = paste0("output/11_ALL_MGs_combined.pdf"),
                   cluster_rows = T,
                   cluster_cols = T,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                   breaks = seq(-1, 1, length.out = 101),cellwidth = 15, cellheight = 15,
                   legend_breaks = seq(-1,1,0.5),
                   annotation_row = temp_dataframe_row,
                   annotation_col = temp_dataframe_col_subset,
                   annotation_colors = list(
                     subtype = c("ER+" = "blue", 
                                 TNBC = "red"
                     )),
                   # cellwidth = 10, 
                   # cellheight = 10,
                   show_rownames = T,
                   show_colnames = T)




# 14: METAGENE CLUSTERSFILTERED: TRIMMED MAPPING OF METAGENES SIMPLE GRID- ---------------------------

# PLOT SAMPLE CID44971
library(raster)
dir.create("output/12_mapping_over_simple_grids/")

for(sample in c(7)){
  temp_sampleID <- paste0(unique(se.list_filtered[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list_filtered[[sample]]@meta.data$rep))
  dir.create(paste0("output/12_mapping_over_simple_grids/",temp_sampleID))
  
  # set number of slices and coordinates
  if(sample == 7){
    temp_num_slices <- c(1:3)
    temp_coordinates_1 <- c(500, # x1
                            1080, # x2
                            1360, #y1
                            1680 # y2
    )
    temp_coordinates_2 <- c(1055, # x1
                            1465, # x2
                            830, #y1
                            1530 # y2
    )
    temp_coordinates_3 <- c(455, # x1
                            950, # x2
                            355, #y1
                            900 # y2
    )
    
  }
  if(sample == 8){
    temp_num_slices <- c(1:2)
    temp_coordinates_1 <- c(705, # x1
                            1250, # x2
                            1185, #y1
                            1650 # y2
    )
    temp_coordinates_2 <- c(895, # x1
                            1410, # x2
                            295, #y1
                            965 # y2
    )
  }
  if(sample == 9){
    temp_num_slices <- c(1:3)
    temp_coordinates_1 <- c(510, # x1
                            965, # x2
                            1310, #y1
                            1660 # y2
    )
    temp_coordinates_2 <- c(925, # x1
                            1355, # x2
                            925, #y1
                            1425 # y2
    )
    temp_coordinates_3 <- c(475, # x1
                            1065, # x2
                            325, #y1
                            1010 # y2
    )
    
  }
  if(sample == 10){
    temp_num_slices <- c(1)
    temp_coordinates_1 <- c(575, # x1
                            1585, # x2
                            285, #y1
                            1535 # y2
    )
  }
  
  # temp_seurat_object <- se.list[[sample]]
  # 
  # temp_seurat_object@reductions$metageneclusters@cell.embeddings[!rownames(temp_seurat_object@reductions$metageneclusters@cell.embeddings) %in% colnames(se.list_filtered[[sample]]),] <- NA
  
  # metagenes
  for(celltype in 1:7){
    
    temp_ggplot <- ST.DimPlot(se.list_filtered[[sample]], 
                              dims = celltype, 
                              reduction = paste0("metageneclusters"), 
                              cols = c("navyblue", "cyan", "yellow", "red", "dark red"),
                              center.zero = FALSE, 
                              value.scale = "all" 
    )
    
    temp_ggplot <- temp_ggplot+scale_fill_continuous(na.value="white")
    
    ggsave(plot=temp_ggplot,
           paste0("output/12_mapping_over_simple_grids/temp.tiff"),
           device = "tiff")
    
    pdf(file = paste0("output/12_mapping_over_simple_grids/",temp_sampleID,"/metagene_",celltype,"legend",".pdf"),
    )
    print(temp_ggplot)
    dev.off()
    
    temp_raster <- brick(paste0("output/12_mapping_over_simple_grids/temp.tiff"))
    
    for(slice in temp_num_slices){
      temp_coordinates <- get(paste0("temp_coordinates_",slice))
      temp_coordinates[2] <- temp_coordinates[2]+20
      
      temp_raster_clipped <- crop(temp_raster, extent(temp_coordinates)
      )
      # PNG function
      temp_png_function <-
        function(x) {
          png(
            file = (x),
            width = 8,
            height = 8,
            res = 600,
            units = 'in'
          )
        }
      
      temp_png_function(paste0("output/12_mapping_over_simple_grids/",temp_sampleID,"/metagene_",celltype,"_slice_", slice,".png"))
      plotRGB(temp_raster_clipped,
              axes = F,
              margins= F,
              maxpixels=(500000*4)
      )
      dev.off()
    }
    
  }
}






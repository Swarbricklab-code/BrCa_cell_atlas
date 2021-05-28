#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# nUMI
# nGene
# CNA values
# Normal cell calls
# SNP array CNVs
# CNV-associated cancer genes optional

project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_131019"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])
cancer_x_threshold_sd_multiplier <- as.numeric(args[3])
cancer_y_threshold_sd_multiplier <- as.numeric(args[4])
normal_x_threshold_sd_multiplier <- as.numeric(args[5])
normal_y_threshold_sd_multiplier <- as.numeric(args[6])
include_group_annotation <- as.logical(args[7])
adjust_normal_thresholds <- as.logical(args[8])
print(include_group_annotation)

#sample_name <- "CID3963"
#include_t_cells <- TRUE
#cancer_x_threshold_sd_multiplier <- 2
#cancer_y_threshold_sd_multiplier <- 1.5
#normal_x_threshold_sd_multiplier <- 1
#normal_y_threshold_sd_multiplier <- 1.25
#include_group_annotation <- FALSE
#adjust_normal_thresholds <- TRUE

# distinguish samples with 1 and > 1 clusters:
mono_samples <- c("CID3586", "CID3921", "CID3941", "CID3948", "CID4067", 
  "CID4290A", "CID4461", "CID4495", "CID4515", "CID4523", "CID4535", 
  "CID44041", "CID44991", "CID45171", "CID4513", "CID4398")
multi_samples <- c("CID3963", "CID4066", "CID4463", "CID4465", 
  "CID4471", "CID4530N", "CID44971")

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included?", as.character(include_t_cells)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(dplyr)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
integrated_dir <- paste0("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/",
	"Jun2019/04_reclustering_analysis/run06_v1.2.1/output/Epithelial/02_Rdata/")

if (adjust_normal_thresholds) {

  adjustment_df <- data.frame(
    sample_id = c("CID4523", "CID4530N"),
    x_int = c(0.1, 0.17),
    y_int = c(0.25, 0.3)
  )
  assign_unassigned <- c("CID4290A", "CID4463", "CID45171", "CID4515")
  all_cancer <- c("CID3586", "CID3963", "CID4067")
  above_y_axis_only <- c("CID4530N")

  if (include_t_cells) {
    in_dir <- paste0(results_dir, 
      "sup_figure/t_cells_included/normal_thresholds_adjusted/", 
      sample_name, "/")
  } else {
    in_dir <- paste0(results_dir, 
      "sup_figure/t_cells_excluded/normal_thresholds_adjusted/", 
      sample_name, "/")
  }

} else {

  if (include_t_cells) {
    in_dir <- paste0(results_dir, "sup_figure/t_cells_included/", sample_name, "/")
  } else {
    in_dir <- paste0(results_dir, "sup_figure/t_cells_excluded/", sample_name, "/")
  }

}

out_dir <- in_dir

Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))

if (adjust_normal_thresholds) {
  non_adjusted_dir <- gsub("normal_thresholds_adjusted/", "", Robject_dir)
  if (
    !file.exists(
      paste0(
        Robject_dir, 
        "3b.epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata"
      )
    )
  ) {
      for (
        fname in c(
          "1a.initial_epithelial_heatmap.Rdata",
          "1b.initial_epithelial_metadata.Rdata",
          "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata",
          "2b.epithelial_metadata_with_cell_type_and_QC.Rdata",
          "3a.epithelial_heatmap_with_cell_type_QC_CNA_and_correlation_values.Rdata",
          "3b.epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata"
        )
      )
      system(paste0("cp ", non_adjusted_dir, fname, " ", Robject_dir))
    }
}


################################################################################
### 0. Define functions ###
################################################################################

gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, 
  "prepare_infercnv_metadata.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
create_array_CNV_annotation <- dget(paste0(func_dir, 
  "create_array_CNV_annotation.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {

  # load InferCNV output:
  print("Loading InferCNV output files...")
  infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

  # create cluster metadata df:
  seurat_10X <- readRDS(paste0(seurat_dir, "03_seurat_object_processed.Rdata"))
  Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  metadata <- prepare_infercnv_metadata(seurat_10X, subset_data = F, 
    as.data.frame(t(infercnv_output)), for_infercnv=F)
  epithelial_metadata <- metadata$metadata
  # determine the epithelial cells and only include these in heatmap:
  print(paste0("Number of heatmap rows before non-epithelial thrown: ", 
  	nrow(infercnv_output)))
  epithelial_ids <- epithelial_metadata$cell_ids[grep("pithelial", 
    epithelial_metadata$cell_type)]
  epithelial_heatmap <- infercnv_output[
    rownames(infercnv_output) %in% epithelial_ids,
  ]
  print(paste0("Number of heatmap rows after non-epithelial thrown: ", 
  	nrow(epithelial_heatmap)))

  # check all epithelial cells are present:
  if (!file.exists(paste0(table_dir, "integrated_object_epithelial_cell_ids.txt"))) {
    integrated_object <- readRDS(paste0(integrated_dir, 
      "/01_seurat_subset_Epithelial.Rdata"))
    ident_names <- gsub("4290", "4290A", names(Idents(integrated_object)))
    if (length(grep(sample_name, ident_names)) > 0) {
      integrated_epi_ids <- gsub("4290", "4290A", 
        names(Idents(integrated_object))[grep(sample_name, ident_names)]
      )
      write.table(integrated_epi_ids, paste0(table_dir, 
        "integrated_object_epithelial_cell_ids.txt"), quote=F, sep="\n", 
        col.names=F, row.names=F)
    }
  } else {
    integrated_epi_ids <- read.table(paste0(table_dir, 
      "integrated_object_epithelial_cell_ids.txt"), header=F, sep="\n", 
      as.is=T)[,1]
  }
  if (exists("integrated_epi_ids")) {
    missing_epis <- integrated_epi_ids[
      !(integrated_epi_ids %in% epithelial_ids)
    ]
    print(paste0(
      "Total no. epithelial cells in integrated data but not in InferCNV output = ",
      length(missing_epis)
    ))
    if (length(missing_epis) > 20) {
      system(paste0("touch ", in_dir, "over_20_epithelial_cells_missing"))
    }
  }

  saveRDS(epithelial_heatmap, paste0(Robject_dir, 
    "/1a.initial_epithelial_heatmap.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, 
    "/1b.initial_epithelial_metadata.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

} else {

  print("Loading heatmap and metadata dfs...")
  epithelial_heatmap <- readRDS(paste0(Robject_dir, 
    "/1a.initial_epithelial_heatmap.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/1b.initial_epithelial_metadata.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

}


################################################################################
### 2. Add QC metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "2b.epithelial_metadata_with_cell_type_and_QC.Rdata"))) {

  # create epithelial_metadata df and only include epithelial cells in epithelial_heatmap:
  print("Creating epithelial metadata df...")
  epithelial_metadata <- epithelial_metadata[rownames(epithelial_heatmap),]
  # order epithelial metadata cell type cluster levels:
  epithelial_metadata$cell_type <- factor(
    epithelial_metadata$cell_type,
    levels = naturalsort(unique(epithelial_metadata$cell_type))
  )
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  # add nUMI and nGene data to epithelial_metadata:
  print("Adding QC metrics to epithelial metadata df...")
  if (!exists("seurat_10X")) {
    seurat_10X <- readRDS(
      paste0(seurat_dir, "03_seurat_object_processed.Rdata")
    )
    Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  }
  QC <- data.frame(
    row.names = names(Idents(seurat_10X)),
    nUMI = seurat_10X@meta.data$nCount_RNA,
    nGene = seurat_10X@meta.data$nFeature_RNA
  )
  QC <- QC[rownames(epithelial_metadata),]
  epithelial_metadata <- cbind(epithelial_metadata, QC)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  saveRDS(
    epithelial_heatmap, paste0(Robject_dir, 
    "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata")
  )
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "2b.epithelial_metadata_with_cell_type_and_QC.Rdata")
  )

} else {

  epithelial_heatmap <- readRDS(
    paste0(Robject_dir, 
    "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata")
  )
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "2b.epithelial_metadata_with_cell_type_and_QC.Rdata")
  )
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
}


################################################################################
### 3. Add CNA and correlation value metadata ###
################################################################################

# create and check normality from density plots for average CNV vector:
CNV_average <- apply(epithelial_heatmap, 2, mean)

CNV_density_plot <- density(CNV_average, bw="SJ")
pdf(paste0(plot_dir, "average_CNV_density_plot.pdf"))
  plot(CNV_density_plot, main=NA, xlab = "CNV values")
dev.off()
png(paste0(plot_dir, "average_CNV_density_plot.png"))
  plot(CNV_density_plot, main=NA, xlab = "CNV values")
dev.off()

if (!file.exists(paste0(Robject_dir, 
  "3b.epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata"))) {
  
  # determine CNA values and add to epithelial_metadata:
  print("Determining CNA values and adding to epithelial metadata df...")
  # scale infercnv values to -1:1, square values and take the mean:
  scaled_df <- as.data.frame(rescale(as.matrix(epithelial_heatmap), c(-1,1)))
  CNA_values <- apply(scaled_df, 1, function(y) {
    #y[is.na(y)] <- 0
    #scaled_y <- rescale(y, c(-1, 1))
    return(mean(y^2))
  })
  CNA_value_df <- data.frame(
    row.names = names(CNA_values),
    CNA_value = CNA_values
  )
  epithelial_metadata <- cbind(epithelial_metadata, CNA_value_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
    # determine correlation with top 5% cancer values and add to epithelial_metadata:
  print(paste0(
    "Determining correlation with top 5% cancer values and adding to epithelial ", 
    "metadata df..."
  ))
  # determine top 5% cancer cells:
  CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
  ordered_CNA_values  <- data.frame(
    row.names = rownames(CNA_value_df)[CNA_order],
    CNA_value = CNA_value_df[CNA_order,]
  )
  top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)
  
  # find average genome-wide CNV predictions across genome:
  top_cancer_CNV_average <- apply(epithelial_heatmap[rownames(top_cancer),], 2, mean)
  # find correlations of each cell's CNVs with top_GIN_CNV_average:
  cancer_correlations <- apply(epithelial_heatmap, 1, function(x) {
    if (length(unique(as.numeric(x))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
        cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
      cor_result <- data.frame(cor$estimate, cor$p.value)
    }
    return(cor_result)
  })
  correlation_df <- do.call("rbind", cancer_correlations)

  # add to epithelial_metadata:
  epithelial_metadata <- cbind(epithelial_metadata, correlation_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  saveRDS(
    epithelial_heatmap, paste0(Robject_dir, 
    "3a.epithelial_heatmap_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "3b.epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )

} else {

  epithelial_heatmap <- readRDS(
    paste0(Robject_dir, 
    "3a.epithelial_heatmap_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )

  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "3b.epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

}


################################################################################
### 4a. Call normals and add to epithelial_metadata ###
################################################################################

# create and add normal cell annotations:
if (
  !file.exists(paste0(Robject_dir, 
    "4b.epithelial_metadata_with_cell_type_QC_CNA_correlation_and_", 
    "normal_call_values.Rdata")) | 
  !(file.exists(paste0(plot_dir,"infercnv_plot.png"))) |
  !(file.exists(paste0(plot_dir,"CNA_density_plot.png"))) |
  !(file.exists(paste0(plot_dir, 
    "normal_call_quad_plot_mean_of_scaled_squares.png")))
) {
  print(paste0(
    "Determing normal epithelial cells and adding to epithelial ",
    "metadata df..."
  ))
  # create density plot of infercnv values:
  CNA_density_plot <- density(epithelial_metadata$CNA_value, bw="SJ")
  pdf(paste0(plot_dir, "CNA_density_plot.pdf"))
    plot(CNA_density_plot, main=NA, xlab = "CNA value")
  dev.off()
  png(paste0(plot_dir, "CNA_density_plot.png"))
    plot(CNA_density_plot, main=NA, xlab = "CNA value")
  dev.off()

  # prepare df for quad plots:
  quad_df <- data.frame(
    row.names = rownames(epithelial_metadata),
    CNA_value = epithelial_metadata$CNA_value, 
    cor.estimate = epithelial_metadata$cor.estimate
  )

  if (sample_name %in% mono_samples) {
    CNA_mean <- mean(epithelial_metadata$CNA_value)
    CNA_std_dev <- sd(epithelial_metadata$CNA_value)
    cor_mean <- mean(epithelial_metadata$cor.estimate)
    cor_std_dev <- sd(epithelial_metadata$cor.estimate)
    x_int1 <- CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev)
    y_int1 <- cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev)
    normal_outliers <- rownames(epithelial_metadata)[
      epithelial_metadata$cor.estimate < y_int1
    ]
    x_int2 <- CNA_mean + (normal_x_threshold_sd_multiplier*cor_std_dev)
    y_int2 <- cor_mean + (normal_y_threshold_sd_multiplier*cor_std_dev)
    cancer_outliers <- rownames(epithelial_metadata)[
      epithelial_metadata$CNA_value > x_int2 & epithelial_metadata$cor.estimate > y_int2
    ]
    if (length(normal_outliers) >= length(cancer_outliers)) {
      x_int <- round(x_int1, 3)
      y_int <- round(y_int1, 3)
    } else if (length(normal_outliers) < length(cancer_outliers)) {
      x_int <- round(x_int2, 3)
      y_int <- round(y_int2, 3)
    }

    if (adjust_normal_thresholds) {
      if (sample_name %in% adjustment_df$sample_id) {
        x_int <- adjustment_df$x_int[adjustment_df$sample_id == sample_name]
        y_int <- adjustment_df$y_int[adjustment_df$sample_id == sample_name]
      }
    }

    # define normal and cancer cells:
    epithelial_metadata$normal_cell_call <- "cancer"
    if (sample_name %in% assign_unassigned) {

      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
      ] <- "normal"
      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value > x_int | epithelial_metadata$cor.estimate > y_int
      ] <- "cancer"

    } else if (above_y_axis_only) {

      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
      ] <- "normal"
      epithelial_metadata$normal_cell_call[
        epithelial_metadata$cor.estimate > y_int
      ] <- "cancer"
      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
      ] <- "unassigned"

    } else {  

      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
      ] <- "normal"
      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int
      ] <- "unassigned"
      epithelial_metadata$normal_cell_call[
        epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
      ] <- "unassigned"
    }
    
    # create quad plot:
    if (!("cancer" %in% unique(epithelial_metadata$normal_cell_call))) {
      p <- ggplot(epithelial_metadata, 
                  aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("#74add1", "#b2182b"), 
                                    labels=c("Normal", "Unassigned"))
      p <- p + xlab("CNA level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    } else if (!("normal" %in% unique(epithelial_metadata$normal_cell_call))) {
      p <- ggplot(epithelial_metadata, 
                  aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#b2182b"), 
                                    labels=c("Cancer", "Unassigned"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    } else if (!("unassigned" %in% unique(epithelial_metadata$normal_cell_call))) {
      p <- ggplot(epithelial_metadata, 
        aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#74add1"), 
                                    labels=c("Cancer", "Normal"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    } else if (length(unique(epithelial_metadata$normal_cell_call)) == 3) {
      p <- ggplot(epithelial_metadata, 
          aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                    labels=c("Cancer", "Normal", "Unassigned"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    }
    png(paste0(plot_dir, 
      "normal_call_quad_plot_mean_of_scaled_squares.png"), 
      width = 450, height = 270)
      print(p)
    dev.off()
    p
    ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares.pdf"),
      width = 18, height = 13, units = c("cm"))
    dev.off()


  ################################################################################
  ### 4b. Call normals and add to epithelial_metadata ###
  ################################################################################

  } else if (sample_name %in% multi_samples) {
    # scale data:
    scaled_quad_df <- scale(quad_df) %>% as.data.frame()
    # run silhouette cluster analysis to determine clusters and thresholds:
    pamk_result <- pamk(scaled_quad_df, krange=1:4)
    pamk_result$nc
    silhouette_result <- pam(scaled_quad_df, pamk_result$nc)
    saveRDS(silhouette_result, paste0(Robject_dir, "silhouette_result.Rdata"))

    # if no. clusters estimated to be > 1, use cluster information to set 
    # normal vs cancer thresholds:
    if (pamk_result$nc > 1) {

      sil_values <- as.data.frame(silhouette_result$silinfo$widths)
      sil_result <- data.frame(row.names=names(silhouette_result$clustering),
        cluster=silhouette_result$clustering,
        sil_width=sil_values$sil_width)

      # add sil_result to epithelial_metadata:
      epithelial_metadata <- cbind(epithelial_metadata, sil_result)
      
      # determine normal and cancer clusters by determining the max CNA values and
      # correlation with top 5% cancer:
      cluster_split <- split(epithelial_metadata, epithelial_metadata$cluster)
      names(cluster_split) <- paste0("cluster_", names(cluster_split))
      # determine order of clusters by adding mean CNA and correlation values:
      cluster_means <- sort(
        unlist(
          lapply(cluster_split, function(x) {
            return(mean(x$CNA_value) + mean(x$cor.estimate))
          })
        )
      )
      # determine second cluster from axes as cancer cluster closest to axes:
      first_cancer_cluster <- names(cluster_means[2])
      first_cancer_df <- eval(parse(text=paste0("cluster_split$", first_cancer_cluster)))
      # make intercepts 1 std dev from mean towards axes:
      # define x-axis as 2 std devs left of mean:
      CNA_mean <- mean(first_cancer_df$CNA_value)
      CNA_std_dev <- sd(first_cancer_df$CNA_value)
      x_int <- round(CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev), 3)
      # define x-axis as 2 std devs left of mean:
      cor_mean <- mean(first_cancer_df$cor.estimate)
      cor_std_dev <- sd(first_cancer_df$cor.estimate)
      y_int <- round(cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev), 3)

      if (adjust_normal_thresholds) {
        if (sample_name %in% adjustment_df$sample_id) {
          x_int <- adjustment_df$x_int[adjustment_df$sample_id == sample_name]
          y_int <- adjustment_df$y_int[adjustment_df$sample_id == sample_name]
        }
      }

      # define normal and cancer cells:
      epithelial_metadata$normal_cell_call <- "cancer"
      if (sample_name %in% assign_unassigned) {
        epithelial_metadata$normal_cell_call <- "cancer"
        epithelial_metadata$normal_cell_call[
          epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
        ] <- "normal"
        epithelial_metadata$normal_cell_call[
          epithelial_metadata$CNA_value > x_int | epithelial_metadata$cor.estimate > y_int
        ] <- "cancer"
      } else {
        epithelial_metadata$normal_cell_call[
          epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
        ] <- "normal"
        epithelial_metadata$normal_cell_call[
          epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int
        ] <- "unassigned"
        epithelial_metadata$normal_cell_call[
          epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
        ] <- "unassigned"
      }
      
      # create quad plot:
      p <- ggplot(epithelial_metadata, 
                  aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                    labels=c("Cancer", "Normal", "Unassigned"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
      quad_plot <- p
      png(paste0(plot_dir, 
          "normal_call_quad_plot_mean_of_scaled_squares.png"), 
          width = 450, height = 270)
          print(quad_plot)
      dev.off()
      quad_plot
      ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares.pdf"))
      dev.off()
      # create quad plot with clusters marked:
      p <- ggplot(epithelial_metadata, 
        aes(x=CNA_value, y=cor.estimate, color=as.factor(cluster)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("#9B59B6", "#b8e186", "#18ffff", "#ef6c00"), 
        labels=c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"))
      p <- p + xlab("CNA value")
      p <- p + ylab("Corr. with top 5% cancer")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
      quad_plot_clusters <- p
      png(paste0(plot_dir, 
        "normal_call_quad_plot_mean_of_scaled_squares_clusters.png"), width = 430, height = 250)
        print(quad_plot_clusters)
      dev.off()
      quad_plot_clusters
      ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares_clusters.pdf"))
      dev.off()
    } else {
      print(
        paste0(
          "Sample ", sample_name, " needs to be recorded as mono-cluster and script rerun!"
        )
      )
    }
  }
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  saveRDS(epithelial_heatmap, paste0(Robject_dir, 
      "4a.epithelial_heatmap_with_cell_type_QC_CNA_correlation_and_", 
      "normal_call_values.Rdata"))
    saveRDS(epithelial_metadata, paste0(Robject_dir, 
      "4b.epithelial_metadata_with_cell_type_QC_CNA_correlation_and_", 
      "normal_call_values.Rdata"))

} else {

  epithelial_heatmap <- readRDS(paste0(Robject_dir, 
    "4a.epithelial_heatmap_with_cell_type_QC_CNA_correlation_and_", 
    "normal_call_values.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "4b.epithelial_metadata_with_cell_type_QC_CNA_correlation_and_", 
    "normal_call_values.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

}


################################################################################
### 5. Add sc50 calls to metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "5b.epithelial_metadata_with_cell_type_QC_CNA_correlation_", 
    "normal_call_values_and_sc50.Rdata"))) {
  print(paste0(
    "Adding sc50 calls to epithelial metadata df..."
  ))
  # read in data and select data for current sample:
  sc50 <- read.table(paste0(ref_dir, "/sc50_calls.txt"), header = T, sep = "\t",
    stringsAsFactors = F)
  sc50$SampleID <- gsub("4290", "4290A", sc50$SampleID)
  sample_sc50 <- sc50[grep(sample_name, sc50$SampleID),]

  if (length(sample_sc50$Sept.SC50_Calls) > 0) {

    # select only cells contained in epithelial_metadata:
    print(paste0("Original sc50 call no. for ", sample_name, " = ", 
      nrow(sample_sc50)))
    rownames(sample_sc50) <- sample_sc50$SampleID
    sample_sc50 <- sample_sc50[rownames(epithelial_metadata),]
     print(paste0("sc50 call no. for ", sample_name, 
      " after selecting only cells in epithelial metadata df = ", nrow(sample_sc50)))

    # cbind to epithelial_metadata:
    epithelial_metadata <- cbind(epithelial_metadata, sample_sc50$Sept.SC50_Calls)
    colnames(epithelial_metadata)[ncol(epithelial_metadata)] <- "sc50"

  } else {

    print(paste0(
      "No sc50 calls found for sample ", sample_name, ", adding NAs for this column..."
    ))
    epithelial_metadata$sc50 <- NA
  }

  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  saveRDS(epithelial_heatmap, paste0(Robject_dir, 
    "5a.epithelial_heatmap_final.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, 
    "5b.epithelial_metadata_final.Rdata"))

} else {

  epithelial_heatmap <- readRDS(paste0(Robject_dir, 
    "5a.epithelial_heatmap_with_cell_type_QC_CNA_correlation_", 
    "normal_call_values_and_sc50.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "5b.epithelial_metadata_with_cell_type_QC_CNA_correlation_", 
    "normal_call_values_and_sc50.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

}


################################################################################
### 6. Order heatmap, metadata and create heatmap annotations ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "6b.epithelial_metadata_final.Rdata"))) {
  # reorder cells starting with normals, unassigned and ending with cancer:
  epithelial_metadata_split <- split(epithelial_metadata, 
    epithelial_metadata$normal_cell_call)
  epithelial_metadata <- do.call("rbind",
    list(
      epithelial_metadata_split$normal,
      epithelial_metadata_split$unassigned,
      epithelial_metadata_split$cancer
    )
  )
  epithelial_heatmap <- epithelial_heatmap[rownames(epithelial_metadata),]
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  saveRDS(epithelial_heatmap, paste0(Robject_dir, "6a.epithelial_heatmap_final.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, "6b.epithelial_metadata_final.Rdata"))

} else {

  epithelial_heatmap <- readRDS(paste0(Robject_dir, "6a.epithelial_heatmap_final.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, "6b.epithelial_metadata_final.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

}

# define group annotation colours:
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
cluster_number <- length(unique(epithelial_metadata$cell_type))
cluster_cols <- col_palette[1:cluster_number]
names(cluster_cols) <- unique(epithelial_metadata$cell_type)
cluster_cols <- cluster_cols[levels(epithelial_metadata$cell_type)]

# create group annotation:
group_annotation_df <- subset(epithelial_metadata, select = cell_type)
group_annotation <- Heatmap(
  group_annotation_df, 
  col = cluster_cols, 
  name = "group_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  show_heatmap_legend = F
)
# create CNA annotation:
CNA_value_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$CNA_value,
    gp = gpar(
      col = "#D95F02", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
CNA_value_annotation@name <- "CNA_value"
# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$nUMI,
    gp = gpar(
      col = "#D8B72E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nUMI_annotation@name <- "nUMI"
nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$nGene, name = "nGene",
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nGene_annotation@name <- "nGene"
# create normal call annotation:
normal_call_annot_df <- subset(epithelial_metadata, select = normal_cell_call)
normal_call_annotation <- Heatmap(normal_call_annot_df, 
  col = c("unassigned" = "#E7E4D3", "normal" = "#1B7837", "cancer" = "#E7298A"), 
  name = "normal_call_annotation", width = unit(6, "mm"), 
  show_row_names = F, show_column_names = F, 
  show_heatmap_legend = F
)

# create array CNV annotation:
all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
colnames(all_array_CNVs) <- gsub("CID4499_1", "CID44991", colnames(all_array_CNVs))
if (any(colnames(all_array_CNVs) %in% sample_name)) {
  if (!file.exists(paste0(Robject_dir, "array_CNV_annotation.Rdata"))) {
    array_CNV_annotation <- create_array_CNV_annotation(epithelial_heatmap, all_array_CNVs)
    saveRDS(array_CNV_annotation, paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  } else {
    array_CNV_annotation <- readRDS(paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  }

  # remove normals and unassigned from epithelial_heatmap:
  cancer_only_heatmap <- epithelial_heatmap[
    epithelial_metadata$cell_ids[epithelial_metadata$normal_cell_call == "cancer"],
  ]

  # correlate infercnv CNV profiles with array profiles:
  average_heatmap <- apply(cancer_only_heatmap, 2, mean)

  cor_pearson <- cor.test(
  	as.numeric(array_CNV_annotation$array_CNVs), 
  	as.numeric(average_heatmap), 
  	method = "pearson",
  	alternative = "greater"
  )
  CNV_correlation_result <- data.frame(cor_pearson$estimate, cor_pearson$p.value)

}


################################################################################
### 7. Generate heatmap ###
################################################################################

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}
# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")

print("Generating final heatmap...")
# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = F,
  heatmap_legend_param = list(labels_gp = gpar(col = "red", fontsize = 12)),
  use_raster = T, raster_device = c("png")
)

if (include_group_annotation) {
  ht_list <- group_annotation + final_heatmap + normal_call_annotation + 
    CNA_value_annotation + nUMI_annotation + nGene_annotation
} else {
  ht_list <- final_heatmap + normal_call_annotation + 
    CNA_value_annotation + nUMI_annotation + nGene_annotation
}

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 8. Plot heatmap ###
################################################################################

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot_non_rounded_cor_p_values.pdf"), 
	height = 13, width = 20) 

  grid.newpage()
  pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
  	just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
        	gp = gpar(lwd = 1, col = "#383838"))
        }
      })
    popViewport()

    pushViewport(viewport(x=x_coord+0.13, y=0.5, width = 0.2, height = 0.1, 
      just = "right"))
      grid.text(sample_name, gp = gpar(fontsize=50))
    popViewport()

    if (exists("grid_array_heatmap")) {
      pushViewport(viewport(x = x_coord + 0.7755, y = 0.86, 
     	width = 0.673, height = 0.13, just = c("right", "bottom")))
        grid.draw(grid_array_heatmap)
      popViewport()

      pushViewport(viewport(x=x_coord + 0.754, y=0.945, width = 0.2, height = 0.1, 
        just = "left"))
        grid.text(
          paste0(
          	"R2 = ", CNV_correlation_result$cor_pearson.estimate
          ), gp = gpar(fontsize=30)
        )
      popViewport()
      pushViewport(viewport(x=x_coord + 0.7015, y=0.91, width = 0.2, height = 0.1, 
        just = "left"))
        grid.text(
          paste0(
          	"p = ", CNV_correlation_result$cor_pearson.p.value
          ), gp = gpar(fontsize=30)
        )
      popViewport()

    }
    
dev.off()

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot.pdf"), 
  height = 13, width = 20) 

  grid.newpage()
    pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
      just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
          gp = gpar(lwd = 1, col = "#383838"))
        }
      })
    popViewport()

    pushViewport(viewport(x=x_coord+0.13, y=0.5, width = 0.2, height = 0.1, 
      just = "right"))
      grid.text(sample_name, gp = gpar(fontsize=50))
    popViewport()

    if (exists("grid_array_heatmap")) {
      pushViewport(viewport(x = x_coord + 0.7755, y = 0.86, 
      width = 0.6725, height = 0.13, just = c("right", "bottom")))
        grid.draw(grid_array_heatmap)
      popViewport()

      pushViewport(viewport(x=x_coord + 0.728, y=0.945, width = 0.2, height = 0.1, 
        just = "left"))
        grid.text(
          paste0(
            "R2 = ", round(CNV_correlation_result$cor_pearson.estimate, 3)
          ), gp = gpar(fontsize=30)
        )
      popViewport()
      pushViewport(viewport(x=x_coord + 0.7015, y=0.91, width = 0.2, height = 0.1, 
        just = "left"))
        grid.text(
          paste0(
            "p = ", round(CNV_correlation_result$cor_pearson.p.value, 3)
          ), gp = gpar(fontsize=30)
        )
      popViewport()
    }
    
dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

# print epithelial_metadata as table:
write.table(epithelial_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))


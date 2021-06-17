# ==================================
# title: "Figure 8"
# project: "brca-atlas paper (2021)"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# date: 17/08/2020
# modified: 10/06/2021
# ==================================

# Overview

# Generates Figure 8 for brca-atlas paper (2021).

# ======================
# FUNCTIONS
# ======================
# ---------------------
# Nice function to re-position legend in blank facet panels
# ---------------------
# from here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
library(gtable)
library(cowplot)

shift_legend <- function(p)
{
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}
# ---------------------


# ==============
# Config
# ==============
message("CONFIG")

# Libraries
suppressPackageStartupMessages({
  source("code/R/plot_utils.R")
  source("code/R/survival_utils.R")  
  library(ComplexHeatmap)
  library(RColorBrewer) 
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggsignif)
  library(randomcoloR)
  library(survival)
  library(survminer)
})

# INPUT

# OUTPUT
output_directory <- "paper/figures/figure_8/"
if(!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# ============================
# Colour Palettes
# ============================

# For ComplexHeatmap
dfColourCodes <- data.frame(name = c("Basal", "Her2", "Her2E", "LumA", "LumB", "Normal", "NC", 
                                     "LumB_Her2", "HER2_ER", "ER", "TNBC", "HER2", "metaplastic",
                                     "Her2_Amp", "Her2_Non_Amp", "TBC",
                                     "LuminalA", "LuminalB", "LuminalB_HER2",
                                     "Normal_breast",
                                     "LumA_SC", "LumB_SC", "Her2_SC", "Basal_SC",                                     
                                     "cancer", "unassigned", "normal",
                                     "METABRIC_discovery", "METABRIC_validation"),
                            value = c("red", "pink", "pink", "darkblue", "cyan", "green", "grey50", 
                                      "purple", "purple", "blue", "red", "pink", "yellow", 
                                      "black", "lightgrey", "yellow",
                                      "darkblue", "cyan", "pink",
                                      "grey50",
                                      "darkblue", "cyan", "pink", "red",
                                      "black", "grey50", "green",
                                      "#7BD18C", "#FFB0BD")
                            )

dfCellTypeColourCodes <- data.frame(name = c("Cancer", "Normal", "T-cells & ILCs", "T-cells", "ILCs", "B-cells", "Plasmablasts", "Myeloid", "CAFs", "PVLs", "Endothelial", "Cycling"),
                                    value = c("#CC9933", "#FFCC66", "#CC6666", "#CC6666", "#CC6666", "#9999FF", "#66CC99", "#66CC66", "#999900", "#6699FF", "#66CCCC", "grey50"))

# For ggplot
subtype_colours_ggplot = c("LumA" = "darkblue",
                           "LumB" = "cyan", 
                           "Her2E" = "pink",
                           "Basal" = "red",
                           "Normal" = "green",
                           "NC" = "grey50")

celltype_colours_ggplot <- c("Cycling" = "grey50",
                             "Cancer" = "#CC9933", 
                            "LumA_SC" = "darkblue",
                            "LumB_SC" = "cyan", 
                            "Her2E_SC" = "pink",
                            "Basal_SC" = "red",
                            "Normal" = "#FFCC66", 
                            "T-cells & ILCs" = "#CC6666",
                            "B-cells" = "#9999FF",
                            "Plasmablasts" = "#66CC99",
                            "Myeloid" = "#66CC66",
                            "CAFs" = "#999900",
                            "PVLs" = "#6699FF",
                            "Endothelial" = "#66CCCC")


# PARAMS
analysis_id <- "Jul2020"
celltype_classification <- "celltype_subset"
normal_epithelial_cells <- "with_normal_epithelial_cells"
cancer_cell_annotation <- "SCTyper_Dec"
cycling_cell_annotation <- "combined"
subsample_fraction <- 0.15
batch_mode <- "Smode"

# PARSE METABRIC cell-fractions
# --------
# discovery
file_metabric_discovery_results <- paste0("analysis/CIBERSORTx/",
                                          analysis_id, "/",
                                          celltype_classification, "/",
                                          normal_epithelial_cells, "/",
                                          cancer_cell_annotation, "/",
                                          "cycling_", cycling_cell_annotation, "/",
                                          "sampled_", subsample_fraction, "/",
                                          "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/discovery/",
                                          "rmbatch", batch_mode, "/",
                                          "refsample_counts/QN_FALSE/perm_0/relative/CIBERSORTx_Adjusted.txt")

# validation
file_metabric_validation_results <- paste0("analysis/CIBERSORTx/",
                                           analysis_id, "/",
                                           celltype_classification, "/",
                                           normal_epithelial_cells, "/",
                                           cancer_cell_annotation, "/",
                                           "cycling_", cycling_cell_annotation, "/",
                                           "sampled_", subsample_fraction, "/",
                                           "cell_fractions/single_cell_TRUE/Gmin_300/Gmax_500/qvalue_0.01/kmax_999/replicates_5/sampling_0.5/fraction_0.75/metabric/validation/",
                                           "rmbatch", batch_mode, "/",
                                           "refsample_counts/QN_FALSE/perm_0/relative/CIBERSORTx_Adjusted.txt")

# meta-data
file_metabric_clinical <- "data/METABRIC/brca_metabric_clinical_data.tsv"
file_metabric_PAM50 <- "data/METABRIC/heloisa/METABRIC_ClassLabels_PAM50.txt"

# PRE-PROCESS
df_metabric_cell_fractions <- NULL

df_metabric <- read.delim(file_metabric_discovery_results, stringsAsFactors = F)
colnames(df_metabric)[1] <- "SampleID"
df_metabric[, "COHORT"] <- "METABRIC_discovery"
df_metabric_cell_fractions <- bind_rows(df_metabric_cell_fractions, df_metabric)

df_metabric <- read.delim(file_metabric_validation_results, stringsAsFactors = F)
colnames(df_metabric)[1] <- "SampleID"
df_metabric[, "COHORT"] <- "METABRIC_validation"
df_metabric_cell_fractions <- bind_rows(df_metabric_cell_fractions, df_metabric)

dfMETABRIC_PAM50 <- read.delim(file_metabric_PAM50, stringsAsFactors = FALSE, comment.char = "#")
colnames(dfMETABRIC_PAM50)[2] <- "SUBTYPE"

# modify the PAM50 labels
dfMETABRIC_PAM50[dfMETABRIC_PAM50$SUBTYPE == "Bas", "SUBTYPE"] <- "Basal"
dfMETABRIC_PAM50[dfMETABRIC_PAM50$SUBTYPE == "Her2", "SUBTYPE"] <- "Her2E"
dfMETABRIC_PAM50[dfMETABRIC_PAM50$SUBTYPE == "Norm", "SUBTYPE"] <- "Normal"

df_metabric_cell_fractions <- dplyr::left_join(df_metabric_cell_fractions,
                                               dfMETABRIC_PAM50,
                                               by = c("SampleID" = "METABRIC_ID"))

cell_type_order <- c("Cancer_LumA_SC", "Cancer_LumB_SC", "Cancer_Her2_SC", "Cancer_Basal_SC",
                     "Cycling",
                     "Myoepithelial", "Luminal_Progenitors", "Mature_Luminal",
                     "CAFs_MSC_iCAFlike_s1", "CAFs_MSC_iCAFlike_s2", "CAFs_myCAF_like_s4", "CAFs_myCAF_like_s5", "CAFs_Transitioning_s3",
                     "Endothelial_ACKR1", "Endothelial_CXCL12", "Endothelial_Lymphatic_LYVE1", "Endothelial_RGS5",
                     "PVL_Immature_s1", "PVL_Immature_s2", "PVL_Differentiated_s3",
                     "T_cells_c0_CD4_CCR7", "T_cells_c1_CD4_IL7R", "T_cells_c2_CD4_Tregs_FOXP3", "T_cells_c3_CD4_Tfh_CXCL13",
                     "T_cells_c4_CD8_ZFP36", "T_cells_c5_CD8_GZMK",  "T_cells_c7_CD8_IFNG", "T_cells_c8_CD8_LAG3",
                     "T_cells_c6_IFIT1", "T_cells_c9_NK_cells_AREG", "T_cells_c10_NKT_cells_FCGR3A",
                     "B_cells_Naive", "B_cells_Memory", "Plasmablasts",
                     "Myeloid_c0_DC_LAMP3", "Myeloid_c1_LAM1_FABP5", "Myeloid_c2_LAM2_APOE", "Myeloid_c3_cDC1_CLEC9A", "Myeloid_c4_DCs_pDC_IRF7", "Myeloid_c5_Macrophage_3_SIGLEC1", "Myeloid_c7_Monocyte_3_FCGR3A", "Myeloid_c8_Monocyte_2_S100A9", "Myeloid_c9_Macrophage_2_CXCL10", "Myeloid_c10_Macrophage_1_EGR1", "Myeloid_c11_cDC2_CD1C", "Myeloid_c12_Monocyte_1_IL1B")
 
# ============================
figNum <- "8A"
message(paste0("FIGURE ", figNum))
# ============================
message("Heatmap of consensus clustered CIBERSORTx deconvolution results across all METABRIC")

# get the eco-type groupings
cohort <- "metabric/combined"
batch_correction <- "smode"
clustering_method <- "consensus"
cluster_cohort_name <- "all"
cluster_file <- paste0("sample_clusters.", cluster_cohort_name, ".txt")

# ecotype groupings
file_sample_clusters <- paste0("analysis/CIBERSORTx/analysis/", 
                               analysis_id, "/",
                               celltype_classification, "/",
                               normal_epithelial_cells, "/",
                               cancer_cell_annotation, "/",
                               "cycling_", cycling_cell_annotation, "/",                                       
                               "sampled_", subsample_fraction, "/",
                               cohort,
                               "/rmbatch_", batch_correction,
                               "/clustering/", clustering_method,
                               "/", cluster_file)

# cell-type groupings
file_celltype_clusters <- paste0("analysis/CIBERSORTx/analysis/", 
                                 analysis_id, "/",
                                 celltype_classification, "/",
                                 normal_epithelial_cells, "/",
                                 cancer_cell_annotation, "/",
                                 "cycling_", cycling_cell_annotation, "/",                               
                                 "sampled_", subsample_fraction, "/",
                                 cohort,
                                 "/rmbatch_", batch_correction,
                                 "/clustering/", clustering_method,
                                 "/", "celltype_clusters.txt")

# Add ecotypes
df_ecotypes <- read.delim(file_sample_clusters, stringsAsFactors = F)
colnames(df_ecotypes)[1] <- "Ecotype"
df_ecotypes$Ecotype <- paste0("E", df_ecotypes$Ecotype)
df_ecotype_cell_fractions <- dplyr::inner_join(df_metabric_cell_fractions, df_ecotypes[, c("SAMPLE_ID", "Ecotype")], by = c("SampleID" = "SAMPLE_ID"))

dfPlot <- df_ecotype_cell_fractions
rownames(dfPlot) <- dfPlot$SampleID

# get cell-type clusters
df_celltypes <- read.delim(file_celltype_clusters, stringsAsFactors = F)

# prepare matrix
dfPlot <- dfPlot[df_ecotypes$SAMPLE_ID, ]
matPlot <- dfPlot[, df_celltypes$celltype]
matPlot <- t(matPlot)
# z-score
matPlot <- scale_rows_zscore(matPlot)

# heatmap
cluster_method <- "average"
cluster_distance <- "euclidean"

meta_data_cols <- rev(c("SUBTYPE", "COHORT", "Ecotype"))
dfColAnno <- dfPlot[, meta_data_cols, drop = F]
colnames(dfColAnno) <- c("Tumour Ecotype", "METABRIC Cohort", "Tumour Subtype")

# row-annotation
# major cell-type
dfRowAnno <- data.frame(major_celltype = rownames(matPlot), stringsAsFactors = F)
rownames(dfRowAnno) <- dfRowAnno$major_celltype
dfRowAnno[grepl(pattern = "Cancer", x = dfRowAnno$major_celltype), 1] <- "Cancer"
dfRowAnno[grepl(pattern = "Myoepithelial", x = dfRowAnno$major_celltype), 1] <- "Normal"
dfRowAnno[grepl(pattern = "Luminal_Progenitor", x = dfRowAnno$major_celltype), 1] <- "Normal"
dfRowAnno[grepl(pattern = "Mature_Luminal", x = dfRowAnno$major_celltype), 1] <- "Normal"
dfRowAnno[grepl(pattern = "T_cells", x = dfRowAnno$major_celltype), 1] <- "T-cells & ILCs"
dfRowAnno[grepl(pattern = "B_cells", x = dfRowAnno$major_celltype), 1] <- "B-cells"
dfRowAnno[grepl(pattern = "Myeloid", x = dfRowAnno$major_celltype), 1] <- "Myeloid"
dfRowAnno[grepl(pattern = "CAFs", x = dfRowAnno$major_celltype), 1] <- "CAFs"
dfRowAnno[grepl(pattern = "PVL", x = dfRowAnno$major_celltype), 1] <- "PVLs"
dfRowAnno[grepl(pattern = "Endothelial", x = dfRowAnno$major_celltype), 1] <- "Endothelial"
colnames(dfRowAnno)[1] <- "Major Celltype"

names_vector <- colnames(dfRowAnno)
values_vector <- c(rep(list(tibble::deframe(dfCellTypeColourCodes)), times = length(colnames(dfRowAnno))))
colour_map_list <- as.list(setNames(values_vector, names_vector))
rowAnno <- rowAnnotation(df = dfRowAnno, col = colour_map_list, simple_anno_size = unit(0.3, "cm"), show_annotation_name = F)

# group split
row_split <- df_celltypes$cluster

# add colours  
set.seed(123)
# NOTE: This isn't giving reproducible colour palletes:
#dfEcotypeColours <- data.frame(name = as.factor(unique(dfColAnno$`Tumour Ecotype`)), value = distinctColorPalette(length(unique(dfColAnno$`Tumour Ecotype`))))
# NOTE: so going to hardcode ecotype colours for the figures
dfEcotypeColours <- data.frame(name = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9"),
                               value = c("#CC99CC", "#99CCCC", "#99CC99", "#CC6666", "#9966CC", "#99CC66", "#CCCCCC", "#9999CC", "#CCCC66")
)
dfColourCodes <- rbind(dfColourCodes, dfEcotypeColours)

set.seed(123)
ht_consensus <- complexHeatmap_wrapper(matData = matPlot, 
                                       dfColAnno = dfColAnno, 
                                       dfAnnotationColourCodes = dfColourCodes, 
                                       min_max = c(-2, 2), 
                                       scale_rows = FALSE, 
                                       name = "Cell-type fraction\n(z-score)",
                                       show_column_names = F,
                                       cluster_columns = T,
                                       clustering_distance_columns = cluster_distance,
                                       clustering_method_columns = cluster_method,
                                       cluster_rows = T,
                                       clustering_distance_rows = cluster_distance,
                                       clustering_method_rows = cluster_method,                             
                                       row_names_gp = gpar(fontsize = 9.5), 
                                       row_names_side = "right",
                                       column_split=dfColAnno$`Tumour Ecotype`,
                                       row_split=row_split,
                                       row_title="cell-types",
                                       right_annotation=rowAnno
)

file_out <- paste0(output_directory, figNum, ".pdf")
fExt <- tools::file_ext(file_out)

height <- 18
width <- 51
if(fExt == "png") {
  png(file = file_out, height = height, width = width, units = "cm", res = 600)
} else if(fExt == "pdf") {
  pdf(file = file_out, height = height*0.39, width = width*0.39)  
} else {
  stop(paste0("Unsupported file type for heatmap: ", fExt))
}
set.seed(123)      
ht_draw <- draw(ht_consensus)
dev.off()

# ============================
figNum <- "8B"
message(paste0("FIGURE ", figNum))
# ============================
message("PAM50 subtype composition of Tumour Ecotypes")

dfPlot <- reshape2::melt(df_ecotype_cell_fractions, id.vars = c("SampleID", "Ecotype"), measure.vars = "SUBTYPE")
dfPlot <- dfPlot %>% 
  dplyr::mutate(value=fct_relevel(value, c("LumA", "LumB", "Her2E", "Basal", "Normal", "NC")))

# Get the proportions (see: https://stackoverflow.com/questions/24576515/relative-frequencies-proportions-with-dplyr/24576703)
dfPlot <- dfPlot %>%
  dplyr::count(Ecotype, value) %>%
  dplyr::group_by(Ecotype) %>%
  dplyr::mutate(freq = n / sum(n))

# stacked barplot
gg <- ggplot(dfPlot, aes(x = Ecotype, fill = value)) +
  geom_bar(aes(y = freq), stat = "identity", position = "stack", colour="black") +   
  scale_y_continuous(labels = scales::percent) + 
  xlab("Ecotype") +
  ylab("Fraction of tumours") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "right", 
        legend.title = element_text(face="bold", size=20),
        legend.text = element_text(size=20),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20),        
        axis.title = element_text(face="bold", size=22)) +
  scale_fill_manual(name = "Tumour Subtype",
                    breaks = names(subtype_colours_ggplot),
                    values = subtype_colours_ggplot,
                    aesthetics = c("fill"))
#print(gg)

fPlot <- paste0(output_directory, figNum, ".pdf")
ggsave(filename = fPlot, plot = gg, height = 5, width = 10)

# ============================
figNum <- "8C"
message(paste0("FIGURE ", figNum))
# ============================
message("Cell-type composition of Tumour Ecotypes")

dfPlot <- reshape2::melt(df_ecotype_cell_fractions, id.vars = c("SampleID", "Ecotype"), measure.vars = cell_type_order)
colnames(dfPlot)[3] <- "subset_celltype"

# get median fraction of cell-types in each ecotype
dfSummary <- dfPlot %>%
  dplyr::group_by(subset_celltype, Ecotype) %>%
  dplyr::summarise(median = median(value))

dfSummary$subset_celltype <- as.character(dfSummary$subset_celltype)
dfSummary$major_celltype <- dfSummary$subset_celltype
dfSummary[grepl(pattern = "Cancer", x = dfSummary$subset_celltype), "major_celltype"] <- "Cancer"
dfSummary[grepl(pattern = "Myoepithelial", x = dfSummary$subset_celltype), "major_celltype"] <- "Normal"
dfSummary[grepl(pattern = "Luminal_Progenitor", x = dfSummary$subset_celltype), "major_celltype"] <- "Normal"
dfSummary[grepl(pattern = "Mature_Luminal", x = dfSummary$subset_celltype), "major_celltype"] <- "Normal"
dfSummary[grepl(pattern = "T_cells", x = dfSummary$subset_celltype), "major_celltype"] <- "T-cells & ILCs"
dfSummary[grepl(pattern = "B_cells", x = dfSummary$subset_celltype), "major_celltype"] <- "B-cells"
dfSummary[grepl(pattern = "Myeloid", x = dfSummary$subset_celltype), "major_celltype"] <- "Myeloid"
dfSummary[grepl(pattern = "CAFs", x = dfSummary$subset_celltype), "major_celltype"] <- "CAFs"
dfSummary[grepl(pattern = "PVL", x = dfSummary$subset_celltype), "major_celltype"] <- "PVLs"
dfSummary[grepl(pattern = "Endothelial", x = dfSummary$subset_celltype), "major_celltype"] <- "Endothelial"

# get sum of median fraction of cell-types in each major cell-type in each ecotype
dfSummary_sum <- dfSummary %>%
  dplyr::group_by(major_celltype, Ecotype) %>%
  dplyr::summarise(median_sum = sum(median))

# re-scale to relative proportions
dfSummary_sum <- dfSummary_sum %>%
  dplyr::group_by(Ecotype) %>%
  dplyr::mutate(median_sum_scaled = median_sum/sum(median_sum))  

major_celltype_order <- c("Cycling", "Cancer", "Normal", "CAFs", "PVLs", "Endothelial", "Myeloid", "Plasmablasts", "B-cells", "T-cells & ILCs")
dfSummary_sum <- dfSummary_sum %>% 
  dplyr::mutate(major_celltype=fct_relevel(major_celltype, major_celltype_order))

gg <- ggplot(dfSummary_sum, aes(x = Ecotype, fill = major_celltype)) +
  geom_bar(aes(y = median_sum_scaled), stat = "identity", position = "stack", colour="black") +   
  xlab("Ecotype") +
  ylab("Cell-type proportion") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "right", 
        legend.title = element_text(face="bold", size=20),
        legend.text = element_text(size=20),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20),        
        axis.title = element_text(face="bold", size=22)) +
  scale_fill_manual(name = "Cell-type",
                    breaks = major_celltype_order,
                    values = celltype_colours_ggplot,
                    aesthetics = c("fill"))
#print(gg)

fPlot <- paste0(output_directory, figNum, ".pdf")
ggsave(filename = fPlot, plot = gg, height = 5, width = 10)

# ============================
figNum <- "8D"
message(paste0("FIGURE ", figNum))
# ============================
message("KM plot of all Ecotypes")

# OS end-time (months)
end_time <- 240

# INPUT
file_metabric_clinical <- "data/METABRIC/brca_metabric_clinical_data.tsv"
file_metabric_PAM50 <- "data/METABRIC/heloisa/METABRIC_ClassLabels_PAM50.txt"

# LOAD COHORT DATA
clin_METABRIC = read.delim(file_metabric_clinical, stringsAsFactors = F)

clin_METABRIC$Overall.Survival.Status <- gsub(":LIVING","",clin_METABRIC$Overall.Survival.Status)
clin_METABRIC$Overall.Survival.Status <- gsub(":DECEASED","",clin_METABRIC$Overall.Survival.Status)
clin_METABRIC$Overall.Survival.Status <- as.numeric(clin_METABRIC$Overall.Survival.Status)
clin_METABRIC$Overall.Survival..Months. <- as.numeric(clin_METABRIC$Overall.Survival..Months.)

# clean up column names
colnames(clin_METABRIC) <- gsub("Overall.Survival..Months.","OS",colnames(clin_METABRIC))
colnames(clin_METABRIC) <- gsub("Overall.Survival.Status","EVENT",colnames(clin_METABRIC))
colnames(clin_METABRIC) <- gsub("Age.at.Diagnosis","AGE",colnames(clin_METABRIC))
colnames(clin_METABRIC) <- gsub("Pam50...Claudin.low.subtype","PAM50_and_CLAUDIN_LOW",colnames(clin_METABRIC))

# Add PAM50 calls
df_metabric_PAM50 <- read.delim(file_metabric_PAM50, stringsAsFactors = F)
colnames(df_metabric_PAM50)[1] <- "SAMPLE_ID"
clin_METABRIC <- dplyr::left_join(clin_METABRIC, df_metabric_PAM50, by = c("Sample.ID" = "SAMPLE_ID"))
colnames(clin_METABRIC) <- toupper(colnames(clin_METABRIC))

# add to clinical info
clin_METABRIC <- dplyr::inner_join(clin_METABRIC, df_ecotypes, by = c("SAMPLE.ID" = "SAMPLE_ID"))
nrow(clin_METABRIC)
rownames(clin_METABRIC) <- clin_METABRIC$SAMPLE_ID

# summarise number of tumours in each ecotype
df_ecotype_tumour_counts <- as.data.frame(table(clin_METABRIC$Ecotype))
colnames(df_ecotype_tumour_counts)[1] <- "Ecotype"
colnames(df_ecotype_tumour_counts)[2] <- "number_of_tumours"
fOut <- paste0(output_directory, figNum, ".ecotype_tumour_counts.txt")
write.table(df_ecotype_tumour_counts, file = fOut, sep = "\t", row.names = FALSE)

gg <- run_survival_analysis_logrank(clinical_data = clin_METABRIC, risk.table = F, conf.int = F, end_time = end_time, palette = as.character(dfEcotypeColours$value), exclude_samples = FALSE)

fPlot <- paste0(output_directory, figNum, ".pdf")
pdf(file = fPlot, width = 9, height = 6)
print(gg)
dev.off()

# ============================
figNum <- "8E"
message(paste0("FIGURE ", figNum))
# ============================
message("KM plot of Ecotype E2 vs E7")

clinical_data <- clin_METABRIC %>% filter(Ecotype %in% c("E2", "E7"))
ecotype_colours <- as.character(dfEcotypeColours$value)[c(2, 7)]
gg <- run_survival_analysis_logrank(clinical_data = clinical_data, risk.table = F, conf.int = F, end_time = end_time, palette = ecotype_colours, exclude_samples = FALSE)

fPlot <- paste0(output_directory, figNum, ".pdf")
pdf(file = fPlot, width = 9, height = 6)
print(gg)
dev.off()

# ============================
figNum <- "8F"
message(paste0("FIGURE ", figNum))
# ============================
message("KM plot of Ecotype E4 vs E7")

clinical_data <- clin_METABRIC %>% filter(Ecotype %in% c("E4", "E7"))
ecotype_colours <- as.character(dfEcotypeColours$value)[c(4, 7)]
gg <- run_survival_analysis_logrank(clinical_data = clinical_data, risk.table = F, conf.int = F, end_time = end_time, palette = ecotype_colours, exclude_samples = F)

fPlot <- paste0(output_directory, figNum, ".pdf")
pdf(file = fPlot, width = 9, height = 6)
print(gg)
dev.off()

# ============================
figNum <- "8G"
message(paste0("FIGURE ", figNum))
# ============================
message("This is the cell-type annotation figure (External file)")

# ============================
message("save image")
# ============================
fOut <- paste0(output_directory, "figure_8.RData")
save.image(file = fOut)

# ============================
message("SESSION INFO")
# ============================
sessionInfo()
# BRCA CELL ATLAS STUDY PLOTTING CODE:
# FIGURE 5 - STROMAL CELLS
# SUNNY Z WU
#
#
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_seuratdev
# R
#
#
# 01: SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(monocle)
library(reshape2)
library(pheatmap)

# DIRECTORIES 
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/05_figures/")
dir.create("05_STROMAL_FIGURES_v6")
setwd("05_STROMAL_FIGURES_v6")

# 02: LOAD ANNOTATED OBJECT ---------------------------------------------

# load all mesenchymal cells 
seurat_10X_integrated <- readRDS("../05_STROMAL_FIGURES_v3/Rdata/Stromal_ANNOTATED_object.Rdata")

# load individual cells 
for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  
  temp_subset <- readRDS(paste0("../05_STROMAL_FIGURES_v3/Rdata/seurat_object_subset_",cluster,".Rdata"))
  n <- paste0(paste0("temp_subset_",cluster))
  assign(n,temp_subset)
  
  temp_subset <- readRDS(paste0("../05_STROMAL_FIGURES_v3/Rdata/seurat_object_subset_",cluster,"_withnormals.Rdata"))
  n <- paste0(paste0("temp_subset_withnormal",cluster))
  assign(n,temp_subset)
  
  HSMM <- readRDS(paste0("../05_STROMAL_FIGURES_v3/Rdata/monocle_object_HSMM_",cluster,".Rdata"))
  n <- paste0(paste0("HSMM_",cluster))
  assign(n,HSMM)
  
  HSMM <- readRDS(paste0("../05_STROMAL_FIGURES_v3/Rdata/monocle_object_HSMM_",cluster,"_withnormals.Rdata"))
  n <- paste0(paste0("HSMM_withnormal",cluster))
  assign(n,HSMM)
  
}

# 03: UMAP ALL STROMAL CELLS ----------------------------------------------------

dir.create("01_ALL_STROMAL_CELLS")
reduction <- "UMAP"

temp_colname <- paste0("celltype_subset")
temp_reduction <- paste0(reduction,"COMBINEDFILTERED",30)

# function dimensions
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 6,
        height = 6,
        res = 450,
        units = 'in'
      )
    }
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 6,
        height = 6,
        useDingbats = F
      )
    }


    temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      label.size = 6,
      pt.size = 0.1,
      label = T,
      repel=T,
      reduction = temp_reduction,
      group.by = temp_colname
    )

    # # get legend
    library(ggpubr)
    leg <- get_legend(temp_dimplot)
    leg <- as_ggplot(leg)

    # remove legend
    temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
      xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))

    temp_pdf_function(paste0("01_ALL_STROMAL_CELLS/",reduction,"_02_dimplot.pdf"))
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


# 04: TSNE INDIVIDUAL STROMAL CELL TYPES ------------------------------------------

# extended data files
reduction <- "TSNE"
res <- 0.8
npcs <- 20

# plotting  dimension functions
temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 5,
      height = 5,
      res = 450,
      units = 'in'
    )
  }
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 5,
      height = 5,
      useDingbats = F
    )
  }


for(cluster in c("CAFs", "PVL", "Endo")){

  temp_subset <- get(paste0("temp_subset_",cluster))
  print(dim(temp_subset))

  temp_dimplot <- DimPlot(
    object = temp_subset,
    label.size = 8,
    pt.size = 0.25,
    reduction = paste0(reduction,"FILTEREDSUBSET",npcs),
    label = T,
    group.by = paste0("FILTEREDSUBSET", npcs, "_res.",res)
  )

  # get legend
  leg <- get_legend(temp_dimplot)
  leg <- as_ggplot(leg)

  # remove legend and change axes labels
  temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
    xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))

  temp_pdf_function(paste0("01_ALL_STROMAL_CELLS/04_dimplot_", cluster, "_", reduction,"_", res,".pdf"))
  print(temp_dimplot)
  dev.off()

  # by State
  temp_dimplot <- DimPlot(
    object = temp_subset,
    label.size = 8,
    pt.size = 0.25,
    reduction = paste0(reduction,"FILTEREDSUBSET",npcs),
    label = T,
    group.by ="State"
  )

  # get legend
  leg <- get_legend(temp_dimplot)
  leg <- as_ggplot(leg)

  # remove legend and change axes labels
  temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
    xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))

  temp_pdf_function(paste0("01_ALL_STROMAL_CELLS/04_dimplot_", cluster, "_", reduction,"_", "State",".pdf"))
  print(temp_dimplot)
  dev.off()

  #celltype annotations
  temp_dimplot <- DimPlot(
    object = temp_subset,
    label.size = 4,
    pt.size = 0.25,
    reduction = paste0(reduction,"FILTEREDSUBSET",npcs),
    label = T,
    group.by ="celltype_subset"
  )

  # remove legend and change axes labels
  temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
    xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))

  temp_pdf_function(paste0("01_ALL_STROMAL_CELLS/05_dimplot_", cluster, "_", reduction,"_", "celltype_subset",".pdf"))
  print(temp_dimplot)
  dev.off()

}


# 05: FEATUREPLOT - MARKERS OF INTEREST --------------------------------------------------

dir.create("02_markers_of_interest")

# plotting dimension functions
temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 6,
      height = 8,
      res = 300,
      units = 'in'
    )
  }
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 6,
      height = 8,
      useDingbats = F
    )
  }

# markers to plot
temp_marker_subset <- c("PDGFRA",
                        "COL1A1",
                        "ACTA2",
                        "PDGFRB",
                        "MCAM",
                        "PECAM1",
                        "CD34",
                        "VWF")

DefaultAssay(seurat_10X_integrated) <- "RNA"
reduction <- "UMAP"
temp_reduction <- paste0(reduction,"COMBINEDFILTERED",30)

for(gene in temp_marker_subset) {
  temp_featureplot <- FeaturePlot(
    object = seurat_10X_integrated,
    features = gene,
    pt.size = 0.01,
    reduction = temp_reduction,
    order = T)
  temp_featureplot <-
    temp_featureplot + theme(plot.title = element_text(size=24,
                                                       face = "italic"),
                             plot.margin = unit(c(0,0,0,0), "cm"),
                             axis.line=element_blank(),
                             axis.text.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks=element_blank(),
                             axis.title.x=element_blank(),
                             axis.title.y=element_blank(),
                             legend.key.size = unit(0.5, 'cm'), #change legend key size
                             legend.key.height = unit(0.5, 'cm'), #change legend key height
                             legend.key.width = unit(0.5, 'cm'), #change legend key width
                             legend.text = element_text(size=18) #change legend text font size
    ) +
    xlab(" ") + ylab(" ")
  
  
  
  n <- paste0("temp_featureplot_",gene)
  assign(n,temp_featureplot)
}

temp_grid <-
  plot_grid(plotlist=mget(paste0("temp_featureplot_",temp_marker_subset)),
            nrow = 4,
            ncol = 2)

temp_pdf_function(paste0("02_markers_of_interest/02_markers_sigs_to_plot.pdf"))
print(temp_grid)
dev.off()


temp_featureplot <- FeaturePlot(
  object = seurat_10X_integrated,
  features = temp_marker_subset,
  pt.size = 0.01,
  reduction = temp_reduction,
  order = T, ncol = 2, coord.fixed = T)
temp_pdf_function(paste0("02_markers_of_interest/02_markers_sigs_to_plot_v2.pdf"))
print(temp_featureplot)
dev.off()

# 06: MONOCLE PSEUDOTIME PLOTS -----------------------------------------------------------

dir.create("03_MONOCLE")

# plot 
for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  HSMM <- get(paste0("HSMM_",cluster))
  colnames(pData(HSMM))[93] <- "Cell_state"
  ## order cells change colors and theta to match your plot
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
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 6,
        height = 6,
        useDingbats = F
      )
    }
  
  if(cluster == "CAFs"){
    temp_row <- 2
  } else {
    temp_row <- 1
  }
  
  type <- "Cell_state"
  temp_ggplot <- plot_cell_trajectory(HSMM, 
                                      color_by = type,
                                      show_branch_points = T,
                                      show_tree = TRUE,
                                      cell_size = 0.5) 
  temp_ggplot <-
    temp_ggplot + theme(legend.text = element_text(size=18),
                        legend.title = element_text(size=18),
                        legend.position = "bottom"
    ) + guides(colour=guide_legend(nrow=temp_row,
                                   byrow=TRUE,
                                   override.aes = list(size=5)))
  
  # temp_ggplot <- temp_ggplot + theme(
  #   legend.title = element_blank(),
  #   legend.position = "none")
  
  
  temp_pdf_function(paste0("03_MONOCLE/",cluster,"_01_plot_cell_trajectory_",type,".pdf"))
  print(temp_ggplot)
  dev.off()
  
  # by pseudotime
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 6,
        height = 5,
        res = 300,
        units = 'in'
      )
    }
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 6,
        height = 5,
        useDingbats = F
      )
    }
  
  
  type <- "Pseudotime"
  temp_ggplot <- plot_cell_trajectory(HSMM, 
                                      color_by = type,
                                      show_branch_points = T,
                                      show_tree = TRUE,
                                      cell_size = 0.5) 
  temp_ggplot <-
    temp_ggplot + theme(legend.text = element_text(size=10),
                        legend.title = element_text(size=18),
                        legend.position = "bottom"
    ) 
  
  
  temp_pdf_function(paste0("03_MONOCLE/",cluster,"_01_plot_cell_trajectory_",type,".pdf"))
  print(temp_ggplot)
  dev.off()
  
  if(cluster == "CAFs"){
    temp_diff_genes_to_plot <- c("ALDH1A1", "CXCL12", "ACTA2", "COL1A1")
    ncol <- 1
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 4,
          height = 8,
          useDingbats = F
        )
      }
  }
  if(cluster == "PVL"){
    temp_diff_genes_to_plot <- c("ALDH1A1", "THY1", "CD36", "MYH11")
    ncol <- 2
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 6,
          height = 4,
          useDingbats = F
        )
      }
  }
  if(cluster == "Endo"){
    temp_diff_genes_to_plot <- c("ACKR1", "RGS5", "DLL4", "CXCL12")
    ncol <- 2
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 6,
          height = 4,
          useDingbats = F
        )
      }
  }
  
  temp_diff_genes_to_plot_filtered <- 
    temp_diff_genes_to_plot[temp_diff_genes_to_plot %in% rownames(HSMM)]
  
  temp_png_function <-
    function(x) {
      png(
        file = (x),
        width = 4,
        height = 11,
        res = 300,
        units = 'in'
      )
    }
  
  temp_pdf_function(paste0("03_MONOCLE/",cluster,"_03_plot_genes_in_branched_pseudotime.pdf"))
  temp_plot <-   
    plot_genes_branched_pseudotime(HSMM[temp_diff_genes_to_plot_filtered,], 
                                   color_by = "Cell_state",
                                   cell_size=0.3,
                                   ncol=ncol)
  temp_plot <- 
    temp_plot + theme(legend.text = element_blank(),
                      legend.title = element_blank(),
                      legend.position = "none"
    )      
  print(temp_plot)
  dev.off()
  
}


# 07: CLUSTER AVERAGED HEATMAP BY STATE --------------------------------------------------

# extended data figures

    dir.create("04_HEATMAP_BY_STATES")
    temp_colname <- "celltype_subset"

      # markers to highlight
      sigGenes_v_CAFs <- c("LEPR", "KLF4", "ALDH1A1", "CXCL12", "C3", "C7", "CFD",
                           "DLK1", "APOE",
                           "COL4A1", "RGS5", "SFRP4",
                          "COL1A1", "FN1", "ACTA2", "TAGLN", "FAP"
                           )
      sigGenes_v_PVL <- c("RGS5", "CD36", "CD44",
                            "ITGA1", "ICAM1", "VCAM1",
                            "MYH11", "CD248", "ACTA2", "CSPG4")
      sigGenes_v_Endo <- c("ACKR1", "DLL4", "VEGFC", "SELE", "SELP",                            "RGS5", "ANGPT2", "FLT1", "IL6")


      # new factorisation
      temp_ordering_Endo <- c("ECs_s1", "ECs_s2","ECs_s3")
      temp_ordering_PVL <- c("PVL_s1", "PVL_s2","PVL_s3")
      temp_ordering_CAFs <- c("CAF_s1", "CAF_s2","CAF_s3",
                              "CAF_s4", "CAF_s5")


    for(cluster in c("CAFs", "PVL", "Endo")){

      sigGenes_v <- get(paste0("sigGenes_v_",cluster))
      temp_ordering <- get(paste0("temp_ordering_",cluster))

      print(cluster)
      temp_subset <- get(paste0("temp_subset_",cluster))
      temp_subset@meta.data[,temp_colname] <-
        factor(temp_subset@meta.data[,temp_colname],
               levels=temp_ordering)

      print(dim(temp_subset))

      DefaultAssay(temp_subset) <- "RNA"
      Idents(temp_subset) <- temp_colname

      temp_cluster.allmarkers <-
        read.csv(paste0("../05_STROMAL_FIGURES_v3/04A_DGE_by_monocle_states/01",cluster,"_MAST_FINDALLMARKERS_RNA_State.csv"))

      temp_cluster.averages <-
        AverageExpression(temp_subset,
                          return.seurat = TRUE,
                          verbose = T,
                          assays="RNA",
                          features=unique(c(as.vector(temp_cluster.allmarkers$gene),sigGenes_v)))

      temp_data_frame <-
        as.matrix(temp_cluster.averages@assays$RNA@data)

      write.csv(temp_data_frame,
                paste0("04_HEATMAP_BY_STATES/",cluster,"_cluster_avg_GE.csv"))

      assign(paste0("temp_data_frame_",cluster),
             temp_data_frame)

    }

    # plot
    library("viridis")
    hmcol <- inferno(24)

    # temp_png_function <-
    #   function(x) {
    #     png(
    #       file = (x),
    #       width = 12,
    #       height = 10,
    #       res = 300,
    #       units = 'in'
    #     )
    #   }

    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 12,
          height = 10,
          useDingbats=F
        )
      }

    for(legend in c(T,F)){
      for(cluster in c("CAFs", "PVL", "Endo")){

        print(cluster)
        # temp_data_frame_hclust <- read.csv(paste0("avg_GE_heatmap/",cluster,"_cluster_avg_GE.csv"),
        #                                    row.names = "X")
        temp_data_frame <- get(paste0("temp_data_frame_",cluster))

        temp_cluster.allmarkers <-
          read.csv(paste0("../05_STROMAL_FIGURES_v3/04A_DGE_by_monocle_states/01",cluster,"_MAST_FINDALLMARKERS_RNA_State.csv"))

        # temp_cluster.allmarkers <- temp_cluster.allmarkers[temp_cluster.allmarkers$avg_logFC > 0.5,]
        sigGenes_v <- get(paste0("sigGenes_v_",cluster))

        temp_cluster.allmarkers <-
          (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(40,
                                                                   -p_val))

        print(table(temp_cluster.allmarkers$cluster))
        print(length(unique(temp_cluster.allmarkers$gene)))

        temp_data_frame_hclust <- temp_data_frame[rownames(temp_data_frame) %in% unique(c(as.vector(temp_cluster.allmarkers$gene),sigGenes_v)),]


        # genes of interest to highlight
        # https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
        {
          rowMeta_df <- data.frame(Sig = rep("No", length(rownames(temp_data_frame_hclust))),
                                   stringsAsFactors = F,
                                   row.names = rownames(temp_data_frame_hclust))
          for (gene_v in sigGenes_v) rowMeta_df[rownames(rowMeta_df) == gene_v, "Sig"] <- gene_v

        }
        # temp_png_function(paste0("04_HEATMAP_BY_STATES/legend_01_pheatmap_",cluster,".png"))
        #
        # heat <- pheatmap(temp_data_frame_hclust,
        #                  color = rev(hmcol),
        #                  cluster_cols = F,
        #                  cluster_rows = T,
        #                  scale = "row",
        #                  clustering_distance_cols = "correlation",
        #                  clustering_distance_rows = "correlation",
        #                  # annotation_row = rowMeta_df,
        #                  fontsize_row = 5,
        #                  show_rownames = F,
        #                  show_colnames = T,
        #                  fontsize_col = 8,
        #                  # annotation_legend = T,
        #                  cellheight= 2,
        #                  cellwidth = 10,
        #                  gaps_col = NULL,
        #                  # annotation_names_col = T,
        #                  angle_col = 45,
        #                  # treeheight_row = 75,
        #                  # treeheight_col = 50,
        #                  legend = T,
        #                  border_color=FALSE,
        #                  main = paste0(cluster)
        # )
        # #
        # # print(heat)
        # dev.off()
        #
        # temp_pdf_function(paste0("avg_GE_heatmap_by_states/legend_01_pheatmap_",cluster,".pdf"))
        # print(heat)
        # dev.off()

        if(cluster == "PVL"){
          # choices: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
          clust_meth <- "maximum"
        } else clust_meth <- "correlation"

        heat <- pheatmap(temp_data_frame_hclust,
                         color = rev(hmcol),
                         cluster_cols = F,
                         cluster_rows = T,
                         scale = "row",
                         # clustering_distance_cols = clust_meth,
                         clustering_distance_rows = clust_meth,
                         # annotation_row = rowMeta_df,
                         fontsize_row = 3,
                         show_rownames = T,
                         show_colnames = T,
                         fontsize_col = 10,
                         # annotation_legend = T,
                         cellheight= 2.5,
                         cellwidth = 15,
                         gaps_col = NULL,
                         # annotation_names_col = T,
                         angle_col = 45,
                         treeheight_row = 15,
                         # treeheight_col = 50,
                         legend = legend,
                         border_color=FALSE,
                         main = paste0(" "),
                         silent = T
        )

        ### Get order of genes after clustering
        genesInHeatOrder_v <- heat$tree_row$labels[heat$tree_row$order]
        whichSigInHeatOrder_v <- which(genesInHeatOrder_v %in% sigGenes_v)

        # print("exists in data frame")
        # for(gene in sigGenes_v){
        #   print(grep(paste0(gene),
        #              genesInHeatOrder_v,
        #              value=T))
        # }

        whichSigInHeatOrderLabels_v <- genesInHeatOrder_v[whichSigInHeatOrder_v]


        sigY <- 1 - (1/length(rownames(temp_data_frame_hclust)) * whichSigInHeatOrder_v)

        ### Change title
        whichMainGrob_v <- which(heat$gtable$layout$name == "main")
        heat$gtable$grobs[[whichMainGrob_v]] <- textGrob(label = paste0(cluster),
                                                         gp = gpar(fontsize = 16))

        ### Remove rows
        whichRowGrob_v <- which(heat$gtable$layout$name == "row_names")
        heat$gtable$grobs[[whichRowGrob_v]] <- textGrob(label = whichSigInHeatOrderLabels_v,
                                                        y = sigY,
                                                        vjust = 1)

        #
        add.flag <- function(pheatmap,
                             kept.labels,
                             repel.degree) {

          # repel.degree = number within [0, 1], which controls how much
          #                space to allocate for repelling labels.
          ## repel.degree = 0: spread out labels over existing range of kept labels
          ## repel.degree = 1: spread out labels over the full y-axis

          heatmap <- pheatmap$gtable

          new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]

          # keep only labels in kept.labels, replace the rest with ""
          new.label$label <- ifelse(new.label$label %in% kept.labels,
                                    new.label$label, "")

          # calculate evenly spaced out y-axis positions
          repelled.y <- function(d, d.select, k = repel.degree){
            # d = vector of distances for labels
            # d.select = vector of T/F for which labels are significant

            # recursive function to get current label positions
            # (note the unit is "npc" for all components of each distance)
            strip.npc <- function(dd){
              if(!"unit.arithmetic" %in% class(dd)) {
                return(as.numeric(dd))
              }

              d1 <- strip.npc(dd$arg1)
              d2 <- strip.npc(dd$arg2)
              fn <- dd$fname
              return(lazyeval::lazy_eval(paste(d1, fn, d2)))
            }

            full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
            selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

            return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                            to = min(selected.range) - k*(min(selected.range) - min(full.range)),
                            length.out = sum(d.select)),
                        "npc"))
          }
          new.y.positions <- repelled.y(new.label$y,
                                        d.select = new.label$label != "")
          new.flag <- segmentsGrob(x0 = new.label$x - unit(0.5, "npc"),
                                   x1 = new.label$x + unit(0.25, "npc"),
                                   y0 = new.label$y[new.label$label != ""],
                                   y1 = new.y.positions)

          # shift position for selected labels
          new.label$x <- new.label$x + unit(1.25, "npc")
          new.label$y[new.label$label != ""] <- new.y.positions

          # add flag to heatmap
          heatmap <- gtable::gtable_add_grob(x = heatmap,
                                             grobs = new.flag,
                                             t = 4,
                                             l = 4
          )

          # replace label positions in heatmap
          heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

          # plot result
          grid.newpage()
          grid.draw(heatmap)

          # return a copy of the heatmap invisibly
          invisible(heatmap)
        }

        # temp_png_function(paste0("04_HEATMAP_BY_STATES/01_pheatmap_",cluster,".png"))
        # add.flag(heat,
        #          kept.labels = sigGenes_v,
        #          repel.degree = 1)
        # dev.off()


        temp_pdf_function(paste0("04_HEATMAP_BY_STATES/02_pheatmap_withlegend_",legend, "_",cluster,".pdf"))
        add.flag(heat,
                 kept.labels = sigGenes_v,
                 repel.degree = 1)
        dev.off()

        # print(sigGenes_v)
      }
    }

# 08: CITE-SEQ LOAD OBJECT  -----------------------------------------------------

dir.create("06_imputed_CITESeq_markers_all_stromal/")

# load cite-seq object for stromal cells
seurat_10X_CITEimputed_STROMAL_raw <- readRDS("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/07_CITE_imputation_Nenad/run03_191203/output/imputed_STROMAL_with_lymphatics_raw.Rdata")
seurat_10X_CITEimputed_STROMAL_raw <- subset(seurat_10X_CITEimputed_STROMAL_raw,
                                             cells = row.names(seurat_10X_integrated@meta.data))
print(dim(seurat_10X_CITEimputed_STROMAL_raw))

# append to individual stromal objects
for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  temp_subset <- get(paste0("temp_subset_",cluster))
  
  for(norm in c("raw")){
    print(norm)
    seurat_10X_integrated_temp <- seurat_10X_CITEimputed_STROMAL_raw
    
    seurat_10X_integrated_subset <- subset(seurat_10X_integrated_temp,
                                           cells = row.names(temp_subset@meta.data))
    
    temp_df <- data.frame(row.names = row.names(temp_subset@meta.data),
                          State = temp_subset@meta.data$State,
                          celltype_subset = temp_subset@meta.data$celltype_subset,
                          celltype_major = temp_subset@meta.data$celltype_major,
                          celltype_minor = temp_subset@meta.data$celltype_minor
    )
    temp_df <- temp_df[row.names(seurat_10X_integrated_subset@meta.data),,drop=F]
    
    if(all.equal(rownames(temp_df), rownames(seurat_10X_integrated_subset@meta.data))){
      seurat_10X_integrated_subset@meta.data$State <- temp_df$State
      seurat_10X_integrated_subset@meta.data$celltype_subset <- temp_df$celltype_subset
      seurat_10X_integrated_subset@meta.data$celltype_major <- temp_df$celltype_major
      seurat_10X_integrated_subset@meta.data$celltype_minor <- temp_df$celltype_minor
    }
    
    m <- paste0("seurat_10X_CITEimputed_",cluster, "_",norm)
    assign(m, seurat_10X_integrated_subset)
    
  }
}

# 09: CITE-SEQ UMAPS -----------------------------------------------------

# plot markers of interest
temp_markers_of_interest_STROMAL <- c("adt_CITE-Podoplanin", 
                                      "adt_CITE-MCAM", 
                                      "adt_CITE-CD31",
                                      "adt_CITE-CD34")

# Stromal Featureplots
reduction <- "UMAP"
temp_reduction <- paste0(reduction,"COMBINEDFILTERED",30)
cluster <- "STROMAL"
norm <- "raw"
# plotting dimensions
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 12,
      height = 3,
      useDingbats = F
    )
  }

temp_markers_of_interest <- get(paste0("temp_markers_of_interest_",cluster))

seurat_10X_integrated_temp <- 
      get(paste0("seurat_10X_CITEimputed_",cluster, "_",norm))
    
for(gene in temp_markers_of_interest) {
      temp_featureplot <- FeaturePlot(
        object = seurat_10X_integrated_temp,
        features = gene,
        pt.size = 0.1,
        reduction = temp_reduction,
        order = T)
      temp_featureplot <- 
        temp_featureplot + theme(plot.title = element_text(size=18,
                                                           face = "italic"),
                                 plot.margin = unit(c(0,0,0,0), "cm")
        ) +
        xlab("UMAP_2") + ylab("UMAP_1")
      
      
      
      n <- paste0("temp_featureplot_",gene)
      assign(n,temp_featureplot)
    }
    
    temp_grid <-
      plot_grid(plotlist=mget(paste0("temp_featureplot_",temp_markers_of_interest)),
                nrow = 1,
                ncol = 4)
    
    temp_pdf_function(paste0("06_imputed_CITESeq_markers_all_stromal/01_FILTERED_imputed_CITE_v2_",
                             cluster, "_markers_of_interest_", "_", norm, 
                             ".pdf"))
    print(temp_grid)
    dev.off()


# 09: CITE-SEQ HEATMAPS -------------------------------------------------------

dir.create("CITE_seq_heatmaps/")
temp_markers_of_interest <- c("adt_CITE-Podoplanin","adt_CITE-B7H3", "adt_CITE-B7-H4", 
                              "adt_CITE-CD10", "adt_CITE-CD34", "adt_CITE-CD40",
                              "adt_CITE-MCAM", "adt_CITE-Thy-1", 
                              "adt_CITE-CD49a",  "adt_CITE-CD49d",
                              "adt_CITE-CD31","adt_CITE-CD34",
                              "adt_CITE-CD49f", "adt_CITE-CD49b",
                              "adt_CITE-MHCII","adt_CITE-CD141", 
                              "adt_CITE-CD73","adt_CITE-CD40"
)
temp_markers_of_interest <- gsub("adt_", "", temp_markers_of_interest)

temp_CITEmerged <- merge(seurat_10X_CITEimputed_CAFs_raw,
                         seurat_10X_CITEimputed_PVL_raw)
temp_CITEmerged <- merge(temp_CITEmerged,
                         seurat_10X_CITEimputed_Endo_raw)
temp_data_frame <- 
  as.matrix(temp_CITEmerged@assays$ADT@data)
temp_data_frame <- as.data.frame(t(temp_data_frame))
temp_data_frame <- temp_data_frame[, colnames(temp_data_frame) %in% temp_markers_of_interest,drop=F]

print(all.equal(rownames(temp_CITEmerged@meta.data),
                row.names(temp_data_frame)))

temp_data_frame$cluster <- temp_CITEmerged@meta.data$celltype_subset
# temp_data_frame$barcode <- rownames(temp_data_frame)

# plot
temp_data_frame_combined.m <- melt(temp_data_frame)

# cluster average
temp_data_frame_combined.m.agg <-
  aggregate(.~cluster+variable,
            temp_data_frame_combined.m,
            mean)
temp_data_frame_combined.m.agg <- arrange(temp_data_frame_combined.m.agg,
                                          cluster)
temp_data_frame_combined.m.agg.dcast <-
  dcast(data = temp_data_frame_combined.m.agg,
        formula = variable~cluster,
        fun.aggregate = sum,
        value.var = "value")
#row names to variables
rownames(temp_data_frame_combined.m.agg.dcast) <- 
  temp_data_frame_combined.m.agg.dcast$variable
temp_data_frame_combined.m.agg.dcast <- 
  temp_data_frame_combined.m.agg.dcast[, ! colnames(temp_data_frame_combined.m.agg.dcast) %in% "variable"]
temp_data_frame_hclust <- as.matrix(temp_data_frame_combined.m.agg.dcast)

# plot
## colour scale
library("viridis")  
hmcol <- viridis(24)
heat <- pheatmap(temp_data_frame_hclust, 
                 filename = "CITE_seq_heatmaps/01_pheatmap.pdf",
                 color = rev(hmcol),
                 cluster_cols = F, 
                 cluster_rows = T, 
                 scale = "row", 
                 fontsize_row = 10,
                 show_colnames = T, 
                 fontsize_col = 15,
                 annotation_legend = T,
                 cellheight=10, 
                 cellwidth =20,
                 gaps_col = NULL, 
                 annotation_names_col = T, 
                 angle_col = 45,
                 treeheight_row = 50, 
                 legend = T,
                 border_color=FALSE)

# 10: GENE ONT PATHWAY FIGURE ------------------------------------------------

dir.create("functional_enrichment/")
library(data.table)

# combined DF
temp_df <- NULL
temp_annotation_df <- NULL
for(cluster in c("CAFs", "PVL", "Endo")){
  
  temp_ont_unique <- NULL
  for(ont in c("enrichKEGG", "enrichGO")) {
    print(ont)
    # ont <- "enrichKEGG"
    type <- "Description"
    temp_go_df <- read.csv(paste0("../05_STROMAL_FIGURES_v3/functional_enrichment/",cluster,"/01_compareCluster_",cluster,"_",ont,".csv"))
    
    temp_go_df_sorted <- 
      temp_go_df[ order( temp_go_df[, "clusterID"],
                         temp_go_df[, type] ),]
    
    if(! ont == "groupGO"){
      temp_go_df_sorted$logp_val <- (-(log10(temp_go_df_sorted$p.adjust)))
    } else {
      temp_go_df_sorted$logp_val <- temp_go_df_sorted$Count
    }
    
    temp_go_df_sorted <- temp_go_df_sorted[ , c("logp_val",
                                                "clusterID",
                                                type)]
    
    rownames(temp_go_df_sorted) <- NULL
    temp_go_combined_sorted_temp <- temp_go_df_sorted
    
    colnames(temp_go_df_sorted) <- c("logp_val", "clusterID", "type")
    
    temp_go_df_sorted$lineage <- cluster
    temp_go_df_sorted$clusterID2 <- paste0(temp_go_df_sorted$clusterID) 
    
    temp_df <- rbind(temp_df,temp_go_df_sorted)
    temp_ont_unique <- rbind(temp_ont_unique,temp_go_df_sorted)
    
    print(unique(temp_go_df_sorted$clusterID))
  }
  
  temp_annotation_df_subset <- data.frame(row.names = as.vector(unique(temp_ont_unique$clusterID2)))
  
  temp_annotation_df_subset$lineage <- cluster
  temp_annotation_df <- rbind(temp_annotation_df,
                              temp_annotation_df_subset)
}

# plot combined DF
temp_go_df_sorted <- temp_df
temp_top_pathways_to_plot <- 20

# plot
temp_combind_summary_df_top <-
  (temp_go_df_sorted %>% group_by(clusterID) %>% top_n(temp_top_pathways_to_plot,logp_val))

# assign value to every pathway for every cluster (0s)
temp_combind_summary_df_dcast <-
  dcast(data = temp_combind_summary_df_top,
        formula = type~clusterID2,
        fun.aggregate = sum,
        value.var = "logp_val")

#row names to variables
rownames(temp_combind_summary_df_dcast) <- 
  temp_combind_summary_df_dcast$type

temp_combind_summary_df_dcast <- 
  temp_combind_summary_df_dcast[, ! colnames(temp_combind_summary_df_dcast) %in% "type"]

temp_data_frame_hclust <- as.matrix(temp_combind_summary_df_dcast)

# pheatmap
library("viridis")  
hmcol <- viridis(24)

# plot
pheatmap(temp_data_frame_hclust, filename = "functional_enrichment/01_pheatmap_combined.pdf",
         color = rev(hmcol),
         cluster_cols = F, 
         cluster_rows = T, 
         scale = "row", 
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         annotation_col = temp_annotation_df,
         fontsize_row = 6,
         show_rownames = T,
         show_colnames = T, 
         fontsize_col = 8,
         annotation_legend = T,
         cellheight= 6, 
         cellwidth = 15,
         gaps_col = NULL, 
         annotation_names_col = T,
         angle_col = 45,
         treeheight_row = 75,
         treeheight_col = 50,
         legend = T,
         main = paste0(" ")
)

# 11: ANTIGEN-PRESENTING CAF GENES ---------------------------------------------------------

dir.create(paste0("02_AUCELL_SIGNATURES/","apCAFs"))

# GOIs 
DefaultAssay(seurat_10X_integrated) <- "RNA"
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 12,
      height = 5,
      useDingbats=F
    )
  }

temp_pdf_function(paste0("02_AUCELL_SIGNATURES/apCAFs/05_GOIS_by_subset.pdf"))
temp_vlnplot <- VlnPlot(object = seurat_10X_integrated,
                        features = c("CLU", "CD74", "CAV1"),
                        pt.size = 0,
                        group.by = "celltype_subset",
                        ncol = 3)
print(temp_vlnplot)
dev.off()

# 12: NORMAL TISSUE COMPARISON --------------------------------------------

dir.create("SUPP_subtype_normal_comparisons")
reduction <- "TSNE"
res <- 0.8
npcs <- 20

# plotting dimensions
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 4, 
      height = 4, 
      res = 450, 
      units = 'in'
    )
  }
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 4,
      height = 4,
      useDingbats=F
    )
  }

# TSNE plots
for(cluster in c("CAFs", "PVL", "Endo")){
  
  temp_subset <- get(paste0("temp_subset_withnormal",cluster))
  
  # tumour vs normal
  temp_subset@meta.data$tumour_normal <- temp_subset@meta.data$subtype
  temp_subset@meta.data$tumour_normal[temp_subset@meta.data$tumour_normal == "TNBC"] <- "Tumour"
  temp_subset@meta.data$tumour_normal[temp_subset@meta.data$tumour_normal == "ER+"] <- "Tumour"
  temp_subset@meta.data$tumour_normal[temp_subset@meta.data$tumour_normal == "HER2+"] <- "Tumour"
  
  print(unique(temp_subset@meta.data$tumour_normal))
  
  for(type in c("tumour_normal", "subtype")){
    temp_dimplot <- DimPlot(
      object = temp_subset,
      label.size = 8,
      pt.size = 0.1,
      reduction = paste0(reduction,"FILTEREDSUBSET",npcs),
      label = F,
      group.by = type
    )
    
    # remove legend
    temp_dimplot <- temp_dimplot + 
      xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))
    
    # remove legend and change axes labels
    temp_dimplot <- temp_dimplot +
      xlab(paste0(reduction,"_2")) + ylab(paste0(reduction,"_1"))
    
    temp_pdf_function(paste0("SUPP_subtype_normal_comparisons/01_dimplot_", cluster,"_",type,".pdf"))
    print(temp_dimplot)
    dev.off()
  }
  
}

# monocle plots
for(cluster in c("CAFs", "PVL", "Endo")){
  print(cluster)
  HSMM <- get(paste0("HSMM_withnormal",cluster))
  
  pData(HSMM)$tumour_normal <- pData(HSMM)$subtype
  pData(HSMM)$tumour_normal[pData(HSMM)$tumour_normal == "TNBC"] <- "Tumour"
  pData(HSMM)$tumour_normal[pData(HSMM)$tumour_normal == "ER+"] <- "Tumour"
  pData(HSMM)$tumour_normal[pData(HSMM)$tumour_normal == "HER2+"] <- "Tumour"
  print(unique(pData(HSMM)$tumour_normal))
  
  pData(HSMM)$clusterID <- 
    paste0(pData(HSMM)$celltype_major,"_c",pData(HSMM)$FILTEREDSUBSET20_res.0.8)
  
  pData(HSMM)$clusterID <- factor(pData(HSMM)$clusterID,
                                  levels=str_sort(unique(pData(HSMM)$clusterID),numeric = T))
  
  
  for(type in c("tumour_normal", "subtype", "clusterID")){
    if(type=="clusterID"){
      temp_size <- 10
      temp_row <- 3
      temp_legend_size <- 3
      temp_png_function <-
        function(x) {
          png(
            file = (x), 
            width = 7, 
            height = 6, 
            res = 450, 
            units = 'in'
          )
        }
      temp_pdf_function <-
        function(x) {
          pdf(
            file = (x),
            width = 7,
            height = 6,
            useDingbats=F
          )
        }
      
      
    } else {
      temp_size <- 15
      temp_row <- 1
      temp_legend_size <- 5
      temp_png_function <-
        function(x) {
          png(
            file = (x), 
            width = 7, 
            height = 5, 
            res = 450, 
            units = 'in'
          )
        }
      temp_pdf_function <-
        function(x) {
          pdf(
            file = (x),
            width = 7,
            height = 5,
            useDingbats=F
          )
        }
      
    }
    
    temp_ggplot <- plot_cell_trajectory(HSMM, 
                                        color_by = type,
                                        show_branch_points = T,
                                        show_tree = TRUE,
                                        cell_size = 0.25) 
    temp_ggplot <-
      temp_ggplot + theme(legend.text = element_text(size=temp_size),
                          legend.title = element_text(size=15),
                          legend.position = "bottom"
      ) + guides(colour=guide_legend(nrow=temp_row,
                                     byrow=TRUE,
                                     override.aes = list(size=temp_legend_size)))
    
    
    temp_pdf_function(paste0("SUPP_subtype_normal_comparisons/02_",cluster,"_plot_cell_trajectory_",type,".pdf"))
    print(temp_ggplot)
    dev.off()
  }
}





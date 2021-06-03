# BRCA CELL ATLAS STUDY PLOTTING CODE:
# FIGURE 1 - OVERVIEW
# SUNNY Z WU
#
#
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_seuratdev
# R
#
# 01: SETUP -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(easyGgplot2)

# Colour blind friendly pallette
# BiocManager::install("colorBlindness")
library(colorBlindness)

# directory
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/05_figures/")
dir.create("01_OVERVIEW_v6_colourblind_friendly")
setwd("01_OVERVIEW_v6_colourblind_friendly")

# 02: LOAD OBJECTS --------------------------------------------------

# brca cell atlas object
seurat_10X <- readRDS("/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/Jul2020_updated_objectfile/Rdata/RDATA_06_BRCA_MINIATLAS_JUL2020_withumap.Rdata")

# MERGED - only stromal and immune cells (for supp figure)
seurat_10X_subset <- readRDS(paste0("/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/analysis/Jul2020_updated_objectfile/Rdata/RDATA_07_stromal_immune_Jul2020_MERGE_analysis.Rdata"))


# 03: COLOUR PALLETE ----------------------------------------------------------

temp_colour_pal_2 <- data.frame(celltype = c("T-cells", "B-cells", "Plasmablasts", "Myeloid", "Epithelial", "Cycling", "Mesenchymal", "Endothelial"),
                                colour = c("#56B4E9", "#F0E442", "#999999", "#009E73",
                                           "#E69F00", "#0072B2", "#D55E00", "#CC79A7"))

temp_colour_pal_2_stromal <- temp_colour_pal_2[!temp_colour_pal_2$celltype %in% c("Epithelial"),,drop=F]

# FACTORISE
seurat_10X@meta.data$celltype_major_fig <- factor(seurat_10X@meta.data$celltype_major_fig,levels = unique(temp_colour_pal_2$celltype))

seurat_10X_subset@meta.data$celltype_major_fig <- factor(seurat_10X_subset@meta.data$celltype_major_fig,levels = unique(temp_colour_pal_2_stromal$celltype))

# FIGURE 1A): UMAP  -----------------------------------------------------------------
  
  temp_colname <- "celltype_major_fig"

  # dimensions
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 8, 
        height = 7, 
        res = 300, 
        units = 'in'
      )
    }
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 8,
        height = 7,
        useDingbats = F
      )
    }
  numpcs <- 20 
  temp_reduction <- paste0("UMAPSIG",numpcs)
    
  # plot UMAP
    temp_dimplot <- DimPlot(
      object = seurat_10X,
      label.size = 8,
      pt.size = 0.01,
      label = F,
      reduction = temp_reduction,
      repel = T,
      group.by = temp_colname,
      cols = as.vector(temp_colour_pal_2$colour)
    )
    
    # get legend
    library(ggpubr)
    leg <- get_legend(temp_dimplot)
    leg <- as_ggplot(leg)
    
    
    # remove legend
    temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
      xlab("UMAP_2") + ylab("UMAP_1")
    
    temp_pdf_function(paste0("01_UMAP_PCs_",numpcs,".pdf"))
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
    temp_pdf_function(paste0("02_UMAP_legend_PCs",numpcs,".pdf"))
    print(leg)
    dev.off()

# FIGURE 1B): MARKER GENES ------------------------------------------------------------

  temp_genes_to_plot <-  c("EPCAM",
                           "MKI67",
                           "CD3D", 
                           "CD68",
                           "MS4A1",
                           "JCHAIN",
                           "PECAM1", 
                           "PDGFRB") 
  temp_reduction <- paste0("UMAPSIG",numpcs)
  for(gene in temp_genes_to_plot) {
        temp_featureplot <- FeaturePlot(
          object = seurat_10X,
          features = gene,
          pt.size = 0.001,
          reduction = temp_reduction,
          order = T)
        temp_featureplot <- 
          temp_featureplot + theme(plot.title = element_text(size=40,
                                                             face="italic"), 
                                   axis.text = element_text(size=25),
                                   axis.title = element_text(size=25),
                                   legend.key.size = unit(1.5, 'cm'), #change legend key size
                                   legend.key.height = unit(1.5, 'cm'), #change legend key height
                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                   legend.title = element_text(size=25), #change legend title font size
                                   legend.text = element_text(size=25) #change legend text font size
          ) +
          xlab("UMAP_2") + ylab("UMAP_1")
        
        
        n <- paste0("temp_featureplot_",gene)
        assign(n,temp_featureplot)
      }
    
    temp_grid <-
      plot_grid(plotlist=mget(paste0("temp_featureplot_",temp_genes_to_plot)),
                nrow = 2,
                ncol = 4)
    
    temp_png_function <-
      function(x) {
        png(
          file = (x),
          width = 24,
          height = 10,
          res = 900,
          units = 'in'
        )
      }
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 24,
          height = 10,
          useDingbats = F
        )
      }
    
    temp_pdf_function(paste0("02_FeaturePlot_as_grid_RNA_npcs_",numpcs, ".pdf"))
    print(temp_grid)
    dev.off()
    
    temp_png_function(paste0("02_FeaturePlot_as_grid_RNA_npcs_",numpcs, ".png"))
    print(temp_grid)
    dev.off()

# FIGURE 1C): PROPORTIONS PER PATIENT ACROSS SUBTYPES ---------------------------------------------

  temp_colname <- "celltype_major_fig"

  temp_cluster_names <- unique(seurat_10X@meta.data[,temp_colname])
  print(temp_cluster_names)
  # Clinical_IHC_colours <- c("ER+" = "blue", "TNBC" = "red", "HER2+" = "pink")

  temp_cellprop_df_all <- NULL
  for(i in c(1:length(temp_Names))) {
    temp_filtered_summary <-
      subset(seurat_10X@meta.data,
             orig.ident == temp_Names[i])

    temp_cellprop_df <-
      data.frame(unclass(table(temp_filtered_summary[,temp_colname])))

    temp_cellprop_df_all <- rbind(temp_cellprop_df_all,
                                  (temp_cellprop_df <-
                                     data.frame(sample = rep(temp_Names[i],
                                                             times = length(row.names(temp_cellprop_df))),
                                                cell_type = row.names(temp_cellprop_df),
                                                value = temp_cellprop_df[,1],
                                                proportions = temp_cellprop_df[,1]/colSums(temp_cellprop_df),
                                                subtype = unique(temp_filtered_summary$subtype)
                                     )))
  }


  # plot
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 12,
        height = 5,
        useDingbats=F
      )
    }

  # reorder celltypes
  temp_cellprop_df_all$cell_type <- factor(temp_cellprop_df_all$cell_type,
                                           levels=temp_colour_pal_2$celltype)

  temp_cellprop_df_all$subtype <- gsub("\\+","",temp_cellprop_df_all$subtype)
  temp_cellprop_df_all$subtype <- factor(temp_cellprop_df_all$subtype,
                                         levels=c("TNBC", "HER2", "ER"))

  temp_ggplot <- ggplot(temp_cellprop_df_all,
                        aes(x=sample,
                            fill=cell_type,
                            y=proportions)) +
    geom_bar(stat="identity") +
    theme(axis.text.x=element_text(angle=45,
                                   size=17.5, hjust=1),
          strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid"),
          strip.text.x = element_text(size=25, face = "bold"),
          axis.text = element_text(size=17.5), 
          axis.title = element_text(size=25, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(0.75, 'cm'), #change legend key size
          legend.key.height = unit(0.75, 'cm'), #change legend key height
          legend.key.width = unit(0.75, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=20) #change legend text font size
    ) + 
    # scale_y_continuovus(position = "right") +
    xlab("Patient ID") +
    ylab("Cell Proportions")  +
    guides(fill=guide_legend(title="Cell Type")) +
    # coord_flip()
    scale_fill_manual(values = as.vector(temp_colour_pal_2$colour)) +
    facet_grid(. ~ subtype,
               scales="free", space = "free") # free x removes white space

  temp_pdf_function(paste0("cell_type_proportions_by_subtype.pdf"))
  print(temp_ggplot)
  dev.off()




# ED FIGURE 1A-D): GENES AND UMIS PER TUMOUR --------------------------------------------------

  temp_df <- seurat_10X@meta.data[,colnames(seurat_10X@meta.data) %in% c("nFeature_RNA", "nCount_RNA", "subtype", "orig.ident", "celltype_major_fig"),drop=F]
  temp_df$barcode <- rownames(temp_df)
  rownames(temp_df) <- NULL

  temp_df$subtype <- gsub("\\+","",temp_df$subtype)
  temp_df$subtype <- factor(temp_df$subtype,levels=c("TNBC", "HER2", "ER"))

  temp_df$celltype_major_fig <- factor(temp_df$celltype_major_fig,
                                           levels=temp_colour_pal_2$celltype)

  temp_df$orig.ident <- factor(temp_df$orig.ident, levels=unique(temp_df$orig.ident))


for(temp_colname in c("orig.ident", "celltype_major_fig")){

      for(i in c("nFeature_RNA", "nCount_RNA")) {
        if(i == "nFeature_RNA"){
          temp_breaks <- c(0,500,1000,2500,5000,10000)
          temp_y <- "Number of Genes Per Cell"
        }
        if(i == "nCount_RNA"){
          temp_breaks <- c(0,500,2000,5000,10000,50000)
          temp_y <- "Number of UMIs Per Cell"
        }

        if(temp_colname == "orig.ident"){
          temp_pdf_function <-
            function(x) {
              pdf(
                file = (x),
                width = 8,
                height = 5,
                useDingbats=F
              )
            }

        temp_ggplot <- ggplot2.violinplot(data=temp_df,
                                          xName=temp_colname,
                                          yName=i,
                           groupName="subtype",
                           legendPosition="top",
                           faceting=TRUE,
                           addMean=TRUE,
                           meanPointShape=23, meanPointSize=1,
                           meanPointColor="black", meanPointFill="black",
                           facetingVarNames="subtype",
                           facetingScales="free",
                           facetingRect=list(background="white", lineType="solid",
                                             lineColor="black", lineSize=1.5),
                           facetingDirection="horizontal"
                           )

        temp_ggplot <- temp_ggplot + theme(axis.text.x = element_text(angle=45, hjust = 1),
                                           strip.text.x = element_text(size=15, face = "bold")) +
          facet_grid(. ~ subtype,
                     scales="free", space = "free") +
          scale_y_continuous(trans='log2', breaks = temp_breaks) +
          scale_fill_manual(values = c("red", "pink","blue")) +
          ylab(temp_y) + xlab("Patient ID")

        }

        if(temp_colname == "celltype_major_fig"){
          temp_pdf_function <-
            function(x) {
              pdf(
                file = (x),
                width = 8,
                height = 5,
                useDingbats=F
              )
            }

          temp_ggplot <- ggplot2.violinplot(data=temp_df,
                                            xName=temp_colname,
                                            yName=i,
                                            groupName=temp_colname,
                                            addMean=TRUE,
                                            meanPointShape=23, meanPointSize=3,
                                            meanPointColor="black", meanPointFill="black",
                                            legendPosition="none"
          )

          temp_ggplot <- temp_ggplot + theme(axis.text.x = element_text(angle=45, hjust = 1),
                                             strip.text.x = element_text(size=25, face = "bold")) +
            scale_y_continuous(trans='log2', breaks = temp_breaks) +
            scale_fill_manual(values = as.vector(temp_colour_pal_2$colour)) +
            ylab(temp_y) + xlab("Cell Type")

        }

        # temp_VlnPlot <- VlnPlot(seurat_10X,
        #                         features = i,
        #                         pt.size = 0,
        #                         group.by = temp_colname,
        #                         log = T)
        # temp_VlnPlot <- temp_VlnPlot + facet_grid(. ~ subtype,
        #                                           scales="free",
        #                                           space = "free")

        temp_pdf_function(paste0("suppfig1_metrics_",temp_colname,"_", i, ".pdf"))
        print(temp_ggplot)
        dev.off()
      }
}


# ED FIGURE 1E-F) STROMAL/IMMUNE CELLS MERGED: --------------------------------------
  
  temp_reduction <- paste0("UMAPMERGE",100)
  temp_colname <- "celltype_major_fig"
  
  # plot
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 8, 
        height = 7, 
        res = 300, 
        units = 'in'
      )
    }
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 8,
        height = 7,
        useDingbats = F
      )
    }
  
  # plot by celltype
  {      
    temp_dimplot <- DimPlot(
      object = seurat_10X_subset,
      label.size = 8,
      pt.size = 0.01,
      label = F,
      reduction = temp_reduction,
      repel = T,
      group.by = temp_colname,
      cols = as.vector(temp_colour_pal_2_stromal$colour)
    )
    
    # get legend
    library(ggpubr)
    leg <- get_legend(temp_dimplot)
    leg <- as_ggplot(leg)
    
    
    # remove legend
    temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
      xlab("UMAP_2") + ylab("UMAP_1")
    
    temp_pdf_function("03_UMAP_stromal_MERGE.pdf")
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
    temp_pdf_function("04_UMAP_stromal_MERGE_legend.pdf")
    print(leg)
    dev.off()
  }
  
  # plot by patient
  temp_pdf_function <-
    function(x) {
      pdf(
        file = (x),
        width = 8,
        height = 7,
        useDingbats = F
      )
    }
  {      
    temp_dimplot <- DimPlot(
      object = seurat_10X_subset,
      label.size = 8,
      pt.size = 0.01,
      label = F,
      reduction = temp_reduction,
      repel = T,
      group.by = "orig.ident",
      # cols = as.vector(temp_colour_pal_2_stromal$colour)
    )
    
    # get legend
    library(ggpubr)
    leg <- get_legend(temp_dimplot)
    leg <- as_ggplot(leg)
    
    
    # remove legend
    temp_dimplot <- temp_dimplot + theme(legend.position = "none") +
      xlab("UMAP_2") + ylab("UMAP_1")
    
    temp_pdf_function("05_UMAP_stromal_MERGE_orig.ident.pdf")
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
    temp_pdf_function("06_UMAP_stromal_MERGE_legend_orig.ident.pdf")
    print(leg)
    dev.off()
  }
  
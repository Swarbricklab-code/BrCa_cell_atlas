library(Seurat)
library(tidyverse)
library(viridis)
library(gridExtra)
library(grid)
library(ggrepel)
library(pheatmap)

options(scipen=10000)

#Function to add flag for pheatmap (genes of interest show only)
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
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
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



#Set wd
gen_dir <- "/Single_cell/breast_atlas/"
setwd(gen_dir)

#Read files
seurat_obj <- readRDS(file = "BrCa_seurat_obj.Rdata")

#Subset out Myeloid clusters and re-ident on cluster level
M_obj <- subset(seurat_obj, ident = c("Myeloid"))
DefaultAssay(M_obj) <- "RNA"
Idents(object = M_obj) <- "celltype_subset"

#Duetere gene list
dutertre_gene_list <- (read.table(paste0("references/dutertre_etal_2019_study_cluster_DGE.csv"), header = T, sep = ","))

#Take top 20 genes for each cluster
dutertre_gene_top20 <- dutertre_gene_list %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_logFC)


#take info of interest
y_newdata <- data.frame(dutertre_gene_top20[c(1,7)])

#Other studies gene list take top 20

#Other studies already taken top 20 gene lists
myeloid_multiple_gene_lists <- (read.table(paste0("references/Myeloid_multiple_study_gene_signature.csv"), header = T, sep = ","))
x_newdata <- data.frame(myeloid_multiple_gene_lists)

#merge two gene of interests and prepare annotations
anno_data <- rbind(x_newdata,y_newdata)
anno_data_unique <- anno_data[!duplicated(anno_data$gene), ]
anno_data_unique <- data.frame(anno_data_unique)
pheatmap_names <- as.data.frame(anno_data_unique$cluster)
rownames(pheatmap_names) <- anno_data_unique$gene
z <- rownames(pheatmap_names)

#pull only genes present within seurat obj
avail_genes <- intersect(x = rownames(M_obj), y = z)


#USe seurat dotplot to generate scaled averaged expression data
dplot <- DotPlot(object = M_obj, features = unique(avail_genes))

#extract matrix and spread, and prepare 
ddata <- as.data.frame(dplot[["data"]])
df <- tidyr::spread(ddata[, c(3,4,5)], id, avg.exp.scaled)
rownames(df) <- df[, 1]
df <- df[, -1]
mat <- as.matrix(df)
mat[is.na(mat)] <- 0

#Change name 
pheatmap_names <- rename(pheatmap_names, c("anno_data_unique$cluster"="Gene_List"))

#order of heatmap
order <- c("LAMP3+ DC","cDC1","cDC2","Cycling Myeloid","LAM 2","LAM 1","Macrophage 2","pDC","Macrophage 3","Macrophage 1","Monocyte 2","Monocyte 3","Monocyte 1")
mat.ordered <- mat[,order]

#colors for pheatmap
cols <- colorRampPalette(brewer.pal(12, "Paired"))
pheatmap_color <- pheatmap_names
mycolors <- cols(length(unique(pheatmap_color$Gene_List)))
names(mycolors) <- unique(pheatmap_color$Gene_List)
mycolors <- list(Gene_List = mycolors)


#generate heatmap
phet<- pheatmap(mat.ordered,
                color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                          "RdYlBu")))(100),
                border_color = NA,
                fontsize = 8,
                cellheight = 2,
                drop_levels = TRUE,
                cluster_cols = F,
                fontsize_row = 7,
                annotation_row = pheatmap_color,
                cellwidth = 7,
                annotation_colors = mycolors)

#Add genes to flag
sigGenes_v <- unique(c("TOP2A","CCR7","APOE","CD1C","LAMP3","SPP1","IL7R",
                       "APOE","MARCO","CXCL11","CXCL10","IL6","FTL","SERPINA1",
                       "CD209","IRF7","LYVE1","SIGLEC1","BIRC3","IRF8","IRF4",
                       "MKI67","TREM2","CLEC9A","EGR1","IL1B","DUSP6","FTL","TNFSF13",
                       "S100A9","FCGR3A","CXCL10","SIGLEC1","APOE","TNFSF10","C1QB",
                       "FABP5","IL1B","XCR1","CADM1","FCER1A","FAM198B","GZMB","IL2RA","IL1B","IDO1"))

#create heatmap repeling genes except ones of interest
flag_phet <- add.flag(phet,
                      kept.labels = sigGenes_v,
                      repel.degree = 0.5)


#Save file
ggsave(flag_phet, dpi = 300, units = "in", width = 20, height = 15,  filename = "/Single_cell/breast_atlas/Myeloid_signature_heatmap.png")
ggsave(flag_phet, dpi = 300, units = "in", width = 20, height = 15,  filename = "/Single_cell/breast_atlas/Myeloid_signature_heatmap.eps")

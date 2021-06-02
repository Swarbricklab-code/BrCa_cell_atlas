library(Seurat)
library(tidyverse)
library(viridis)
library(gridExtra)
library(grid)
library(ggrepel)
library(pheatmap)

options(scipen=10000)

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

#Duetere gene list
myeloid_multiple_gene_lists <- (read.table(paste0("references/myeloid_multiple_study_gene_signature.csv"), header = T, sep = ","))
multiple_myeloid_top20 <- myeloid_multiple_gene_lists %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = gene)

#take info of interest
x_newdata <- data.frame(multiple_myeloid_top20)

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

library(Seurat)
library(tidyverse)
library(viridis)
library(pheatmap)

genes_OI <- read.table("/Single_cell/gene_list/immune_modulator_gene_list.csv", header = T, sep = ",")
seurat_obj <- readRDS(file = paste0(d_in_obj, "BrCa_seurat_obj.Rdata"))
Idents(object = seurat_obj) <- "cell_subset"
DefaultAssay(seurat_obj) <- "RNA"

dplot <- DotPlot(object = seurat_obj, features = c(as.character(unique(genes_OI$gene_name))), cols = world_col, split.by = "subtype")

ddata <- as.data.frame(dplot[["data"]])
df_spread <- tidyr::spread(ddata[, c(3,4,5)], id, avg.exp.scaled)
rownames(df_spread) <- df_spread[, 1]
df_spread <- df_spread[, -1]
mat <- as.matrix(df_spread)

phet<- pheatmap(mat,
                color = inferno(10),
                fontsize = 6,
                cellheight = 6,
                border_color = T,
                fontsize_row = 6,cutree_cols = 10, cutree_rows = 5,
                cellwidth = 6)


ggsave(phet, filename = "pheatmap_all_BrCa.eps", width = 17, height = 16, dpi = 300, device = "eps")

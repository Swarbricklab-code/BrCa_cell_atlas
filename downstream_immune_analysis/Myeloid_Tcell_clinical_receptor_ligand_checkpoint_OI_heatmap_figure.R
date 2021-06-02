library(RColorBrewer)
library(Seurat)
library(pheatmap)

#Set wd
gen_dir <- "/Single_cell/breast_atlas/"
setwd(gen_dir)

###Colours###
world_col <- c("#b2182b", "#1F618D","#f4a582", "#9B59B6","#85929E",
    "#1c9099","#74add1","#053061","#1b7837","#b8e186","#bebada",
    "#fed976","#e7298a","#47E3FF","#FFBF47","#A93226","#270e26",
    "orange","#b8bc53","#5628ce","#fa909c","#8ff331",
    "#FF6347","#6347FF","#556270","#4ECDC4","#C7F464",
    "#FF6B6B","#C44D58","#E3FF47","#FF4787","#771155","#AA4488",
    "#CC99BB","#114477","#4477AA","#77AADD","#117777","#44AAAA",
    "#77CCCC","#777711","#AAAA44","#DDDD77","#774411","#AA7744",
    "#DDAA77","#771122","#AA4455","#DD7788")

###Read files###
#Seurat object
seurat_obj <- readRDS(file = "BrCa_seurat_obj.Rdata")

#Gene list of interest
genes_OI <- read.table("gene_list/clinical_immunegene_of_interest.csv"), header = T, sep = ",")

#read cluster annotation across different tiers of breast cancer atlas
BrCa_subtype_annotation <- read.table("breast_cancer_tiered_cluster_annotation.csv"), header = T, sep = ",")


###Subset out Seurat obj Myeloid and T-cells only
Idents(object = seurat_obj) <- "celltype_major"
T_M <- subset(seurat_obj, ident = c("Myeloid","T-cells"))
Idents(object = T_M) <- "celltype_subset"
DefaultAssay(T_M) <- "RNA"


###Use seurat dotplot function to generate averaged expression values and extract into matrix
dplot <- DotPlot(object = T_M, features = genes_OI$gene_name, split.by = "clinical_subtype")
ddata <- as.data.frame(dplot[["data"]])
dfT <- tidyr::spread(ddata[, c(3,4,5)], id, avg.exp.scaled)
rownames(dfT) <- dfT[, 1]
dfT <- dfT[, -1]
mat <- as.matrix(dfT)


#pheatmap row annotation
pheatmap_names_row <- data.frame(genes_OI$Inhibitory_or_Stimulatory, genes_OI$Receptor_Ligand, genes_OI$Likely_cell_type)
rownames(pheatmap_names_row) <- genes_OI$gene_name
colnames(pheatmap_names_row) <- c("Mod","Flag","ExpCellType")

#pheatmap column annotation
pheatmap_names_col <- data.frame(BrCa_subtype_annotation[])
rownames(pheatmap_names_col) <- rownames(BrCa_subtype_annotation)
colnames(pheatmap_names_col) <- c("Major_CellType","Mincor_CellType","Subset_CellType","BrCa_type")

#######Colouring heatmap
pheatmap_color_row <- pheatmap_names_row
pheatmap_color_col <- pheatmap_names_col

#For rows colors chosen
mycolors1 <- c("#d73027", "#4575b4")
mycolors2 <- c("#542788", "#1b7837")
mycolors3 <- c("#7fc97f", "#beaed4", "#c51b7d")

#Colors For column headers colors array
ncolours1 <- colorRampPalette(brewer.pal(8, "Set2"))
ncolours2 <- colorRampPalette(brewer.pal(12, "Paired"))
ncolours3 <- world_col

#chosen colors
mycolors4 <- ncolours1[1:length(unique(BrCa_subtype_annotation$Major_CellType))]
mycolors5 <- ncolours2[1:length(unique(BrCa_subtype_annotation$Minor_CellType))]
mycolors6 <- ncolours3[1:length(unique(BrCa_subtype_annotation$Subset_CellType))]
mycolors7 <- c("blue","pink","red")


#Row Assignment
names(mycolors1) <- unique(pheatmap_color_row$Mod)
names(mycolors2) <- unique(pheatmap_color_row$Flag)
names(mycolors3) <- unique(pheatmap_color_row$ExpCellType)

#Column assignment
names(mycolors4) <- unique(pheatmap_color_col$Major_CellType)
names(mycolors5) <- unique(pheatmap_color_col$Minor_CellType)
names(mycolors6) <- unique(pheatmap_color_col$Subset_CellType)
names(mycolors7) <- unique(pheatmap_color_col$BrCa_type)




anno <- list(Mod = mycolors,
             Flag = mycolors2,
             ExpCellType = mycolors3,
             Major_CellType = mycolors4,
             Minor_CellType = mycolors5,
             Subset_CellType = mycolors6,
             BrCa_type = mycolors7)

#Generate heatmap
phet <- pheatmap(mat,
                color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                          "RdBu")))(100),
                border_color = F,
                fontsize = 5,
                cellheight = 5,
                drop_levels = TRUE,
                cluster_cols = T,
                fontsize_row = 5,
                annotation_row = pheatmap_names_row,
                annotation_colors = anno,
                annotation_col = pheatmap_names_col,
                cellwidth = 5)
#save heatmap
ggsave(phet, filename = "BrCa_TcellsNmcells_heatmap.pdf", width = 10, height = 10)

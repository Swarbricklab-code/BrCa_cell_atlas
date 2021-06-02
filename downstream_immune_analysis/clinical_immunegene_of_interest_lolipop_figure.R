library(Seurat)
library(MAST)
library(insight)
library(ggpubr)
library(ggplot2)
library(future)
library(future.apply)

options(future.globals.maxSize = 25000 * 1048^2)
plan("multisession", workers = 8)
options(scipen=10000)

#Set wd
gen_dir <- "/Single_cell/breast_atlas/"
setwd(gen_dir)

#Read files
seurat_obj <- readRDS(file = "BrCa_seurat_obj.Rdata")
genes_list_file <- read.table("gene_list/Clinical_immunegene_of_interest.csv"), header = T, sep = ",")


#Subset out Myeloid and T-cells, split by celltype then run DGE for each BrCa subtype comparison
T_M <- subset(seurat_obj, ident = c("Myeloid","T-cells"))
DefaultAssay(T_M) <- "RNA"
Idents(object = T_M) <- "celltype_subset"
sample_list <- SplitObject(T_M, split.by = "celltype_subset")

for(i in 1:length(sample_list)){
  tryCatch({
    name <- names(sample_list[i])
    aaa <- FindMarkers(sample_list[[i]], ident.1 = "HER2+",ident.2 = "ER+", test.use = "MAST",assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    bbb <- FindMarkers(sample_list[[i]], ident.1 = "TNBC",ident.2 = "HER2+", test.use = "MAST",assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    ccc <- FindMarkers(sample_list[[i]], ident.1 = "TNBC",ident.2 = "ER+", test.use = "MAST", assay = "RNA", logfc.threshold = 0.1, min.cells.feature = 1)
    write.csv(aaa, file = paste0("out/DGE/two_way/ER_HER2/",name,"_MAST_across_HER2_ER.csv"), row.names = TRUE)
    write.csv(bbb, file = paste0("out/DGE/two_way/TNBC_HER2/",name,"_MAST_across_TNBC_HER2.csv"), row.names = TRUE)
    write.csv(ccc, file = paste0("out/DGE/two_way/TNBC_ER/",name,"_MAST_across_TNBC_ER.csv"), row.names = TRUE)
  }, error=function(e){print("error")})
}


#For TNBC vs HER2 comparison, generate dataframe of the genes of interest for each celltype comparison
path_file <- "/Single_cell/breast_atlas/out/DGE/two_way/TNBC_HER2/"
dge_files <- list.files(path = path_file)
gene_OI <- genes_list_file$gene_name
file_check <- read.csv(file = paste0(path_file,dge_files[1]),row.names = 1)
checkpoints_significant <- data.frame(matrix(ncol = 5, nrow = 0))


for (i in 1:length(dge_files)){
  tryCatch({
    file_check <- read.csv(file = paste0(path_file,dge_files[i]),row.names = 1)
    name_to_use <- gsub("_MAST_across_TNBC_HER2.csv", "", dge_files[i])
    inter_df <- file_check[which(rownames(file_check) %in% gene_OI),]
    inter_df$celltype <- name_to_use
    inter_df$comparison <- "TNBC_HER2"
    inter_df$gene <- paste0("TNBC vs HER2 ",rownames(inter_df))
    checkpoints_significant <- rbind(checkpoints_significant,inter_df)
  }, error=function(e){print("error")})
}

TNBC_HER2 <- checkpoints_significant


#For TNBC vs ER comparison, generate dataframe of the genes of interest for each celltype comparison
path_file <- "/Single_cell/breast_atlas/out/DGE/two_way/TNBC_ER/"
dge_files <- list.files(path = path_file)
gene_OI <- genes_list_file$gene_name
file_check <- read.csv(file = paste0(path_file,dge_files[1]),row.names = 1)
checkpoints_significant <- data.frame(matrix(ncol = 5, nrow = 0))


for (i in 1:length(dge_files)){
  tryCatch({
    file_check <- read.csv(file = paste0(path_file,dge_files[i]),row.names = 1)
    name_to_use <- gsub("_MAST_across_TNBC_ER.csv", "", dge_files[i])
    inter_df <- file_check[which(rownames(file_check) %in% gene_OI),]
    inter_df$celltype <- name_to_use
    inter_df$comparison <- "TNBC_ER"
    inter_df$gene <- paste0("TNBC vs ER ",rownames(inter_df))
    checkpoints_significant <- rbind(checkpoints_significant,inter_df)
  }, error=function(e){print("error")})
}

TNBC_ER <- checkpoints_significant


#For HER vs ER comparison, generate dataframe of the genes of interest for each celltype comparison
path_file <- "/Single_cell/breast_atlas/out/DGE/two_way/ER_HER2/"
dge_files <- list.files(path = path_file)
gene_OI <- genes_list_file$gene_name
file_check <- read.csv(file = paste0(path_file,dge_files[1]),row.names = 1)
checkpoints_significant <- data.frame(matrix(ncol = 5, nrow = 0))


for (i in 1:length(dge_files)){
  tryCatch({
    file_check <- read.csv(file = paste0(path_file,dge_files[i]),row.names = 1)
    name_to_use <- gsub("_MAST_across_HER2_ER.csv", "", dge_files[i])
    inter_df <- file_check[which(rownames(file_check) %in% gene_OI),]
    inter_df$celltype <- name_to_use
    inter_df$comparison <- "HER2_ER"
    inter_df$gene <- paste0("HER2 vs ER ",rownames(inter_df))
    checkpoints_significant <- rbind(checkpoints_significant,inter_df)
  }, error=function(e){print("error")})
}

ER_HER2 <- checkpoints_significant

#Aggregate all data.frames into one and add NS for any none signficant values
all_data <- rbind(TNBC_HER2,TNBC_ER,ER_HER2)
all_data$sign <- insight::format_p(all_data$p_val_adj,stars_only = TRUE)
all_data <- all_data %>% mutate(sign = ifelse(all_data$p_val_adj>0.05,"ns",sign))
write.csv(all_data, file = "/Single_cell/breast_atlas/out/DGE/BrCa_paired_comparison_DGE_MAST.csv", row.names = T)

#Arrange by high to low Log threshold
plot_data <- all_data
plot_data$gene <- as.factor(plot_data$gene)
plot_data$comparison <- as.factor(plot_data$comparison)
#plot_data <- plot_data %>% group_by(comparison) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
plot_data <- plot_data %>% arrange(desc(avg_log2FC))

#Remove all not signicant
plot_data2 <- plot_data %>% filter(plot_data$sign != "ns")


#best figure, but other options are out there
a <- ggdotchart(plot_data2, x = "gene", y = "avg_log2FC",
                color = "celltype",                          # Color by groups
                palette = map_colours, # Custom color palette
                sorting = "descending",
                add = "segments",                # Sort value in descending order                    # Add segments from y = 0 to dots
                add.params = list(color = "lightgray", size = 1), # Change segment color and size
                group = "comparison",                                # Order by groups
                dot.size = 6,                                 # Large dot size
                # label = round(All_signf$p_val_adj,1),                        # Add mpg values as dot labels
                label = "sign",             # Add mpg values as dot labels
                font.label = list(color = "white", size = 7,vjust = 0.5),               # Adjust label parameters
                ggtheme = theme_pubr()) +
  theme_cleveland() + geom_hline(yintercept = 0, linetype = 1, color = "black") +
  coord_cartesian(ylim = c(-1.4, 1.4)) + facet_wrap(~comparison, nrow=1, scales="free_x")

a <- a+font("xy.text", size = 9, color = "black", face = "bold.italic")

ggsave(a, dpi = 300, units = "in", width = 20, height = 10,  filename = "/Single_cell/breast_atlas/DGE_paired_comparison_clinical_immunegene_of_interest_genes.png")
ggsave(a, dpi = 300, units = "in", width = 20, height = 10,  filename = "/Single_cell/breast_atlas/DGE_paired_comparison_clinical_immunegene_of_interest_genes.eps")



#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_131019"
numcores <- 7
sample_name <- args[1]
subset_data <- as.logical(args[2])
include_t_cells <- as.logical(args[3])
#sample_name <- "CID4461"
#subset_data <- FALSE
#include_t_cells <- TRUE


print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Number cores = ", numcores))
print(paste0("Include T cells? ", as.character(include_t_cells)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(Seurat)
library(infercnv, lib.loc=lib_loc)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(dplyr)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
setwd(in_dir)

if (include_t_cells) {
  out_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/")
} else {
  out_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/")
}

input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))
plot_dir <- paste0(out_dir, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "/tables/")
system(paste0("mkdir -p ", table_dir))
integrated_dir <- paste0("/share/ScratchGeneral/sunwu/projects/", 
  "MINI_ATLAS_PROJECT/Jun2019/04_reclustering_analysis/",
  "run06_v1.2.1/output/Epithelial/02_Rdata/")

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name))


################################################################################
### 0. Define functions ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))


################################################################################
### 1. Generate input matrix and metadata files ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1

# create raw matrix input file and subset if necessary:
count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
if (subset_data) {
  count_df <- count_df[1:500, 1:500]
}
print(
  paste0(
    "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
  )
)

# create metadata df:
print("Creating inferCNV metadata file...")
infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, subset_data=subset_data, 
  count_df, for_infercnv=T)
seurat_10X <- infercnv_metadata$seurat
print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))

# only keep cells in metadata df:
print(paste0("No cells in count df before filtering for those in metadata df = ", 
    ncol(count_df)))
count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
print(paste0("No cells in count df after filtering for those in metadata df = ", 
    ncol(count_df)))

# generate cluster metric plots for epithelial cluster:
epithelial_clusters <- grep("pithelial", unique(infercnv_metadata$metadata$cell_type), value=T)
print(paste0("Epithelial cluster = ", epithelial_clusters))
if (!file.exists(paste0(plot_dir, "metrics_by_epithelial_cluster.png"))) {
  png(paste0(plot_dir, "metrics_by_epithelial_cluster.png"),
    width=14, height=8, res=300, units='in')
    temp_violinplot <- VlnPlot(
      object = seurat_10X,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
      pt.size = 1.5,
      idents = epithelial_clusters
    )
    print(temp_violinplot)
  dev.off()
}

# remove cluster information for epithelial cells:
infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)] <- 
gsub("_[0-9].*$", "", 
  infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)])
# remove CAFs from analysis:
infercnv_metadata$metadata <- infercnv_metadata$metadata[
  grep("CAF", infercnv_metadata$metadata$cell_type, invert=T),
]
# if necessary, remove T cells from analysis:
if (!include_t_cells) {
  infercnv_metadata$metadata <- infercnv_metadata$metadata[
    grep("[t,T][-_][c,C]ell", infercnv_metadata$metadata$cell_type, invert=T),
  ]
}
# collapse all stromal cells into 'stromal' cell type:
infercnv_metadata$metadata$cell_type[
  grep("pithelial", infercnv_metadata$metadata$cell_type, invert=T)
] <- "Stromal"

# if no epithelial clusters present, abort:
if (length(epithelial_clusters) < 1) {
  print(paste0("No epithelial/myoepithelial clusters detected, aborting..."))
} else {
  # write count, metadata files and new seurat object:
  if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
    print("Creating inferCNV raw counts file...")
    write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  }

  if (!file.exists(paste0(input_dir, "metadata.txt"))) {
    write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
      quote=F, sep="\t", col.names=F, row.names=F)
    write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
      "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
    saveRDS(seurat_10X, paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  }
  
  # define normals which will act as InferCNV reference cells:
  normals <- grep(
    "[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned|[u,U]nknown|[t,T]umour|[t,T]umor", 
    unique(infercnv_metadata$metadata$cell_type[
      infercnv_metadata$metadata$cell_ids %in% colnames(count_df)
    ]), value=T, 
    invert=T
  )


  ################################################################################
  ### 2. Run InferCNV ###
  ################################################################################
  
  print(paste0("Normal is: ", normals))

  print("Creating inferCNV object...")
  raw_path <- paste0(input_dir, "input_matrix.txt")
  annotation_path <- paste0(input_dir, "metadata.txt")
  gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
  initial_infercnv_object <- CreateInfercnvObject(
    raw_counts_matrix=raw_path,
    annotations_file=annotation_path,
    delim="\t",
    gene_order_file=gene_path,
    ref_group_names=normals
  )

  print("InferCNV object created, running inferCNV...")
  system.time(
    infercnv_output <- try(
      infercnv::run(
        initial_infercnv_object,
        num_threads=numcores-1,
        out_dir=out_dir,
        cutoff=0.1,
        window_length=101,
        max_centered_threshold=3,
        cluster_by_groups=F,
        plot_steps=F,
        denoise=T,
        sd_amplifier=1.3,
        analysis_mode = "samples"
      )
    )
  )
}

# remove temp files:
system(paste0("rm ", out_dir, "/*infercnv_obj"))
system(paste0("rm ", out_dir, "/*dat"))
system(paste0("rm ", out_dir, "/*preliminary*"))
system(paste0("rm ", out_dir, "/run.final.infercnv_obj"))


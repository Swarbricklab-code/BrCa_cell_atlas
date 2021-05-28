#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

subproject_name <- "brca_mini_atlas_131019"

args = commandArgs(trailingOnly=TRUE)
include_t_cells <- as.logical(args[1])
subset_data <- as.logical(args[2])
subset_samples <- as.logical(args[3])
na_colour <- args[4]
include_normals <- as.logical(args[5])
missing_genes_min_proportion <- as.numeric(args[6])	# only missing genes with presence in at 
													# least x% 
													# missing_genes_inclusion_threshold samples 
													# will be included
cluster_by <- args[7]
type_order <- strsplit(args[8], split = "\\.")[[1]]
cluster_method <- args[9]

include_t_cells <- T
subset_data <- F
subset_samples <- T
na_colour <- "white"
include_normals <- F
missing_genes_min_proportion <- 0.5	# only missing genes with presence in at least x% 
										# missing_genes_inclusion_threshold samples will be 
										#included
cluster_by <- "sc50"
type_order <- strsplit(
  "LumA_SC50.LumB_SC50.Her2_SC50.Basal_SC50.NA", split = "\\."
)[[1]]
cluster_method <- "ward.D2"

heatmap_prefix <- paste0(
  "combined_infercnv_", na_colour, "_missing_values_heatmap"
)

if (include_normals) {
  heatmap_prefix <- paste0(heatmap_prefix, "_normals_included"
}

if (subset_samples) {

  # determine order of samples by subtype:
  ER <- c("CID3941")
  HER2_ER <- c("CID3586")
  HER2 <- c("CID3921")
  TNBC <- c("CID44041")
  sample_names <- c(ER, HER2, HER2_ER, TNBC)

} else {
  # determine order of samples by subtype:
  ER <- c("CID3941", "CID3948", "CID4067", "CID4290A", "CID4461", "CID4463",
    "CID4471", "CID4530N", "CID4535")
  HER2_ER <- c("CID3586", "CID4066")
  HER2 <- c("CID3921", "CID3963", "CID45171")
  TNBC <- c("CID4465", "CID44041", "CID4495", "CID44971", "CID44991", "CID4515",
    "CID4513", "CID4523")
  sample_names <- c(ER, HER2, HER2_ER, TNBC)

}

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample names = ", sample_names))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Subset samples? ", as.character(subset_samples)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(Seurat)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(Polychrome)

project_name <- "identify_epithelial"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

if (include_t_cells) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/")
  out_dir <- paste0(results_dir, 
  	"infercnv/t_cells_included/combined_infercnv/genes_in_at_least_",
  	missing_genes_min_proportion, "_of_samples/clustered_by_", cluster_by, 
    "_using_", cluster_method, "/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/")
  out_dir <- paste0(results_dir, 
  	"infercnv/t_cells_included/combined_infercnv/genes_in_at_least_",
  	missing_genes_min_proportion, "_of_samples/clustered_by_", cluster_by, 
    "_using_", cluster_method,"/")
}

if (subset_data & subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_double_sub/")
  plot_dir <- paste0(out_dir, "plots_double_sub/")
} else if (subset_data) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
} else if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sample_sub/")
  plot_dir <- paste0(out_dir, "plots_sample_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
}

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))

print("Running InferCNV identify normals pipeline on: ")
print(sample_names)


################################################################################
### 0. Define functions and colours ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
create_PAM50_CNV_annotation <- dget(paste0(func_dir, 
  "create_PAM50_CNV_annotation.R"))
get_sample_colours <- function(sample_ids, sample_colours) {
  cols <- all_sample_colours[sample_ids,]$colour
  names(cols) <- sample_ids
  return(cols)
}
sc50_colours <- c("LumA_SC50" = "darkblue", "LumB_SC50" = "cyan", 
  "Her2_SC50" = "pink", "Basal_SC50" = "red", "NA" = "#74add1")
Clinical_IHC_colours <- c("HER2_ER" = "purple", "ER" = "blue", "TNBC" = "red", 
  "HER2" = "pink")

set.seed(9641)
all_sample_colours <- data.frame(
  sample_id = c("CID3586", "CID44041", "CID4461", "CID4463", "CID4465", 
    "CID4471", "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", 
    "CID3921","CID45171", "CID4523", "CID4530N", "CID4535", "CID3941", 
    "CID3948","CID3963", "CID4066", "CID4067", "CID4290A", "CID4398", 
    "CID3838", "CID3946", "CID4040"),
  colour = createPalette(
    26, c("#2A95E8", "#E5629C"), range = c(10, 60), M = 100000
  ),
  stringsAsFactors = F
)
rownames(all_sample_colours) <- all_sample_colours$sample_id
sample_colours <- get_sample_colours(sample_names)


################################################################################
### 1. Load heatmap data and combine ###
################################################################################

print("Loading heatmap and metadata dfs for each sample...")

if (!file.exists(paste0(Robject_dir, "1a.initial_combined_heatmap.Rdata")) | 
  !file.exists(paste0(Robject_dir, "1b.initial_combined_heatmap_metadata.Rdata"))) {

  # create tally table of infercnv genes:
  all_genes <- read.table(
    paste0(ref_dir, "/infercnv_gene_order.txt"),
    header = F,
    sep = "\t",
    as.is = T
   )$V1

  gene_tally <- data.frame(
    row.names = all_genes,
    genes = all_genes,
    tally = rep(0, length(all_genes))
  )
  
  for (s in 1:length(sample_names)) {
    print(sample_names[s])

    # load heatmap and add genes to tally:
    sample_Robject_dir <- paste0(in_path, sample_names[[s]], 
      "/Rdata/")
    heatmap_df <- readRDS(
      paste0(sample_Robject_dir, "heatmap_df.Rdata")
    )
    gene_tally[colnames(heatmap_df),]$tally <- gene_tally[colnames(heatmap_df),]$tally+1
    
    # load metadata:
    heatmap_metadata <- readRDS(
      paste0(sample_Robject_dir, "heatmap_metadata.Rdata")
    )
    if ("sc50" %in% colnames(heatmap_metadata)) {
      heatmap_metadata <- subset(heatmap_metadata, 
        select = c(cell_ids, cell_type, nUMI, nGene, CNA_value, cor.estimate, 
        cor.p.value, normal_cell_call, sc50))
    } else {
      heatmap_metadata <- subset(heatmap_metadata, 
        select = c(cell_ids, cell_type, nUMI, nGene, CNA_value, cor.estimate, 
        cor.p.value, normal_cell_call))
      heatmap_metadata$sc50 <- NA
    }
    

    print(paste0(
      "Are heatmap_metadata rownames in the same order as heatmap_df?? ",
      identical(rownames(heatmap_df), rownames(heatmap_metadata))
    ))
  
    if (s==1) {
      heatmap_dfs <- list(heatmap_df)
      names(heatmap_dfs) <- sample_names[s]
      heatmap_metadatas <- list(heatmap_metadata)
      names(heatmap_metadata) <- sample_names[s]
    } else {
      heatmap_dfs[[s]] <- heatmap_df
      names(heatmap_dfs)[s] <- sample_names[s]
      heatmap_metadatas[[s]] <- heatmap_metadata
      names(heatmap_metadatas)[s] <- sample_names[s]
    }
  }
  
  # determine which genes are in enough samples to be included:
  gene_threshold <- floor(missing_genes_min_proportion*length(sample_names))
  gene_list <- as.character(gene_tally$gene[gene_tally$tally >= gene_threshold])
  
  # add all missing genes as columns to all heatmap dfs:
  print("Adding missing genes and collating heatmap and metadata dfs...")
  complete_heatmap_dfs <- lapply(heatmap_dfs, function(x) {
    
    missing_genes <- gene_list[!(gene_list %in% colnames(x))]
    missing_genes_df <- data.frame(matrix(NA, nrow = nrow(x), 
      ncol = length(missing_genes)))
    colnames(missing_genes_df) <- missing_genes
  
    complete_df <- cbind(x, missing_genes_df)
    m <- match(gene_list, colnames(complete_df))
    complete_df <- complete_df[,m]
  
    return(complete_df)
  
  })
  
  # collate all dfs and check rows line up:
  group_heatmap_df <- do.call(rbind, complete_heatmap_dfs)
  rownames(group_heatmap_df) <- gsub("^.*\\.C", "C", rownames(group_heatmap_df))
  print("Are all rows present in group heatmap df?")
  identical(nrow(group_heatmap_df), sum(unlist(lapply(complete_heatmap_dfs, nrow))))
  print("Are all genes present in group heatmap df?")
  identical(colnames(group_heatmap_df), gene_list)
  
  # subset group_heatmap_df if needed:
  if (subset_data) {
    group_heatmap_df <- group_heatmap_df[,1:300]
  }
  
  # collate heatmap metadata:
  group_heatmap_metadata <- do.call("rbind", heatmap_metadatas)
  rownames(group_heatmap_metadata) <- gsub(
    "^.*\\.", "", rownames(group_heatmap_metadata)
  )

  print(paste0(
    "Are heatmap_metadata rownames still in the same order as heatmap_df?? ",
    identical(rownames(group_heatmap_df), rownames(group_heatmap_metadata))
  ))


  ################################################################################
  ### 2. Format heatmap data ###
  ################################################################################

  # select only group_heatmap_metadata rows in group_heatmap_df and order:
  group_heatmap_metadata <- group_heatmap_metadata[rownames(group_heatmap_df),]
  m <- match(rownames(group_heatmap_df), rownames(group_heatmap_metadata))
  group_heatmap_metadata <- group_heatmap_metadata[m,]
  
  print(paste0(
    "Are heatmap_metadata rownames still in the same order as heatmap_df?? ",
    identical(rownames(group_heatmap_df), rownames(group_heatmap_metadata))
  ))
  
  # if necessary, remove all but cancer cells from group_heatmap_df and 
  # group_metadata_df:
  if (!include_normals) {
    print("Removing normal epithelial cells from heatmap_df and heatmap_metadata...")
    print(paste0("No. of group_heatmap_df rows before removing normals = ", 
      nrow(group_heatmap_df)))
    cancer_ids <- rownames(group_heatmap_metadata)[
      group_heatmap_metadata$normal_cell_call == "cancer"
    ]
    print(paste0("No. of cancer cells = ", length(cancer_ids)))
    group_heatmap_df <- group_heatmap_df[cancer_ids,]
    group_heatmap_metadata <- group_heatmap_metadata[cancer_ids,]
    print(paste0("No. of group_heatmap_df rows after removing normals = ", 
      nrow(group_heatmap_df)))
    
    print(paste0(
      "Are heatmap_metadata rownames still in the same order as heatmap_df?? ",
      identical(rownames(group_heatmap_df), rownames(group_heatmap_metadata))
    ))
  }

  # split heatmap by sample and remove any samples with <2 cells:
  sample_rownames <- split(
    rownames(group_heatmap_df), gsub("_.*$", "", rownames(group_heatmap_df))
  )
  split_heatmap <- lapply(sample_rownames, function(x) group_heatmap_df[x,])
  for (m in 1:length(split_heatmap)) {
    if (nrow(split_heatmap[[m]]) < 2) {
      if (exists("remove_heatmaps")) {
        remove_heatmaps <- c(remove_heatmaps, names(split_heatmap)[m])
      } else {
        remove_heatmaps <- c(names(split_heatmap)[m])
      }
    }
  }
  if (exists("remove_heatmaps")) {
    split_heatmap <- split_heatmap[!(names(split_heatmap) %in% remove_heatmaps)]
    sample_names <- sample_names[!(sample_names %in% remove_heatmaps)]
    for (r in remove_heatmaps) {
      group_heatmap_metadata <- group_heatmap_metadata[
        grep(r, rownames(group_heatmap_metadata), invert=T),
      ]
      group_heatmap_df <- group_heatmap_df[
        grep(r, rownames(group_heatmap_df), invert=T),
      ]
    }
    print(paste0(
      "Are heatmap_metadata rownames still in the same order as heatmap_df?? ",
      identical(rownames(group_heatmap_df), rownames(group_heatmap_metadata))
    ))
  }

  # add sample and subtype columns to heatmap_metadata:
  group_heatmap_metadata$sample <- gsub(
    "_.*$", "", group_heatmap_metadata$cell_id
  )
  group_heatmap_metadata$subtype <- NA
  group_heatmap_metadata$subtype[
    group_heatmap_metadata$sample %in% ER
  ] <- "ER"
  group_heatmap_metadata$subtype[
    group_heatmap_metadata$sample %in% HER2_ER
  ] <- "HER2_ER"
  group_heatmap_metadata$subtype[
    group_heatmap_metadata$sample %in% HER2
  ] <- "HER2"
  group_heatmap_metadata$subtype[
    group_heatmap_metadata$sample %in% TNBC
  ] <- "TNBC"
  
  print(paste0(
      "Are heatmap_metadata rownames still in the same order as heatmap_df?? ",
      identical(rownames(heatmap_df), rownames(heatmap_metadata))
  ))
  
  saveRDS(group_heatmap_df, paste0(Robject_dir, "1a.initial_group_heatmap_df.Rdata"))
  saveRDS(
    group_heatmap_metadata, 
    paste0(Robject_dir, "1b.initial_group_heatmap_metadata.Rdata")
  )

} else {

  group_heatmap_df <- readRDS(paste0(Robject_dir, "1a.initial_group_heatmap_df.Rdata"))
  group_heatmap_metadata <- readRDS(
    paste0(Robject_dir, "1b.initial_group_heatmap_metadata.Rdata")
  )
  print(paste0(
      "Are heatmap_metadata rownames still in the same order as heatmap_df?? ",
      identical(rownames(heatmap_df), rownames(heatmap_metadata))
  ))

}


################################################################################
### 3. Order and recluster data ###
################################################################################

if (!file.exists(paste0(Robject_dir, "2a.reclustered_combined_heatmap.Rdata")) | 
  !file.exists(paste0(Robject_dir, "2b.reclustered_combined_heatmap_metadata.Rdata"))) {

  # re-add sc50 cells:
  sc50_calls <- read.table(paste0(ref_dir, "/sc50_calls.txt"), sep = "\t", 
    header = T, as.is=T)
  colnames(sc50_calls)[1] <- "cell_ids"
  sc50_calls <- subset(sc50_calls, select = c(cell_ids, Calls))
  sc50_calls$Calls <- paste0(sc50_calls$Calls, "50")
  rownames(sc50_calls) <- sc50_calls$cell_ids

  if (subset_samples) {
  	sc50_calls <- sc50_calls[sc50_calls$cell_ids %in% rownames(group_heatmap_metadata),]
  }

  # determine which cells are not present in sc50_calls:
  missing_cells <- group_heatmap_metadata$cell_ids[
    !(group_heatmap_metadata$cell_ids %in% sc50_calls$cell_ids)
  ]
  missing_cell_df <- group_heatmap_metadata[missing_cells,]

  # merge sc50 calls with group_heatmap_metadata:
  group_heatmap_metadata <- subset(group_heatmap_metadata, select = -sc50)
  group_heatmap_metadata <- merge(group_heatmap_metadata, sc50_calls, 
  	by = "cell_ids")
  rownames(group_heatmap_metadata) <- group_heatmap_metadata$cell_ids
  colnames(group_heatmap_metadata) <- gsub(
    "Calls",
    "sc50",
    colnames(group_heatmap_metadata),
  )

  # add missing cells back:
  missing_cell_df$sc50 <- NA
  group_heatmap_metadata <- rbind(
    group_heatmap_metadata,
    missing_cell_df
  )
  print(paste0(
      "Are all cells in group_heatmap_df still in group_heatmap_metadata??? ",
      all(rownames(group_heatmap_df) %in% rownames(group_heatmap_metadata))
  ))

  # make NA values as.character:
  temp_cluster_by <- as.character(
    group_heatmap_metadata[
      , match(cluster_by, colnames(group_heatmap_metadata))
    ]
  )
  temp_cluster_by[is.na(temp_cluster_by)] <- "NA"
  group_heatmap_metadata[
    , match(cluster_by, colnames(group_heatmap_metadata))
  ] <- factor(temp_cluster_by)

  # split metadata and heatmap by cluster_by:
  split_heatmap_metadata <- split(
    group_heatmap_metadata, 
    eval(parse(text = paste0("group_heatmap_metadata$", cluster_by)))
  )
  # reorder split_heatmap_metadata:
  split_heatmap_metadata <- split_heatmap_metadata[
    match(type_order, names(split_heatmap_metadata))
  ]
  split_heatmap_metadata <- split_heatmap_metadata[
    !is.na(names(split_heatmap_metadata))
  ]

  split_heatmap_df <- lapply(split_heatmap_metadata, function(x) {
    return(group_heatmap_df[x$cell_ids,])
  })

  # recluster by cluster_by:
  recluster_split_df <- function(df, method = cluster_method) {
    # remove NAs:
    na_index <- apply(df, 2, function(x) {
      return(all(!(is.na(x))))
    })
    no_na_df <- df[,na_index]
    # scale by rows to mean=0, sd=1:
    scaled_df <- as.matrix(scale(no_na_df))
    # cluster by rows:
    cluster_object <- hclust(dist(scaled_df), method = method)
    cluster_order <- cluster_object$order
    return(df[cluster_object$order,])
  }

  if (!file.exists(paste0(Robject_dir, cluster_method, 
  	"_reclustered_dfs.Rdata"))) {
  	reclustered_dfs <- lapply(
      split_heatmap_df, recluster_split_df, cluster_method
    )
    saveRDS(reclustered_dfs, paste0(Robject_dir, cluster_method, 
  	  "_reclustered_dfs.Rdata"))
  } else {
  	reclustered_dfs <- readRDS(paste0(Robject_dir, cluster_method, 
  	"_reclustered_dfs.Rdata"))
  }
  
  # bind heatmap df back together and order metadata accordingly:
  group_heatmap_df <- do.call("rbind", reclustered_dfs)
  rownames(group_heatmap_df) <- gsub("^.*\\.", "", rownames(group_heatmap_df))

  group_heatmap_metadata <- group_heatmap_metadata[
    rownames(group_heatmap_df),
  ]
  # reorder metadata factor level:
  group_heatmap_metadata[
    , match(cluster_by, colnames(group_heatmap_metadata))
  ] <- factor(
    group_heatmap_metadata[
      , match(cluster_by, colnames(group_heatmap_metadata))
    ], levels = type_order
  )

  print(paste0("Are group_heatmap_df and group_heatmap_metadata rownames identical? ",
    identical(rownames(group_heatmap_df), rownames(group_heatmap_metadata))))

  # create lists of split and combined heatmap dfs/metadatas:
  all_heatmap_dfs <- c(reclustered_dfs, list(group_heatmap_df))
  names(all_heatmap_dfs)[length(all_heatmap_dfs)] <- "all"
  all_heatmap_metadatas <- c(split_heatmap_metadata, list(group_heatmap_metadata))
  names(all_heatmap_metadatas)[length(all_heatmap_metadatas)] <- "all"

  saveRDS(group_heatmap_df, paste0(Robject_dir, "2a.reclustered_combined_heatmap.Rdata"))
  saveRDS(
    group_heatmap_metadata, 
    paste0(Robject_dir, "2b.reclustered_combined_heatmap_metadata.Rdata")
  )
  saveRDS(all_heatmap_dfs, paste0(Robject_dir, 
  	"2c.reclustered_combined_heatmaps_full_and_split.Rdata"))
  saveRDS(
    all_heatmap_metadatas, 
    paste0(Robject_dir, "2d.reclustered_combined_heatmap_metadata_full_and_split.Rdata")
  )

} else {

  group_heatmap_df <- readRDS(paste0(Robject_dir, "2a.reclustered_combined_heatmap.Rdata"))
  group_heatmap_metadata <- readRDS(
    paste0(Robject_dir, "2b.reclustered_combined_heatmap_metadata.Rdata")
  )
  all_heatmap_dfs <- readRDS(paste0(Robject_dir, 
  	"2c.reclustered_combined_heatmaps_full_and_split.Rdata"))
  all_heatmap_metadatas <- readRDS(paste0(Robject_dir, 
  	"2d.reclustered_combined_heatmap_metadata_full_and_split.Rdata")
  )

}


################################################################################
### 3a. Create heatmap row annotations ###
################################################################################

for (j in 1:length(all_heatmap_dfs)) {

  # create subtype annotation:
  subtype_annotation_df <- subset(all_heatmap_metadatas[[j]], select = subtype)
  subtype_annotation_df$subtype <- factor(subtype_annotation_df$subtype, 
    levels = c("ER", "HER2_ER", "HER2", "TNBC"))
  print(paste0("Are subtype_annotation_df and all_heatmap_dfs[[j]] rownames identical? ",
    identical(rownames(all_heatmap_dfs[[j]]), rownames(subtype_annotation_df))))
  subtype_annotation <- Heatmap(
    subtype_annotation_df, 
    col = Clinical_IHC_colours, 
    name = "subtype_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(
        title = "Clinical \nsubtype", title_gp = gpar(fontsize = 20, fontface = "bold"), 
        labels_gp = gpar(fontsize = 20), 
        at = as.character(levels(subtype_annotation_df$subtype))
    )
  )
  
  # create sample annotation:
  sample_annotation_df <- subset(all_heatmap_metadatas[[j]], select = sample)
  sample_annotation_df$sample <- factor(sample_annotation_df$sample,
    levels=names(sample_colours))
  print(paste0("Are sample_annotation_df and all_heatmap_dfs[[j]] rownames identical? ",
    identical(rownames(all_heatmap_dfs[[j]]), rownames(sample_annotation_df))))
  # create annotation:
  sample_annotation <- Heatmap(
    sample_annotation_df, 
    col = sample_colours, 
    name = "sample_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(
        title = "Sample", title_gp = gpar(fontsize = 20, fontface = "bold"), 
        labels_gp = gpar(fontsize = 20), 
        at = as.character(levels(sample_annotation_df$sample))
      )
  )
  

  ################################################################################
  ### 3b. Create heatmap row annotations ###
  ################################################################################
  
  # create sc50 annotation:
  sc50_annotation_df <- subset(all_heatmap_metadatas[[j]], select = sc50)
  #temp_sc50 <- as.character(sc50_annotation_df$sc50)
  #temp_sc50[is.na(temp_sc50)] <- "NA"
  if (names(all_heatmap_metadatas[j]) == "all") {
  	sc50_annotation_df$sc50<- factor(temp_sc50, 
      levels = c("LumA_SC50", "LumB_SC50", "Her2_SC50", "Basal_SC50", "NA"))
    print(paste0("Are sc50_annotation_df and all_heatmap_dfs[[j]] rownames identical? ",
      identical(rownames(all_heatmap_dfs[[j]]), rownames(sc50_annotation_df))))
  }
  
  # create annotation:
  sc50_annotation <- Heatmap(
    sc50_annotation_df, 
    col = sc50_colours, 
    name = "sc50_annotation",
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(
        title = "SC50 call", title_gp = gpar(fontsize = 20, fontface = "bold"), 
        labels_gp = gpar(fontsize = 20), 
        at = as.character(levels(sc50_annotation_df$sc50))
    )
  )
  
  # create normal call annotation if needed:
  if (include_normals) {
    normal_call_annot_df <- subset(all_heatmap_metadatas[[j]], select = normal_cell_call)
    normal_call_annot_df$normal_cell_call <- gsub(
      "normal", "Normal", normal_call_annot_df$normal_cell_call
    )
    normal_call_annot_df$normal_cell_call <- gsub(
      "unassigned", "Unassigned", normal_call_annot_df$normal_cell_call
    )
    normal_call_annot_df$normal_cell_call <- gsub(
      "cancer", "Cancer", normal_call_annot_df$normal_cell_call
    )
    ######
    normal_call_annot_df$normal_cell_call[
      grep("CID3586", rownames(normal_call_annot_df))
    ] <- gsub(
    	"Normal", "Unassigned",
    	normal_call_annot_df$normal_cell_call[
        grep("CID3586", rownames(normal_call_annot_df))
      ]
    )
    ######
    normal_call_annot_df$normal_cell_call <- factor(normal_call_annot_df$normal_cell_call, 
    levels = c("Normal", "Unassigned", "Cancer"))
    normal_call_annotation <- Heatmap(normal_call_annot_df, 
      col = c("Unassigned" = "#E7E4D3", "Normal" = "#1B7837", "Cancer" = "#E7298A"), 
      name = "normal_call_annotation", width = unit(6, "mm"), 
      show_row_names = F, show_column_names = F, 
      show_heatmap_legend = F,
      heatmap_legend_param = list(
        at = as.character(levels(normal_call_annot_df$normal_cell_call))
      )
    )
  }
  
  # create QC annotations:
  nUMI_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      all_heatmap_metadatas[[j]]$nUMI,
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
      all_heatmap_metadatas[[j]]$nGene, name = "nGene",
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
  
  
  ################################################################################
  ### 4. Create heatmap column annotations ###
  ################################################################################
  
  # create CNA annotation:
  CNA_value_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      all_heatmap_metadatas[[j]]$CNA_value,
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
  
  # determine co-ordinates of vertical lines at chromosome borders:
  chr_data <- fetch_chromosome_boundaries(all_heatmap_dfs[[j]], ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
  
  # create PAM50 annotations:
  PAM50_subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")
  METABRIC_CNV_frequencies <- read.table(paste0(ref_dir, 
    "infercnv_metabric_cnv.txt"), header=T, as.is=T, fill=T)
  for ( i in 1:length(PAM50_subtypes) ) {
      print(paste0("Generating ", PAM50_subtypes[i], " CNV plot..."))
      if (i==1) {
        metabric_plots <- 
          list(create_PAM50_CNV_annotation(all_heatmap_dfs[[j]], METABRIC_CNV_frequencies, 
          PAM50_subtypes[i], chr_data$ends, chr_data$lengths))
      } else {
        metabric_plots[[i]] <- 
          create_PAM50_CNV_annotation(all_heatmap_dfs[[j]], METABRIC_CNV_frequencies, 
          PAM50_subtypes[i], chr_data$ends, chr_data$lengths)
      }
  }
  names(metabric_plots) <- PAM50_subtypes
  
  ## create heatmap annotation for gain and loss-associated genes, collated by Niantao
  ## read in CNV_genes
  #CNV_genes <- read.table(paste0(ref_dir, 
  #  "./infercnv_brca_genes_associated_with_CNVs.txt"), header = T, as.is = T)
  ## create CNV_genes annotation:
  #print("Annotating CNV-associated genes...")
  #CNV_genes_annotation <- create_CNV_genes_annotation(all_heatmap_dfs[[j]], CNV_genes)
  
  
  #############################################################################
  ### 5. Prepare heatmap and annotations ###
  ################################################################################
  
  if (!file.exists(paste0(Robject_dir, "final_group_heatmap.Rdata"))) {
    if (!exists("split_heatmap")) {
      sample_rownames <- split(
        rownames(all_heatmap_dfs[[j]]), gsub("_.*$", "", rownames(all_heatmap_dfs[[j]]))
      )
      split_heatmap <- lapply(sample_rownames, function(x) all_heatmap_dfs[[j]][x,])
    }
    
    # scale all samples the same:
    print(paste0("Summaries of split heatmap df before rescaling to min/max values: "))
    print(lapply(split_heatmap, function(x) summary(unlist(x))))
    
    rescaled_split_heatmap <- lapply(split_heatmap, function(x) {
      return(as.data.frame(rescale(as.matrix(x), c(-1, 1))))
    })
    
    print(paste0("Summaries of split heatmap df after rescaling: "))
    print(lapply(rescaled_split_heatmap, function(x) summary(unlist(x))))
    
    
    rescaled_heatmap_df <- do.call("rbind", rescaled_split_heatmap)
    rownames(rescaled_heatmap_df) <- gsub("^.*\\.", "", rownames(rescaled_heatmap_df))
    rescaled_heatmap_df <- rescaled_heatmap_df[rownames(all_heatmap_metadatas[[j]]),]
    print(paste0("Are rescaled_heatmap_df and all_heatmap_metadatas[[j]] rownames identical? ",
      identical(rownames(rescaled_heatmap_df), rownames(all_heatmap_metadatas[[j]]))))
    
    saveRDS(rescaled_heatmap_df, paste0(Robject_dir, "final_heatmap.Rdata"))
  
  } else {
    rescaled_heatmap_df <- readRDS(paste0(Robject_dir, "final_heatmap.Rdata"))
  }
  
  # prepare df for plotting:
  plot_object <- rescaled_heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  # define heatmap colours:
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 0, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")
  
  print("Generating final heatmap...")
  
  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"),
    na_col = na_colour,
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = FALSE,
  #  bottom_annotation = CNV_genes_annotation, bottom_annotation_height = unit(2, "cm"),
  #    gap = unit(1, "cm"),
    heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
    grid_height = unit(2, "cm"), grid_width = unit(1, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 20, fontface = "bold"), labels_gp = gpar(fontsize = 18)),
    use_raster = T, raster_device = c("png")
  )
  
  # determine co-ordinates of horizontal lines at group borders:
  if (names(all_heatmap_metadatas)[j] == "all") {
  	spl_groups <- split(
      eval(parse(text=paste0("all_heatmap_metadatas[[j]]$", cluster_by))),
      eval(parse(text=paste0("all_heatmap_metadatas[[j]]$", cluster_by)))
    )
    spl_groups <- spl_groups[unique(
      eval(parse(text=paste0("all_heatmap_metadatas[[j]]$", cluster_by)))
    )]
    if (length(spl_groups) > 1) {
      for ( n in 1:(length(spl_groups)-1) ) {
        if (n==1) {
          hlines <- c(length(spl_groups[[n]])/length(
          	eval(parse(text=paste0("all_heatmap_metadatas[[j]]$", cluster_by)))
          ))
        } else {
          hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(
          	eval(parse(text=paste0("all_heatmap_metadatas[[j]]$", cluster_by)))
          )
        }
      }
      hlines <- 1-hlines
    } else {
      hlines <- c(length(spl_groups[[1]])/length(
      	eval(parse(text=paste0("all_heatmap_metadatas[[j]]$", cluster_by)))
      ))
    }
  } else {
  	# determine co-ordinates of horizontal lines at group borders:
    spl_groups <- split(group_heatmap_metadata$sample, 
      group_heatmap_metadata$sample)
    spl_groups <- spl_groups[unique(group_heatmap_metadata$sample)]
    if (length(spl_groups) > 1) {
      for ( n in 1:(length(spl_groups)-1) ) {
        if (n==1) {
          hlines <- c(length(spl_groups[[n]])/length(group_heatmap_metadata$sample))
        } else {
          hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(group_heatmap_metadata$sample)
        }
      }
      hlines <- 1-hlines
    } else {
      hlines <- c(length(spl_groups[[1]])/length(group_heatmap_metadata$sample))
    }
  }

  
  
  if (include_normals) {
    ht_list <- subtype_annotation + sample_annotation + sc50_annotation + 
      final_heatmap + normal_call_annotation + CNA_value_annotation + 
      nUMI_annotation + nGene_annotation
  } else {
    ht_list <- subtype_annotation + sample_annotation + sc50_annotation + 
      final_heatmap + CNA_value_annotation + nUMI_annotation + nGene_annotation
  }
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
  )
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster 
  # name (0.00604 units per character):
  longest_cluster_name <- max(
  	nchar(unique(as.character(all_heatmap_metadatas[[j]]$subtype)))
  )
  x_coord <- longest_cluster_name*0.0037

  
  ################################################################################
  ### 6. Plot heatmap and annotations ###
  ################################################################################
  
  if (include_normals) {
    pdf(paste0(plot_dir, heatmap_prefix, "_rescaled_with_normals_",
      names(all_heatmap_metadatas[j]), ".pdf"), height = 21, 
      width = 20)
  } else {
    pdf(paste0(plot_dir, heatmap_prefix, "_rescaled_",
      names(all_heatmap_metadatas[j]), ".pdf"), height = 21, width = 20)
  }
  grid.newpage()
  
    pushViewport(viewport(x = 0, y = 0.3, 
                        width = 0.99, height = 0.68, just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 1, 
          col = "#383838"))
        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-2.1, "mm"), gp=gpar(fontsize=20))
      }
      for ( m in 1:length(hlines) ) {
        grid.lines(c(0, 1), c(hlines[m], hlines[m]), gp = gpar(lwd = 1, col = "#383838"))
      }
    })
    popViewport()
    pushViewport(viewport(x=x_coord + 0.92, y=0.24, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("CNA", rot=65, gp=gpar(fontsize=20))
    popViewport()
    pushViewport(viewport(x=x_coord + 0.94, y=0.24, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("nUMI", rot=65, gp=gpar(fontsize=20))
    popViewport()
    pushViewport(viewport(x=x_coord + 0.957, y=0.24, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("nGene", rot=65, gp=gpar(fontsize=20))
    popViewport()
  
    # plot Normal subtype CNV frequencies:
    pushViewport(viewport(x = x_coord+0.051 + 0.015, y = 0.058,
                          width = 0.82 + 0.015, height = 0.05, just = c("left", "top")))
    grid.draw(metabric_plots[[5]])
    popViewport()
    
    # plot Basal subtype CNV frequencies:
    pushViewport(viewport(x = x_coord+0.056 + 0.015, y = 0.107,
                          width = 0.82+0.024, height = 0.05, just = c("left", "top")))
    grid.draw(metabric_plots[[4]])
    popViewport()
    
    # plot Her2 subtype CNV frequencies:
    pushViewport(viewport(x = x_coord+0.06 + 0.015, y = 0.156, 
                          width = 0.82+0.023, height = 0.05, just = c("left", "top")))
    grid.draw(metabric_plots[[3]])
    popViewport()
    
    # plot LumB subtype CNV frequencies:
    pushViewport(viewport(x = x_coord+0.054 + 0.015, y = 0.205, 
                          width = 0.82+0.026, height = 0.05, just = c("left", "top")))
    grid.draw(metabric_plots[[2]])
    popViewport()
    
    # plot LumA subtype CNV frequencies:
    pushViewport(viewport(x = x_coord+0.054 + 0.015, y = 0.254, 
                          width = 0.82+0.026, height = 0.05, just = c("left", "top")))
    grid.draw(metabric_plots[[1]])
    popViewport()
      
  dev.off()
  
  print(paste0("Heatmap created, output in ", plot_dir))
}


#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))



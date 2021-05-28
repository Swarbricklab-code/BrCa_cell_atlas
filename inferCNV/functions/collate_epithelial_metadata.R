#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_131019"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/",
  subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")

include_t_cells <- TRUE

if (include_t_cells) {
  in_path <- paste0(results_dir, "sup_figure/t_cells_included/")
  adj_path <- paste0(
    results_dir, "sup_figure/t_cells_included/normal_thresholds_adjusted/"
  )
} else {
  in_path <- paste0(results_dir, "sup_figure/t_cells_excluded/")
  adj_path <- paste0(
    results_dir, "sup_figure/t_cells_excluded/normal_thresholds_adjusted/"
  )
}

sample_names <- list.files(in_path, pattern = "CID")
sample_names <- sample_names[!(sample_names %in% "CID4530")]

adj_sample_names <- list.files(adj_path, pattern = "CID")

for (s in 1:length(sample_names)) {

  in_dir <- paste0(in_path, sample_names[s], "/tables/")

  print(paste0("Loading ", sample_names[s], " metadata from ", in_dir, "..."))

  metadata <- read.table(paste0(in_dir, "epithelial_metadata.txt"), header = T,
    sep = "\t", as.is = T)

  metadata <- subset(
    metadata, select = c("cell_ids", "cell_type", "nUMI", 
    "nGene", "CNA_value", "cor.estimate", "cor.p.value", 
    "normal_cell_call", "sc50")
  )

  if (sample_names[s] %in% adj_sample_names) {
    adj_dir <- paste0(adj_path, sample_names[s], "/tables/")
    adj_metadata <- read.table(paste0(adj_dir, "epithelial_metadata.txt"), header = T,
      sep = "\t", as.is = T)
    rownames(adj_metadata) <- adj_metadata$cell_ids
    adj_metadata <- adj_metadata[metadata$cell_ids,]
    colnames(adj_metadata) <- gsub(
      "normal_cell_call", "adj_normal_cell_call",
      colnames(adj_metadata))
    metadata <- cbind(metadata, adj_metadata$adj_normal_cell_call)
    colnames(metadata)[ncol(metadata)] <- "adj_normal_cell_call"
  } else {
    metadata <- cbind(metadata, metadata$normal_cell_call)
    colnames(metadata)[ncol(metadata)] <- "adj_normal_cell_call"
  }

  print(paste0("No. columns = ", ncol(metadata)))

  writeLines("\n")

  if (s==1) {
    all_metadata <- metadata
  } else {
    all_metadata <- rbind(all_metadata, metadata)
  }

}

sc50_calls <- read.table(paste0(ref_dir, "sc50_calls.txt"),
  header = T, sep = "\t", as.is = T)

rownames(all_metadata) <- all_metadata$cell_ids
metadata_sub1 <- all_metadata[sc50_calls$SampleID,]
metadata_sub1$sc50_nov <- sc50_calls$Calls

metadata_sub2 <- all_metadata[
  !(all_metadata$cell_ids %in% sc50_calls$SampleID),
]
metadata_sub2$sc50_nov <- NA

all_metadata_sc50 <- rbind(metadata_sub1, metadata_sub2)

colnames(all_metadata_sc50) <- gsub(
  "sc50$", "sc50_oct", colnames(all_metadata_sc50)
)

write.table(
  all_metadata_sc50, 
  paste0(in_path, "all_infercnv_metadata.txt"),
  col.names = T,
  row.names = F,
  quote = F
)














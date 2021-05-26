#### SEURAT CCA HPC R SCRIPT
#     written by; Sunny Z Wu
#     last modified; 20190628
#
#     This script aligns multiple seurat V3 objects using SEURAT V3
#
### INPUTS (IN ORDER);
        # 01 PROJECT NAME
        # 02 SPECIES (human/mouse)
        # 03 WORKING DIRECTORY (if submitting job from directory use $(pwd))
        # 04 INPUT CSV FILE
            # - input_file.csv
              # sample_id,subtype,molecular_subtype,tissue_site,Ki67,previous_Tx,collection,chemistry,path_to_individual_output_directory_containing_all_samples
              # CIDx,ER+,LuminalA,Primary,50%,NA,Surgery,5_prime,/share/ScratchGeneral/sunwu/MINI_ATLAS_PROJECT/Feb2019_final_primary_set/individual_seurat_objects/output_old/
        # 05 GENES TO PLOT FILE
            # - genes.csv
              # panel,CAFs,immune
              # ACTB,PDGFRB,ITGAM
              # EPCAM,PDGFRA,LY6G6D
              # KRT18,COL1A1,LY6C1
              # KRT5,COL1A2,CD33
              # ESR1,ACTA2,FUT4
              # PGR,MCAM,CD15
              # ERBB2,PDPN,CD68
              # COL1A1,FAP,ADGRE1
              # PDGFRB,CD34,CD274
              # PDGFRA,CXCL12,PDCD1LG2
              # ACTA2,THY1,PDCD1
              # CD34,S100A4,LGALS9
              # PECAM1,,CD80
              # PTPRC,,CD86
              # CD3D,,CD276
              # CD8A,,TNFRSF18
              # FOXP3,,TNFSF18
              # CD19,,
              # JCHAIN,,
              # CD68,,
        # 06 CCA PARAMS FILE
            # seurat_CCA_params_file.csv  
              # SEURAT_PROCESSING_PARAMS_FILE,OPTION
              # #CCA_PARAMS
              # Number_of_CCs_to_compute,50
              # Number_genes_for_integration_anchors_CCA,4000
              # CCs_clustering_1,20
              # CCs_clustering_2,30
              # CCs_clustering_3,50
              # #CLUSTERING_RESOLUTIONS,--
              # RES_1,0.4
              # RES_2,0.8
              # RES_3,1.2
              # regress_cell_cycle?,T
              # #DIFFERENTIAL_GENE_EXPRESSION,--
              # Perform_differential_gene_expression,F
              # minimum_fraction_detected_min.pct,0.5
              # minimum_fraction_diff_min.diff.pct,0.05
              # minimum_threshold_diff_min.thresh.use,0.75
              # #EXPORT_MATRICES,--
              # export_raw_UMI_matrix?,F
              # export_raw_normalised_matrix?,F
# 
### REQUIRES R v.3.5.0
### QSUB ARGUMENTS
#     qsub 
#     -cwd 
#     -pe smp 32 
#     -l h_vmem=200G 
#     -P TumourProgression 
#     -b y 
#     -j y 
#     -V 
#     -N sCCA_sampleID
#     "R CMD BATCH 
#     --no-save 
#     '--args 
#     projectID
#     human
#     /working/directory/path/
#     CCA_input_file.csv
#     gene_set_file.csv
#     seurat_CCA_params_file.csv'
#     /path/to/this/script.R" 
#
### QSUB ONE LINER:
# qsub -cwd -pe smp 32 -l h_vmem=200G -b y -j y -V -P TumourProgression "R CMD BATCH --no-save ‘—-args multiCCA_project human /working/directory/ ./input_file.csv ./genes_file.csv ./CCA_params.csv’ ~/scripts/seurat.R" 

# PARSE ARGUMENTS  ------------------------------------------------------------

# arguments from command line
temp_start_time <- Sys.time()
temp_args <-
  commandArgs(trailingOnly = T)

# 01 PROJECT/SAMPLE NAME
temp_project_name <- temp_args[1]
# 02 SPECIES (human/mouse)
temp_species_type <- temp_args[2]
# 03 WORKING DIRECTORY (if submitting job from directory use $(pwd))
temp_wd <- temp_args[3]
# 04 CCA SAMPLE INPUT FILE
temp_input_file_path <- temp_args[4]
# 05 GENES TO PLOT FILE
temp_gene_plot_file <- temp_args[5]
# 06 SEURAT CCA PARAMS FILE
temp_params_file <- temp_args[6]
# 07 PATH TO OBJECTS OUTPUT FILES
temp_path_to_individiual_objects <- temp_args[7]

  # temp_species_type <- "human"
  # temp_input_file_path <- "./input.file.csv"
  # temp_gene_plot_file <- "seurat_gene_input_file.csv"
  # temp_params_file <- "seurat_CCA_params.csv"
    
# LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(plyr) # always load before plyr to avoid problems in depreciated functions
library(dplyr)
library(Matrix)
library(cowplot)
library(tidyr)
library(eply)
library(rowr)

# SET UP AND FUNCTIONS ------------------------------------------------------------------

# set WD
setwd(temp_wd)

# error processing
options(error = expression(NULL))

# sub-directory outputs
dir.create("Output")
dir.create("Output/Gene_expression_and_stats")
dir.create("Output/Figures")
dir.create("Output/Rdata")

# PNG function
temp_png_function <- 
  function(x) {
  png(
    file = (x), 
    width = 14, 
    height = 8, 
    res = 300, 
    units = 'in'
  )
}

temp_png_function_small <-
  function(x) {
    png(
      file = (x), 
      width = 5, 
      height = 5, 
      res = 300, 
      units = 'in'
    )
  }

temp_png_function_big <-
  function(x) {
    png(
      file = (x), 
      width = 22, 
      height = 12, 
      res = 300, 
      units = 'in'
    )
  }

temp_png_function_longer <-
  function(x) {
    png(
      file = (x), 
      width = 14, 
      height = 12, 
      res = 300, 
      units = 'in'
    )
  }

# PDF function 
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 14,
      height = 8.5, 
      onefile = T 
    )
  }

# FUNCTION TO ADD METADATA
temp_function_add_metadata <-
  function(x, y) {
    temp_seurat_name <-
      as.character(temp_input_file[(x), (y)])
  }

# FUNCTION TO ADD PATH
temp_function_add_path <-
  function(x) {
    
    temp_name <- temp_object_list$SAMPLENAME[(x)]
    
    paste0(temp_path_to_individiual_objects, 
           "seurat_",
           temp_name,
           "/Output/Rdata/03_seurat_object_processed.RData"
           )
  }

# first letter up for mouse gene conversion
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# ANALYSIS PARAMETERS -----------------------------------------------------

temp_params <- 
  read.csv(temp_params_file, 
           row.names = "SEURAT_PROCESSING_PARAMS_FILE")

# run on cluster parallel 
temp_run_parallel <-
  as.logical(temp_params["Run_parallel",]) 

temp_run_parallel_cores <- 
  as.numeric(as.character((temp_params["Number_of_cores_requested",]))) 

if(temp_run_parallel){
  library(future)
}

# number of CCs to compute
temp_num_cc_compute <- 
  as.numeric(as.character((temp_params["Number_of_CCs_to_compute",]))) 

# number of variable genes for CCA
temp_variable_genes_for_CCA <- 
  as.character((temp_params["Number_genes_for_integration_anchors_CCA",]))

if(!temp_variable_genes_for_CCA == "ALL") {
  temp_variable_genes_for_CCA <- as.numeric(temp_variable_genes_for_CCA)
}

# integrate all genes by intersect?
temp_integrate_all_genes_by_intersect <- 
  as.logical(temp_params["Integrate_all_genes_by_intersect",]) 

# PCA analysis and t-SNE resolution (reduce for less Jackstraw processing time)
temp_jackstraw_cutoff <- 
  as.logical(temp_params["Jackstraw_PC_cutoff",]) 

temp_PCs_to_compute <-
  as.numeric(as.character((temp_params["Number_of_PCs_to_compute",]))) 
temp_PC_A <- 
  as.numeric(as.character((temp_params["jackstraw_pval_cutoff_1",])))  
temp_PC_B <- 
  as.numeric(as.character((temp_params["jackstraw_pval_cutoff_2",])))  
temp_PC_C <- 
  as.numeric(as.character((temp_params["jackstraw_pval_cutoff_3",])))  

# MANUAL PCA CUTOFFS
if(temp_jackstraw_cutoff == F) {
  temp_PC_A <- 
    (1:as.numeric(as.character((temp_params["PC_CUTOFF_1",])))) 
  temp_PC_B <-
    (1:as.numeric(as.character((temp_params["PC_CUTOFF_2",])))) 
  temp_PC_C <- 
    (1:as.numeric(as.character((temp_params["PC_CUTOFF_3",])))) 
}

# CLUSTERING RESOLUTIONS
temp_res_1 <- 
  as.numeric(as.character((temp_params["RES_1",])))  
temp_res_2 <- 
  as.numeric(as.character((temp_params["RES_2",])))  
temp_res_3 <- 
  as.numeric(as.character((temp_params["RES_3",])))  

#Differential gene expression params
# Do differential expression?
temp_do_differential_expression <- 
  as.logical(temp_params["Perform_differential_gene_expression",]) 
temp_min.pct <- 
  as.numeric(as.character((temp_params["minimum_fraction_detected_min.pct",]))) 
temp_min.diff.pct <- 
  as.numeric(as.character((temp_params["minimum_fraction_diff_min.diff.pct",])))  
temp_thresh.use <- 
  as.numeric(as.character((temp_params["minimum_threshold_diff_min.thresh.use",])))  

#exporting matrix and associated metadata (T or F)
temp_export_normalised_matrix <- 
  as.logical(temp_params["Export_raw_UMI_matrix",]) 
temp_export_raw_UMI_matrix <- 
  as.logical(temp_params["Export_raw_normalised_matrix",]) 

# GARNETT
  # plot garnett calls from individual cell type annotations
temp_plot_garnett_calls <- 
  as.logical(temp_params["plot_garnett",]) 
  # Performed integrated garnett cluster calls
temp_integrated_garnett_calls <- 
  as.logical(temp_params["integrated_cluster_calls_from_raw_garnett",]) 

# INTEGRATE WITHOUT SEURAT CCA INTEGRATION (SEURAT MERGE PCA)
temp_run_PCA_without_integration <- 
  as.logical(temp_params["additional_seurat_merge_analysis",]) 

# MANUAL PCA CUTOFFS
if(temp_run_PCA_without_integration) {
  temp_PCMERGE_A <- 
    (1:as.numeric(as.character((temp_params["PC_MERGE_CUTOFF_1",])))) 
  temp_PCMERGE_B <-
    (1:as.numeric(as.character((temp_params["PC_MERGE_CUTOFF_2",])))) 
  temp_PCMERGE_C <- 
    (1:as.numeric(as.character((temp_params["PC_MERGE_CUTOFF_3",])))) 
}

# LOAD DATA ------------------------------------------------------------

# input file
temp_input_file_raw <- 
  read.csv(temp_input_file_path)

# temp_path_to_individiual_objects <- as.character(temp_input_file_raw$path_to_individual_output_directory[1])

temp_input_file <- temp_input_file_raw[ , !names(temp_input_file_raw) %in% ("path_to_individual_output_directory")]
temp_colnames <- names(temp_input_file)

# number of seurat objects to align
temp_number_of_objects <-
  length(temp_input_file$sample_id)

# read sample IDs and path to seurat objects
temp_object_list <- list(NULL)
temp_object_list$SAMPLENAME <-
  lapply((1:temp_number_of_objects),
         temp_function_add_metadata,
         y = "sample_id")

for (i in c(1:temp_number_of_objects)) {
  n <- paste0("temp_seurat_name_",
              i)
  assign(n,
         as.character(temp_object_list$SAMPLENAME[i]))
}

# add path to seurat object
temp_object_list$PATH <-
  lapply((1:temp_number_of_objects),
         temp_function_add_path)


for(i in temp_colnames[!temp_colnames %in% "sample_id"]) {
  
  temp_object_list[[paste0(i)]] <-
    lapply((1:temp_number_of_objects),
           temp_function_add_metadata,
           y = i)
}


# read Rdata for each input seurat object
for (i in c(1:temp_number_of_objects)) {
  n <- paste0("temp_seurat_object_", i)
  assign(n,
         readRDS(as.character(temp_object_list$PATH[i])))
}

# ADD CLINICAL METADATA  --------------------------------------------------

for (i in c(1:temp_number_of_objects)) {
  
  temp_seurat_metadata <- get(paste0("temp_seurat_object_", 
                                            i))
  
  for(type in temp_colnames[!temp_colnames %in% "sample_id"]) {
    temp_seurat_metadata@meta.data[[type]] <- as.character(temp_object_list[[type]][i])
  }
  
  n <- paste0("temp_seurat_object_", 
              i)
  assign(n, 
         temp_seurat_metadata)
  
  rm(temp_seurat_metadata)
}

# FIND INTEGRATION ANCHORS -----------------------------------------------------------------

# Note that FindIntegrationAnchors will, by default, use the top 2000 variable genes
# that are identified in common between all individual datasets. If temp_variable_genes_for_CCA is
# set higher than 2000 genes (for a larger batch corrected matrix), then the following
# code will rerun FindVariableFeatures for all individual datasets to identify the set 
# number of genes.
if(temp_variable_genes_for_CCA > 2000){
  for(i in c(1:temp_number_of_objects)){
    
    temp_seurat_object <- get(paste0("temp_seurat_object_", 
                                       i))
    
    temp_seurat_object <- FindVariableFeatures(
      object = temp_seurat_object,
      do.plot = F,
      nfeatures=temp_variable_genes_for_CCA
    )
    
    n <- paste0("temp_seurat_object_", 
                i)
    assign(n, 
           temp_seurat_object)
    
    rm(temp_seurat_object)
    
  }
}

temp_sample_list <-
  mget(ls(pattern = "temp_seurat_object_*"))

# Find integration anchors
if(temp_run_parallel){
  plan("multiprocess",
       workers = temp_run_parallel_cores)
  options(future.globals.maxSize = 10 * 1024^3)
}

temp_anchor_start <- Sys.time()
seurat_10X_anchors <- FindIntegrationAnchors(object.list = temp_sample_list, 
                                       dims = 1:temp_num_cc_compute,
                                       anchor.features= temp_variable_genes_for_CCA)
temp_anchor_finish <- Sys.time()
if(temp_run_parallel){
  plan("sequential")
}

# union of genes
# temp_genes_union <- NULL
# for(i in c(1:length(temp_sample_list))){
#   temp_genes <- rownames(temp_sample_list[[i]])
#   temp_genes_union <- union(temp_genes_union, temp_genes)
# }

# intersect
if(temp_integrate_all_genes_by_intersect) {
for(i in c(1:length(temp_sample_list))){
  temp_genes <- rownames(temp_sample_list[[i]])
  n <-
    paste0("temp_genes_sample_", i)
  assign(n,
         temp_genes)
}
temp_genes_list <-
  mget(ls(pattern = "temp_genes_sample_*"))
temp_genes_interesect <-
  Reduce(intersect, temp_genes_list)
}

# Return batch corrected seurat object
if(temp_integrate_all_genes_by_intersect == F) {
  temp_integration_start <- Sys.time()
  seurat_10X_integrated <- IntegrateData(anchorset = seurat_10X_anchors, 
                                         dims = 1:temp_num_cc_compute)
  temp_integration_finish <- Sys.time()
}

if(temp_integrate_all_genes_by_intersect) {
  temp_integration_start <- Sys.time()
  seurat_10X_integrated <- IntegrateData(anchorset = seurat_10X_anchors, 
                                         dims = 1:temp_num_cc_compute,
                                         features.to.integrate = temp_genes_interesect)
  temp_integration_finish <- Sys.time()
}  

DefaultAssay(object = seurat_10X_integrated) <- "integrated"

# scale data
seurat_10X_integrated <- ScaleData(object = seurat_10X_integrated, 
                                   verbose = FALSE)

# PCA --------------------------------------------------------------

seurat_10X_integrated <- RunPCA(object = seurat_10X_integrated, 
                                npcs = temp_PCs_to_compute, 
                                verbose = FALSE)


# SAVE OBJECT PRE-DIMENSTIONAL REDUCTION ----------------------------------

saveRDS(seurat_10X_integrated,
        "Output/Rdata/01_seurat_CCA_aligned.Rdata")

# JACKSTRAW ANALYSIS ------------------------------------------------------

if(temp_jackstraw_cutoff) {
  
  if(temp_run_parallel){
    plan("multiprocess", 
         workers = temp_run_parallel_cores)
    options(future.globals.maxSize = 10 * 1024^3)
  }
  
  temp_jackstraw_start <- Sys.time()
  print(temp_jackstraw_start)
  
  seurat_10X_integrated <- 
    JackStraw(object = seurat_10X_integrated, 
              dims=temp_PCs_to_compute,
              verbose=T)
  
  seurat_10X_integrated <- 
    ScoreJackStraw(seurat_10X_integrated,
                   dims=1:temp_PCs_to_compute)
  
  temp_jackstraw_finish <- Sys.time()
  print(temp_jackstraw_finish)
  if(temp_run_parallel){
    plan("sequential")
  }
}

if(temp_jackstraw_cutoff) {

  temp_PC <- seurat_10X_integrated@reductions$pca@jackstraw$overall.p.values
  temp_PC <- as.data.frame(temp_PC)
  temp_PC <- arrange(temp_PC,
                     Score)
  
  temp_PCA_df <- NULL
  for(i in c("A", "B", "C")) {
    
    temp_jackstraw_cutoff <- 
      get(paste0("temp_PC_",i))
    
    print(paste0("jackstraw significance cutoff ",i," = ",temp_jackstraw_cutoff))
    
    temp_PC_filtered <-
      temp_PC[temp_PC$Score < temp_jackstraw_cutoff,]
    
    temp_PC_filtered <- 
      temp_PC_filtered$PC
    
    print(paste0("Number of PCs used = ", length(temp_PC_filtered)))
    
    temp_PCA_df_subset <- 
      data.frame(PCA = temp_PC_filtered)
    names(temp_PCA_df_subset) <- paste0("PCA_cutoff_",temp_jackstraw_cutoff)
    
    temp_PCA_df <- 
      cbind.fill(temp_PCA_df,
                 temp_PCA_df_subset,
                 fill=NA)
    
    n <- 
      paste0("temp_PC_",i)
    assign(n, 
           temp_PC_filtered)
  }
  temp_PCA_df <- 
    temp_PCA_df[,!colnames(temp_PCA_df) %in% "init", drop=F]
  write.csv(temp_PCA_df, 
            "Output/Gene_expression_and_stats/02_PCAs_used_clustering.csv", 
            row.names = F)
}

saveRDS(seurat_10X_integrated,
        "Output/Rdata/02_seurat_CCA_aligned_jackstraw.Rdata")

# DIMENSIONAL REDUCTION TSNE & UMAP  --------------------------------------------------------------

# TSNE
dir.create("Output/Figures/TSNE/")
for(i in c("A", "B", "C")) {
  
  temp_PC <- get(paste0("temp_PC_",i))
  seurat_10X_integrated <-
    RunTSNE(object = seurat_10X_integrated,
            dims = temp_PC,
            reduction.key = paste0("TSNE",i,"_"),
            reduction.name = paste0("TSNE",i))
}

# UMAP
dir.create("Output/Figures/UMAP/")
for(i in c("A", "B", "C")) {
  
  temp_PC <- get(paste0("temp_PC_",i))
  seurat_10X_integrated <- 
    RunUMAP(seurat_10X_integrated, 
            dims = temp_PC,
            reduction.key = paste0("UMAP",i,"_"),
            reduction.name = paste0("UMAP",i),
            verbose = F
    )
}

# plot
for(i in c("A", "B", "C")) {
  for(dr in c("TSNE", "UMAP")) {
    
    temp_png_function(paste0("Output/Figures/",dr,"/",dr,"_PC_",
                             i,
                             ".png"))
    temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      label.size = 4,
      group.by = "orig.ident",
      pt.size = 1,
      reduction = paste0(dr,i)
    )
    print(temp_dimplot)
    dev.off()
  }
}

# plot split UMAP by sample ID 
for(i in c("A", "B", "C")) {
  if(temp_number_of_objects > 10){
    temp_pt_size <- 0.5
  } else temp_pt_size <- 1
  for(dr in c("TSNE", "UMAP")) {
    
    temp_png_function_big(paste0("Output/Figures/",dr,"/SPLIT_",dr,"_PC_",
                             i,
                             ".png"))
    temp_dimplot <- DimPlot(
      object = seurat_10X_integrated,
      split.by = "orig.ident",
      pt.size = temp_pt_size,
      reduction = paste0(dr,i)
    )
    print(temp_dimplot)
    dev.off()
  }
}


# plot by metadata
dir.create("Output/Figures/Metadata_plots/")
for(type in temp_colnames[!temp_colnames %in% "sample_id"]) {
  dir.create(paste0("Output/Figures/Metadata_plots/", type, "/"))
  
  for(i in c("A", "B", "C")) {
    if(temp_number_of_objects > 10){
      temp_pt_size <- 1
    } else temp_pt_size <- 2
    for(dr in c("TSNE", "UMAP")) {
      
      temp_png_function_big(paste0("Output/Figures/Metadata_plots/",type, "/",
                                   dr,"_PC_",i,
                                   ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated,
        split.by = type,
        group.by = type,
        pt.size = temp_pt_size,
        reduction = paste0(dr,i)
      )
      print(temp_dimplot)
      dev.off()
    }
    }
  }

# plot Garnett calls + the highest res calculated seurat cluster calls.
if(temp_plot_garnett_calls) {
  dir.create(paste0("Output/Figures/Metadata_plots/", "Garnett", "/"))
  temp_colnames <- 
    tail(grep("garnett_seurat_cluster_call_", colnames(seurat_10X_integrated@meta.data), value=T),n=2)
  for(garnettcalltype in c("garnett_call", "garnett_call_ext", "garnett_call_major", "garnett_call_ext_major",temp_colnames)) {
    dir.create(paste0("Output/Figures/Metadata_plots/", "Garnett", "/", garnettcalltype))
    for(i in c("A", "B", "C")) {
    if(temp_number_of_objects > 10){
      temp_pt_size <- 0.25
    } else temp_pt_size <- 0.5
    
    for(dr in c("TSNE", "UMAP")) {
      temp_png_function(paste0("Output/Figures/Metadata_plots/",
                                   "Garnett", "/", garnettcalltype, "/",
                                   dr, "_PC_", i,
                                   ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated,
        group.by = garnettcalltype,
        pt.size = temp_pt_size,
        reduction = paste0(dr,i),
        label=T
      )
      print(temp_dimplot)
      dev.off()
      
      temp_png_function_big(paste0("Output/Figures/Metadata_plots/",
                                   "Garnett", "/", garnettcalltype, "/SPLIT_",
                                   dr, "_PC_", i,
                                   ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated,
        split.by = garnettcalltype,
        pt.size = temp_pt_size,
        reduction = paste0(dr,i),
        label=T
      )
      print(temp_dimplot)
      dev.off()
    }
    }
  }
}

# READ INPUT GENE SETS -----------------------------------------------------------------

# detected genes in seurat object
temp_genes.detected <- rownames(seurat_10X_integrated)

# read in gene input file
temp_gene_plot <- read.csv(temp_gene_plot_file)

# read marker gene sets
for(i in c(1:ncol(temp_gene_plot))){
  n <- paste0("temp_markers_",
              i)
  assign(n,
         na.omit(temp_gene_plot[i]))
  
  m <- paste0("temp_markers_set_name_",
              i)
  assign(m,
         colnames(temp_gene_plot[i]))
}

# species conversion 
# mouse
if (temp_species_type == "mouse") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_convert <- 
      get(paste0("temp_markers_",
                 i))
    colnames(temp_genes_to_convert) <- "gene"
    
    n <- paste0("temp_mouse_markers_",
                i)
    assign(n,
           firstup(tolower(temp_genes_to_convert$gene))
    )
  }
}
# human
if (temp_species_type == "human") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_convert <- 
      get(paste0("temp_markers_",
                 i))
    colnames(temp_genes_to_convert) <- "gene"
    
    n <- paste0("temp_human_markers_",
                i)
    assign(n,
           temp_genes_to_convert$gene)
  }
}

# filter gene sets for genes that are detected
# mouse
if (temp_species_type == "mouse") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_filter <- 
      get(paste0("temp_mouse_markers_",
                 i))
    temp_genes_to_filter <- 
      data.frame(genes = temp_genes_to_filter)
    temp_genes_to_filter <- 
      subset(temp_genes_to_filter,
             genes %in% temp_genes.detected)
    n <- paste0("temp_mouse_markers_",
                i)
    assign(n,
           temp_genes_to_filter)
  }
}

# human
if (temp_species_type == "human") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_filter <- 
      get(paste0("temp_human_markers_",
                 i))
    temp_genes_to_filter <- 
      data.frame(genes = temp_genes_to_filter)
    temp_genes_to_filter <- 
      subset(temp_genes_to_filter,
             genes %in% temp_genes.detected)
    n <- paste0("temp_human_markers_",
                i)
    assign(n,
           temp_genes_to_filter)
  }
}





# EXPRESSION OF GENE MARKERS AND METRICS --------------------------------------------------

dir.create("Output/Figures/Gene_markers/")
dir.create("Output/Figures/Gene_markers/TSNE/")
dir.create("Output/Figures/Gene_markers/UMAP/")

for(j in c(1:ncol(temp_gene_plot))) {
  for(i in c("A", "B", "C")) {
    
    temp_markers_set_name <- get(paste0("temp_markers_set_name_",
                                        j))
    
    temp_markers_to_use <- get(paste0("temp_",
                                      temp_species_type,
                                      "_markers_",
                                      j))
    
    temp_markers_to_use <- (temp_markers_to_use$genes)
    temp_markers_to_use <- levels(droplevels(temp_markers_to_use))
    
    for(dr in c("TSNE", "UMAP")) {
      
      temp_png_function_big(paste0("Output/Figures/Gene_markers/",dr,"/FeaturePlot_0",
                                   j, 
                                   "_",
                                   temp_markers_set_name,
                                   "_markers_", 
                                   i,
                                   ".png"))
      
      temp_featureplot <- FeaturePlot(
        object = seurat_10X_integrated,
        features = temp_markers_to_use,
        pt.size = 1,
        reduction = paste0(dr,i),
        min.cutoff = 0)
      print(temp_featureplot)
      dev.off()
    }
  }
}

# CLUSTERING  --------------------------------------------------------------

temp_PC_res_dataframe <- 
  data.frame(row.names = row.names(seurat_10X_integrated@meta.data))

for(i in c("A", "B", "C")) {
  for(j in c(temp_res_1, temp_res_2, temp_res_3)) {
    
    temp_res_use <- as.numeric(gsub("res.",
                                    "",
                                    j))
    temp_dims.use <- get(paste0("temp_PC_",i))
    temp_graph_name <- paste0("int_PCA_", i)
    
    seurat_10X_integrated <- 
      FindNeighbors(object = seurat_10X_integrated,
                    dims = temp_dims.use,
                    graph.name = temp_graph_name
      )
    if(temp_run_parallel){
      plan("multiprocess", 
           workers = temp_run_parallel_cores)
      options(future.globals.maxSize = 10 * 1024^3)
    }
    seurat_10X_integrated <-
      FindClusters(
        object = seurat_10X_integrated,
        resolution = temp_res_use,
        graph.name = temp_graph_name
      )
    if(temp_run_parallel){
      plan("sequential")
    }
    temp_colname <- paste0(temp_graph_name, "_res.", temp_res_use)
    
    temp_PC_res_dataframe <- cbind(temp_PC_res_dataframe,
                                   (seurat_10X_integrated@meta.data[,temp_colname]))
    
    colnames(temp_PC_res_dataframe)[(ncol(temp_PC_res_dataframe))] <- temp_colname
    
    for(dr in c("TSNE", "UMAP")) {
      temp_png_function(paste0("Output/Figures/",dr,"/",dr,"_PC_",
                               i, "_", "res", temp_res_use,
                               ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated,
        label.size = 4,
        pt.size = 1,
        reduction = paste0(dr,i),
        group.by = temp_colname,
        label = T
      )
      print(temp_dimplot)
      dev.off()
    }
    
  }
}

seurat_10X_integrated <- 
  AddMetaData(seurat_10X_integrated,
              metadata = temp_PC_res_dataframe)


# DIFFERENTIAL EXPRESSION & HEATMAPS --------------------------------------

if(temp_do_differential_expression) {
  dir.create("Output/Figures/Heatmaps/")
  
    for(i in c("A", "B", "C")) {
      for(j in c(temp_res_1, temp_res_2, temp_res_3)) {
    
          temp_res_use <- as.numeric(gsub("res.",
                                          "",
                                          j))
          temp_graph_name <- paste0("int_PCA_", i)
          
          temp_colname <- paste0(temp_graph_name, "_res.", temp_res_use)
          
          Idents(object = seurat_10X_integrated) <- temp_colname
          
            if(temp_run_parallel == F){
              
              Idents(object = seurat_10X_integrated) <- temp_colname
              
              options(digits = 4)
              temp_cluster.allmarkers <- FindAllMarkers(
                only.pos = T,
                object = seurat_10X_integrated,
                min.pct = temp_min.pct, 
                logfc.threshold = temp_thresh.use,
                min.diff.pct = temp_min.diff.pct, 
                test.use = 'MAST', 
                print.bar = F
              )
              
              temp_cluster.allmarkers <- arrange(temp_cluster.allmarkers,
                                                 (cluster),
                                                 desc(avg_logFC))
              
              write.csv(temp_cluster.allmarkers,
                        file = paste0("Output/Gene_expression_and_stats/04_FindAllMarkers_PC_",
                                      i,
                                      "_res.",
                                      temp_res_use, 
                                      ".csv"))
            }
            if(temp_run_parallel){
              plan("multiprocess", 
                   workers = temp_run_parallel_cores)
              options(future.globals.maxSize = 10 * 1024^3)
              
              temp_cluster.allmarkers <- NULL
              for(clusterID in unique(seurat_10X_integrated@meta.data[,temp_colname])) {
                temp_cluster.allmarkers_subset <- FindMarkers(
                  ident.1 = clusterID,
                  only.pos = T,
                  object = seurat_10X_integrated,
                  min.pct = temp_min.pct, 
                  logfc.threshold = temp_thresh.use,
                  min.diff.pct = temp_min.diff.pct, 
                  test.use = 'MAST', 
                  print.bar = F,
                  min.cells.group=1
                )  
                temp_cluster.allmarkers_subset$cluster <- 
                  paste0(clusterID)
                temp_cluster.allmarkers_subset$gene <- 
                  row.names(temp_cluster.allmarkers_subset)
                temp_cluster.allmarkers <- 
                  rbind(temp_cluster.allmarkers,
                        temp_cluster.allmarkers_subset)
              }
              temp_cluster.allmarkers <- arrange(temp_cluster.allmarkers,
                                                 (cluster),
                                                 desc(avg_logFC))
              
              write.csv(temp_cluster.allmarkers,
                        file = paste0("Output/Gene_expression_and_stats/04_FindAllMarkers_PC_",
                                      i,
                                      "_res.",
                                      temp_res_use, 
                                      ".csv"))
            }
            if(temp_run_parallel){
              plan("sequential")
            }
            
            temp_genes_for_heatmap <- 
              (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(10,
                                                                       avg_logFC))
            
            seurat_10X_integrated@meta.data[,temp_colname] <- 
              as.numeric(seurat_10X_integrated@meta.data[,temp_colname])
            
            seurat_10X_integrated@meta.data[,temp_colname] <-
              factor(seurat_10X_integrated@meta.data[,temp_colname],
                     levels=sort(unique(seurat_10X_integrated@meta.data[,temp_colname])))
            
            temp_png_function_big(paste0("Output/Figures/Heatmaps/DoHeatMap_PC_",
                                         i,
                                         "_res.",
                                         temp_res_use,
                                         ".png"))
            print(DoHeatmap(
              object = seurat_10X_integrated,
              features = temp_genes_for_heatmap$gene,
              group.by = temp_colname,
              size = 0.3
            ))
            dev.off()
    }
  }
}

# INTEGRATED GARNETT CLUSTER CALLS ----------------------------------------

if(temp_integrated_garnett_calls) {
  temp_new_idents_append_COMBINED <- 
    data.frame(row.names = row.names(seurat_10X_integrated@meta.data))
  temp_table <- seurat_10X_integrated@meta.data
  
  for(i in c("A", "B", "C")) {
      for(j in c(temp_res_1, temp_res_2, temp_res_3)) {
        
        temp_res_use <- as.numeric(gsub("res.",
                                        "",
                                        j))
        
        temp_colname <- paste0("int_PCA_", i, "_res.", temp_res_use)
        
        temp_new_idents <- data.frame(row.names = row.names(temp_table),
                                      original_garnett_call = temp_table[,"garnett_call"],
                                      old_ident = temp_table[,temp_colname],
                                      new_ident_major = NA,
                                      new_ident_subset = NA)
        temp_idents <- 
          as.numeric(unique(seurat_10X_integrated@meta.data[,temp_colname]))

        for(ident in unique(temp_new_idents$old_ident)) {
          # print(ident)
          
          # major lineages
          temp_table_subset_major <- 
            temp_new_idents[temp_new_idents$old_ident == ident,]
          
          # epithelial collapse
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "Epithelial cycling"] <- "Epithelial"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "Epithelial Luminal Mature"] <- "Epithelial"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "Epithelial Luminal Progenitor"] <- "Epithelial"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "Basal_Myoepithelial"] <- "Epithelial"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "Epithelial"] <- "Epithelial"
          
          # t-cell collapse
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "T-Regulatory cells"] <- "T cells"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "CD4 T cells"] <- "T cells"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "T-cells cycling"] <- "T cells"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "CD8 T cells"] <- "T cells"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "T cells"] <- "T cells"
          temp_table_subset_major$original_garnett_call[temp_table_subset_major$original_garnett_call == "NK Cells"] <- "T cells"
          
          # table      
          temp_subset_table_major_lineage <- 
            as.data.frame(table(temp_table_subset_major$original_garnett_call))
          # drop Unknown 
          temp_subset_table_major_lineage <- 
            temp_subset_table_major_lineage[!temp_subset_table_major_lineage$Var1 == "Unknown",]
          
          # ID
          temp_id_major <- 
            temp_subset_table_major_lineage$Var1[which.max(temp_subset_table_major_lineage$Freq)]
          temp_id_major <- 
            levels(droplevels(temp_id_major))
          temp_new_idents$new_ident_major[temp_new_idents$old_ident == ident] <- temp_id_major
          
          # subset lineages
          temp_table_subset_subset <- 
            temp_new_idents[temp_new_idents$old_ident == ident,]
          
          # table
          temp_subset_subset_lineage_table <- 
            as.data.frame(table(temp_table_subset_subset$original_garnett_call))
          
          temp_subset_subset_lineage_table <- 
            temp_subset_subset_lineage_table[!temp_subset_subset_lineage_table$Var1 == "Epithelial",]
          
          temp_subset_subset_lineage_table <- 
            temp_subset_subset_lineage_table[!temp_subset_subset_lineage_table$Var1 == "T cells",]
          
          # drop unknown
          temp_subset_subset_lineage_table <- 
            temp_subset_subset_lineage_table[!temp_subset_subset_lineage_table$Var1 == "Unknown",]
          
          # ID
          temp_id_subset <- 
            temp_subset_subset_lineage_table$Var1[which.max(temp_subset_subset_lineage_table$Freq)]
          temp_id_subset <- 
            levels(droplevels(temp_id_subset))
          temp_new_idents$new_ident_subset[temp_new_idents$old_ident == ident] <- temp_id_subset
          
        }
        
        temp_new_idents_append <- data.frame(row.names = row.names(temp_new_idents),
                                             temp_col_major = temp_new_idents$new_ident_major,
                                             temp_col_subset = temp_new_idents$new_ident_subset)
        
        colnames(temp_new_idents_append) <- c(paste0("int_cluster_major_PC_",i,"_res.",temp_res_use),
                                              paste0("int_cluster_subset_PC_",i,"_res.",temp_res_use))
        
        temp_new_idents_append_COMBINED <- cbind(temp_new_idents_append_COMBINED,
                                                 temp_new_idents_append)
      }
    }
    
  seurat_10X_integrated <- AddMetaData(seurat_10X_integrated, 
                              metadata = temp_new_idents_append_COMBINED)
}

# plot integrated cluster calls
if(temp_integrated_garnett_calls) {
  
    dir.create(paste0("Output/Figures/Integrated_cluster_calls/"))
  
    for(garnettcalltype in c("int_cluster_major_PC_", "int_cluster_subset_PC_")) {
      
      if(garnettcalltype == "int_cluster_major_PC_") {
        temp_folder_name <- "Major_lineage"
      }
      if(garnettcalltype == "int_cluster_subset_PC_") {
        temp_folder_name <- "Subset_lineage"
      }
      dir.create(paste0("Output/Figures/Integrated_cluster_calls/", temp_folder_name))
      
      for(i in c("A", "B", "C")) {
        if(temp_number_of_objects > 10){
          temp_pt_size <- 0.25
        } else temp_pt_size <- 0.5
          for(j in c(temp_res_1, temp_res_2, temp_res_3)) {
            for(dr in c("TSNE", "UMAP")) {
          temp_png_function(paste0("Output/Figures/Integrated_cluster_calls/",temp_folder_name, "/",
                                       dr, "_PC_", i, "_res.", j,
                                       ".png"))
          temp_dimplot <- DimPlot(
            object = seurat_10X_integrated,
            group.by = paste0(garnettcalltype, i, "_res.", j),
            pt.size = temp_pt_size,
            label = T,
            reduction = paste0(dr,i)
          )
          print(temp_dimplot)
          dev.off()
          
          #split UMAP
          temp_png_function(paste0("Output/Figures/Integrated_cluster_calls/",temp_folder_name, "/SPLIT_",
                                   dr, "_PC_", i, "_res.", j,
                                   ".png"))
          temp_dimplot <- DimPlot(
            object = seurat_10X_integrated,
            split.by = paste0(garnettcalltype, i, "_res.", j),
            pt.size = temp_pt_size,
            reduction = paste0(dr,i)
          )
          print(temp_dimplot)
          dev.off()
        }
      }
    }
    }
  }

# SEURAT MERGE INTEGRATION ------------------------------------------------

# reprocess data without integration (RNA Assay)
if(temp_run_PCA_without_integration){
  
  DefaultAssay(seurat_10X_integrated) <- "RNA"
  
  seurat_10X_integrated <- FindVariableFeatures(
    object = seurat_10X_integrated,
    do.plot = F
  )
  
  seurat_10X_integrated <- ScaleData(
    object = seurat_10X_integrated,
    vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
    display.progress = F
  )
  
  seurat_10X_integrated <- RunPCA(
    object = seurat_10X_integrated,
    npcs = temp_PCs_to_compute
  )
  
}

# recluster using different dims and low resolution 0.4
if(temp_run_PCA_without_integration){
  for(i in c("A", "B", "C")) {
    
    temp_dims.use <- get(paste0("temp_PCMERGE_",i))
    
    seurat_10X_integrated <-
      RunTSNE(object = seurat_10X_integrated,
              dims = temp_dims.use,
              reduction.key = paste0("TSNEMERGE",i,"_"),
              reduction.name = paste0("TSNEMERGE",i))
    
    seurat_10X_integrated <- 
      RunUMAP(seurat_10X_integrated, 
              dims = temp_dims.use,
              reduction.key = paste0("UMAPMERGE",i,"_"),
              reduction.name = paste0("UMAPMERGE",i),
              verbose = F
      )
    
    temp_graph_name <- paste0("int_PCAMERGE_", i)
    
    seurat_10X_integrated <- 
      FindNeighbors(object = seurat_10X_integrated,
                    dims = temp_dims.use,
                    graph.name = temp_graph_name
      )
    
    # clustering at 0.4
    seurat_10X_integrated <-
      FindClusters(
        object = seurat_10X_integrated,
        graph.name = temp_graph_name
      )
    
  }
  
}

# plot
if(temp_run_PCA_without_integration){
  dir.create("Output/Figures/Merge_analysis/")
  dir.create(paste0("Output/Figures/Merge_analysis/","TSNE","/"))
  dir.create(paste0("Output/Figures/Merge_analysis/","UMAP","/"))
  
  for(i in c("A", "B", "C")) {
    for(dr in c("TSNE", "UMAP")) {
      
      temp_png_function(paste0("Output/Figures/Merge_analysis/",dr,"/",dr,"_PC_",
                               i,
                               ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated,
        group.by = "orig.ident",
        pt.size = 1,
        reduction = paste0(dr,"MERGE",i)
      )
      print(temp_dimplot)
      dev.off()
      
      temp_png_function(paste0("Output/Figures/Merge_analysis/",dr,"/",dr,"_PC_",
                               i, "_", "res", 0.8,
                               ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X_integrated,
        label.size = 4,
        pt.size = 1,
        reduction = paste0(dr,"MERGE",i),
        group.by = paste0(temp_graph_name, "_res.", 0.8),
        label = T
      )
      print(temp_dimplot)
      dev.off()
      
    }
  }
  
}


# SAVE CCA ALIGNED OBJECT --------------------------------------------------------------------

Idents(object = seurat_10X_integrated) <- paste0("int_PCA_C_res.",temp_res_1)

saveRDS(seurat_10X_integrated,
        "Output/Rdata/03_seurat_CCA_aligned_processed.Rdata")

# FILTERED SUMMARY STATS AND PLOTS ----------------------------------------

dir.create("Output/Figures/Metrics/")

# filtered metrics
temp_filtered_data_summary <-
  as.data.frame(cbind((as.matrix(summary(
    seurat_10X_integrated@meta.data$nFeature_RNA
  ))),
  (as.matrix(summary(
    seurat_10X_integrated@meta.data$nCount_RNA
  ))),
  (as.matrix(
    summary(seurat_10X_integrated@meta.data$percent.mito)
  ))))

colnames(temp_filtered_data_summary) <-
  c("#GENES",
    "#UMIS",
    "#MITO%")

temp_filtered_data_summary$CELLNUMBER <- ncol(seurat_10X_integrated)
temp_filtered_data_summary$CELLNUMBER[2:6] <- "NA" 

write.csv(temp_filtered_data_summary, 
          "Output/Gene_expression_and_stats/01_EMPTYDROPS_FILTERED_METRICS.csv", 
          quote = F, 
          row.names = T)

Idents(object = seurat_10X_integrated) <- "orig.ident"

# Plots
temp_png_function("Output/Figures/Metrics/01_FILTERED_VlnPlot_nGene_nUMI_percent_mito.png")
VlnPlot(
  object = seurat_10X_integrated,
  features = c("nFeature_RNA",
               "nCount_RNA",
               "percent.mito"),
  pt.size = 0.25,
  group.by = "orig.ident"
)
dev.off()

temp_png_function("Output/Figures/Metrics/02_FILTERED_GenePlot_nGene_nUMI.png")
FeatureScatter(object = seurat_10X_integrated, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               pt.size = 1,
               group.by = "orig.ident")
dev.off()
temp_png_function("Output/Figures/Metrics/03_FILTERED_GenePlot_nUMI_percent_mito.png")
FeatureScatter(object = seurat_10X_integrated, 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mito",
               pt.size = 1, 
               group.by = "orig.ident")
dev.off()
temp_png_function("Output/Figures/Metrics/04_FILTERED_GenePlot_nGene_percent_mito.png")
FeatureScatter(object = seurat_10X_integrated, 
               feature1 = "nFeature_RNA", 
               feature2 = "percent.mito",
               pt.size = 1, 
               group.by = "orig.ident")
dev.off()

# EXPORTING EXPRESSION MATRIX ------------------------------------------------------

#Raw UMI matrix
if (temp_export_raw_UMI_matrix) {
  dir.create("Output/Matrices_and_metadata")
  temp_matrix <- GetAssayData(object = seurat_10X_integrated, 
                              assay = "RNA", 
                              slot = "data")
  
  temp_matrix.sparse <- Matrix(temp_matrix , sparse = T )
  
  write.table(
    temp_matrix.sparse,
    file = "Output/Matrices_and_metadata/01_normalized_expression_sparse_matrix.txt",
    sep = "\t",
    col.names = T
  )
}

#Normalised matrix
if (temp_export_normalised_matrix) {
  dir.create("Output/Matrices_and_metadata")
  temp_matrix <- GetAssayData(object = seurat_10X_integrated, 
                              assay = "RNA", 
                              slot = "counts")
  
  temp_matrix.sparse <- Matrix(temp_matrix, 
                               sparse = T )
  write.table(
    temp_matrix.sparse,
    file = "Output/Matrices_and_metadata/02_raw_UMI_count_sparse_matrix.txt",
    sep = "\t",
    col.names = T
  )
}

#Associated Metadata
if (temp_export_raw_UMI_matrix) {
temp_table <- seurat_10X_integrated@meta.data
temp_table_charactered <- apply(temp_table,2,as.character)
temp_table_charactered <- cbind(barcode = row.names(temp_table),
                                temp_table_charactered)


write.table(temp_table_charactered,
            file = "Output/Matrices_and_metadata/03_metadata.txt",
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F
)
}

# FINISH ------------------------------------------------------------------

temp_finish_time <- Sys.time()

temp_timing_df <- data.frame(total_process_time_mins = as.numeric(temp_finish_time - temp_start_time, units = "mins"),
                             find_anchors_time_mins = as.numeric(temp_anchor_finish - temp_anchor_start, units = "mins"),
                             integrate_data_time_mins = as.numeric(temp_integration_finish - temp_integration_start, units = "mins"))

write.csv(temp_timing_df,
          "Output/Gene_expression_and_stats/03_run_times.csv")

print(temp_timing_df)

#clean environment
rm(list = ls(pattern = "temp"))
rm(seurat_10X_integrated)




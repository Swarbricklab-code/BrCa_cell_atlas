# BUILDING GARNETT CELL TYPE CALLER
# Sunny Wu
# 20190609
# 
# USING NEW MARKER LIST SIGNATURES FROM:
# XCELL PAPER 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1
#
# brca signatures 
# https://pubmed.ncbi.nlm.nih.gov/19648928/
#
#
# SETUP -------------------------------------------------------------------

library(garnett)
library(Seurat)
library(org.Hs.eg.db)

# DIRECTORIES -------------------------------------------------------------

dir.create("Output")
setwd("Output")

# LOAD DATA ---------------------------------------------------------------

seurat_10X <- 
  readRDS("01_seurat_CCA_aligned_processed.Rdata")

DefaultAssay(object = seurat_10X) <- "RNA"

# CREATE MONOCLE OBJECT ---------------------------------------------------

temp_raw_data <- GetAssayData(object = seurat_10X, 
                              slot = "counts")

# temp_raw_data <- as.data.frame(temp_raw_data)
# temp_raw_data <- as(as.matrix(temp_raw_data), "sparseMatrix")

pd <- new("AnnotatedDataFrame", 
          data = seurat_10X@meta.data)

fData <- data.frame(gene_short_name = row.names(temp_raw_data), row.names = row.names(temp_raw_data))
fd <- new("AnnotatedDataFrame", data = fData)

lowerDetectionLimit <- 0

if(all(temp_raw_data == floor(temp_raw_data))) {
  expressionFamily <- negbinomial.size()
} else if(any(data < 0)){
  expressionFamily <- uninormal()
} else {
  expressionFamily <- tobit()
}

monocle_cds <- newCellDataSet(temp_raw_data,
                              phenoData = pd, 
                              featureData = fd,
                              lowerDetectionLimit=lowerDetectionLimit,
                              expressionFamily=expressionFamily)

monocle_cds <- estimateSizeFactors(monocle_cds)


saveRDS(monocle_cds,
        "monocle_cds_mini_atlas23.Rdata")

# CHECK MARKERS --------------------------------------------------------

# monocle_cds <- readRDS("monocle_cds_mini_atlas23.Rdata")

# PDF function 
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 8,
      height = 120, 
      onefile = T 
    )
  }

marker_file_path <- "./breast_cancer_TME_markers_v4_XCell_revised.txt"
marker_check <- check_markers(monocle_cds, 
                              marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL",
                              classifier_gene_id_type= "SYMBOL")
temp_pdf_function("Output/04_marker_check_refined.pdf")
plot_markers(marker_check, label_size=2.5)
dev.off()

# TRAIN CLASSIFIER ON ATLAS DATA --------------------------------------------------------

start_time <- Sys.time()
set.seed(260)
temp_classifier <- train_cell_classifier(cds = monocle_cds,
                                         marker_file = marker_file_path,
                                         db= org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         marker_file_gene_id_type = "SYMBOL",
                                         classifier_gene_id_type = "SYMBOL",
                                         cores = 32,
                                         num_unknown=5000,
                                         max_training_samples=5000)
finish_time <- Sys.time()

# SAVE CLASSIFIER ---------------------------------------------------------

saveRDS(temp_classifier,
        "garnett_classifier.Rdata")

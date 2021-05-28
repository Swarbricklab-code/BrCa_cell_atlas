prepare_infercnv_metadata <- function(seurat_object, subset_data = FALSE, count_df, for_infercnv=TRUE) {
  # annotate cell identities using garnett_call_ext_major metadata column:
  if (length(grep(seurat_object@meta.data$garnett_call_ext_major[1], 
    Idents(seurat_object)[1])) == 0) {

    annotated_idents <- gsub(
      " ", "_", paste0(seurat_object@meta.data$garnett_call_ext_major, " ", 
      Idents(seurat_object))
    )
    # remove sample id if present:
    annotated_idents <- gsub("CID.*_", "", annotated_idents)
    Idents(seurat_object) <- factor(annotated_idents)
  
    # create infercnv metadata:
    temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
      cell_type = Idents(seurat_object), stringsAsFactors = F)
    temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]
  
    # throw unassigned cells and CAFs (as they may be malignant):
    temp_metadata <- temp_metadata[grep("[u,U]nknown|[u,U]nassigned|CAF", 
      temp_metadata$cell_type, invert=T),]
  
    # only include cells present in count_df:
    print(paste0("No cells in metadata df before filtering for those in count df = ", 
      nrow(temp_metadata)))
    temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
    print(paste0("No cells in metadata df after filtering for those in count df = ", 
      nrow(temp_metadata)))

    if (for_infercnv) {
      # label cells in clusters of < 2 cells as 'outliers' as these will break InferCNV:
      temp_metadata$cell_type <- as.character(temp_metadata$cell_type)
      cluster_list <- split(temp_metadata, temp_metadata$cell_type)
      metadata_outliers_labelled <- do.call(
        "rbind", lapply(cluster_list, function(x) {
          if (nrow(x) < 2) {
            x$cell_type <- gsub("_[0-9].*$", "_outlier", x$cell_type)
          }
          return(x)
        })
      )
  
      # remove outlier cell types with <2 cells:
      cluster_list2 <- split(metadata_outliers_labelled, 
        metadata_outliers_labelled$cell_type)
      i=1
      metadata_final <- do.call(
        "rbind", lapply(cluster_list2, function(x) {
          if (nrow(x) < 2) {
            print(paste0(cluster_list[[i]]$cell_type[1], 
              " removed as contained <2 cells"))
            i <<- i+1
            return(NULL)
          } else {
            i <<- i+1
            return(x)
          }
        })
      )
    } else {
      metadata_final <- temp_metadata
    }
    rownames(metadata_final) <- metadata_final$cell_ids 
    
    # record number per cell type:
    number_per_cell_type <- as.data.frame(table(metadata_final$cell_type))
  
    # create and label results list:
    result_list <- list(metadata_final, number_per_cell_type, seurat_object)
    names(result_list) <- c("metadata", "number_per_group", "seurat")
  
    return(result_list)
  } else {
    print("Cell types already annotated...")
  }
}
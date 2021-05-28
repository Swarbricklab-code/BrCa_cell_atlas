create_array_CNV_annotation <- function(df, CNV_df) {
  # define sample name from df rownames:
  sample_name <- gsub("_.*$", "", rownames(df))

  # select CNV information for sample:
  array_CNVs <- CNV_df[,colnames(CNV_df) %in% sample_name]
  names(array_CNVs) <- rownames(CNV_df)

  # convert allele number to total copy number values and cap 
  # at +3 as occurs in InferCNV:
  array_CNVs <- array_CNVs/2
  array_CNVs[array_CNVs > 3] <- 3
  array_CNVs <- array_CNVs[colnames(df)]
  print(paste0(
  	"No. genes without array CNV information = ", length(which(is.na(array_CNVs)))
  ))
  names(array_CNVs) <- colnames(df)
  array_CNVs <- as.data.frame(array_CNVs)
  array_CNVs <- as.data.frame(t(array_CNVs))

  # create array CNV heatmap annotation:
  array_CNV_heatmap <- Heatmap(
  	array_CNVs, name=paste0("array_heatmap"),
  	col = colorRamp2(c(min(array_CNVs[!is.na(array_CNVs)]), 1, 
    max(array_CNVs[!is.na(array_CNVs)])), 
    c("#00106B", "white", "#680700"), space = "sRGB"), na_col="white",
    cluster_columns = F, cluster_rows = F, show_row_dend = FALSE,
    show_row_names = F, show_column_names = F,
    #heatmap_legend_param = list(title = "SNP array\nCNV", legend_direction = "horizontal"),
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png")
  )

  # determine average infercnv value vector:


  res <- list(array_CNVs, array_CNV_heatmap)
  names(res) <- c("array_CNVs", "array_CNV_heatmap")
  return(res)
}
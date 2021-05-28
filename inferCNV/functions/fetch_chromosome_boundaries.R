fetch_chromosome_boundaries <- function(df, ref_dir) {
  # load in gene annotations:
  gene_order <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"), 
    header=F, as.is=T)
  # subset to only include genes in df:
  genes <- gene_order[gene_order$V1 %in% colnames(df),1:2]
  # remove genes on chrM and chrY from df:
  genes_to_keep <- genes$V1[genes$V2 == "chrM|chrY"]
  df <- df[,!(colnames(df) %in% genes_to_keep)]
  # split and determine lengths of each chromosome:
  chr_split <- split(genes$V2, genes$V2)
  chr_lengths <- unlist(lapply(chr_split, length))
  chr_lengths <- chr_lengths[!(names(chr_lengths) %in% c("chrY", "chrM"))]
  # order chromosomes:
  chr_lengths <- c(
    chr_lengths[
      order(
        as.numeric(
          as.character(
            gsub(
              "chr", "", names(chr_lengths[!(names(chr_lengths) %in% c("chrX"))])
            )
          )
        )
      )
      ], chr_lengths[names(chr_lengths) %in% c("chrX")])
  
  for ( l in 1:length(chr_lengths) ) {
    if (l==1) {
      chr_ends <- c(chr_lengths[l])
    } else {
      chr_ends[l] <- chr_ends[l-1] + chr_lengths[l]
    }
  }
  names(chr_ends) <- names(chr_lengths)
  # find total length of plot:
  total_length <- ncol(df)
  # define vertical line co-ordinates for each chromosomal end:
  end_pos <- chr_ends/total_length
  # define chromosome label co-ordinates for midpoint of each chromosome:
  # find centre of each chromosome:
  for ( i in 1:length(chr_lengths) ) {
    if (i==1) {
      lab_pos <- c(chr_lengths[i]/2)
    } else {
      lab_pos[i] <- chr_ends[i-1] + (chr_lengths[i]/2)
    }
  }
  lab_pos <- lab_pos/chr_ends[length(chr_ends)]
  names(lab_pos) <- names(chr_lengths)
  result_list <- list(lengths = chr_lengths, ends = chr_ends, end_pos = end_pos, 
    lab_pos = lab_pos)

  return(result_list)
}

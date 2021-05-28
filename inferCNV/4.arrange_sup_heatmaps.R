#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Arrange CNV heatmaps with samples, arrays, annotations and chromosomes
# labelled, 12 per page ####

Rstudio <- TRUE

project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_131019"
include_t_cells <- TRUE

# state samples with matched arrays and subsets in order:
array_samples <- c("CID4463", "CID3921", "CID4066", "CID3963", "CID44991",
                   "CID4515", "CID4513", "CID4523")
luminal <- c("CID3941", "CID3948", "CID4067", "CID4290A", "CID4461",
             "CID4465", "CID4471", "CID4530N", "CID4535")
HER2 <- c("CID3586", "CID45171")
TNBC <- c("CID44041", "CID4495", "CID44971", "CID4515")
sample_names <- c(array_samples, luminal, HER2, TNBC)
sample_names_1 <- sample_names[1:18]
sample_names_2 <- sample_names[19:22]

# array_samples <- c("CID4066")
# sample_names <- array_samples

print(paste0("Subproject name = ", subproject_name))

library(png)
library(grid)
library(gridExtra)
library(rlist)

if (Rstudio) {
  home_dir <- "/Users/jamestorpy/clusterHome/"
  library(rlist)
} else {
  home_dir <- "/share/ScratchGeneral/jamtor/"
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  library(rlist, lib.loc=lib_loc)
}

project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/",
                      subproject_name, "/")
results_dir <- paste0(project_dir, "results/")

if (include_t_cells) {
  in_path <- paste0(results_dir, "sup_figure/t_cells_included/")
} else {
  in_path <- paste0(results_dir, "sup_figure/t_cells_excluded/")
}

print(paste0("In path = ", in_path))


################################################################################
### 1. Load heatmap images for each sample as list of rasterGrobs ###
################################################################################

for (s in 1:length(sample_names_1)) {
  
  print(paste0("Loading ", sample_names_1[s], " heatmap png..."))
  in_dir <- paste0(in_path, "/", sample_names_1[s], "/plots/")
  
  if (s==1) {
    heatmap_images_1 <- list(
      rasterGrob(
        readPNG(
          paste0(in_dir, "infercnv_plot_no_labels.png")
        ),
        width = 1,
        height = unit(5, "cm")
      )
    )
    
  } else {
    
    heatmap_images_1 <- list.append(
      heatmap_images_1,
      rasterGrob(
        readPNG(
          paste0(in_dir, "infercnv_plot_no_labels.png")
        ),
        width = 1,
        height = unit(5, "cm")
      )
    )
  }
  names(heatmap_images_1)[s] <- sample_names_1[s]
  
}

for (s in 1:length(sample_names_2)) {
  
  print(paste0("Loading ", sample_names_2[s], " heatmap png..."))
  in_dir <- paste0(in_path, "/", sample_names_2[s], "/plots/")
  
  if (s==1) {
    heatmap_images_2 <- list(
      rasterGrob(
        readPNG(
          paste0(in_dir, "infercnv_plot_no_labels.png")
        ),
        width = 1,
        height = unit(5, "cm")
      )
    )
    
  } else {
    
    heatmap_images_2 <- list.append(
      heatmap_images_2,
      rasterGrob(
        readPNG(
          paste0(in_dir, "infercnv_plot_no_labels.png")
        ),
        width = 1,
        height = unit(5, "cm")
      )
    )
  }
  names(heatmap_images_2)[s] <- sample_names_2[s]
  
}


################################################################################
### 2. Arrange and print heatmaps onto 2 pages ###
################################################################################


pdf(paste0(in_path, "combined_figure_1.pdf"), height = 11.69, width = 8.27)
  do.call("grid.arrange", c(heatmap_images_1, ncol=3))
dev.off()

pdf(paste0(in_path, "combined_figure_2.pdf"), height = 11.69, width = 8.27)
  do.call("grid.arrange", c(heatmap_images_2, ncol=3))
dev.off()

#convert pdf to png:
system(paste0("for p in ", in_path, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
              "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
              "convert -density 150 ", in_path, "$f -quality 90 ", in_path, "$new; done"))


# PROCESSING OF VISIUM DATASET SUSING STUTILITY
# Sunny Z. Wu
# 
#
# Run in screen 
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_spatial
# R
# 
#
# 01: SETUP---------------------------------------

# setup
library(zeallot)
library(STutility)
library(ggplot2)
library(magrittr)
library(Seurat)
library(magick)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# DIRECTORY
dir.create("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/data_processing/")
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/data_processing/")

# 02: PREPARE INPUT DATA --------------------------------------------------

temp_csv <- read.csv("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/Meta_Data.csv")
infiles <- list.files("/share/ScratchGeneral/sunwu/projects/spatial/data_Zenodo_miniatlas/", full.names = T, recursive = T)
infiles_h5 <- list.files("/share/ScratchGeneral/sunwu/projects/spatial/data_Zenodo_miniatlas/raw_feature_bc_matrix_h5", full.names = T, recursive = T)

infoTable <- NULL
for(sample in unique(temp_csv$Clinical_Case)){
  print(sample)
  infiles_sample <- grep(sample,infiles,value = T)
  infiles_sample_h5 <- grep(sample,infiles_h5,value = T)
  
  temp_df <- data.frame(samples=as.character(infiles_sample_h5),
                       imgs=as.character(grep("tissue_hires_image",infiles_sample,value = T)),
                       spotfiles=as.character(grep("tissue_positions",infiles_sample,value = T)),
                       json=as.character(grep("json",infiles_sample,value = T)),
                       patientid=sample,
                       subtype=temp_csv[temp_csv$Clinical_Case==sample,"subtype"])
  infoTable <- rbind(infoTable, temp_df)
}

for(col in colnames(infoTable)){
  infoTable[,col] <- as.character(infoTable[,col])
}

write.csv(infoTable,
          "infoTable.csv")



# 03: LOAD AND PROCESS DATA ---------------------------------------------------------------

se.list <- lapply(unique(infoTable$patientid), function(s) {
  print(s)
  se <- InputFromTable(infoTable[infoTable$patientid==s,,drop=F])
  se <- SCTransform(se) %>%
    LoadImages(time.resolve = FALSE) %>%
    RunNMF(nfactors = 20) %>% 
    RunUMAP(reduction = "NMF", dims = 1:20) %>%
    FindNeighbors(reduction = "NMF", dims = 1:20) %>% 
    FindClusters()
})

# rename obj files
names(se.list)[1] <- "1142243F"
names(se.list)[2] <- "1160920F"
names(se.list)[3] <- "CID4290"
names(se.list)[4] <- "CID4465"
names(se.list)[5] <- "CID44971"
names(se.list)[6] <- "CID4535"

for(sample in unique(infoTable$patientid)){
  print(sample)
  print(head(se.list[[sample]]@meta.data))
}

# 04: ADD PATH ANNOTATIONS ----------------------------------------------------

temp_path <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/spatial/Garvan_PhaseII_spaceranger_loupe_annotations/"
temp_csv_combined <- NULL
for(sample in unique(infoTable$patientid)){
  if(sample == "1142243F"){
    temp_csv <- read.csv(paste0(temp_path,"182964_1142243F_new.csv"))
  }
  if(sample == "1160920F"){
    temp_csv <- read.csv(paste0(temp_path,"182968_1160920F_new.csv"))
  }
  if(sample == "CID4290"){
    temp_csv <- read.csv(paste0(temp_path,"35_E2_4290_new.csv"))
  }
  if(sample == "CID4465"){
    temp_csv <- read.csv(paste0(temp_path,"35_C2_4465_new.csv"))
  }
  if(sample == "CID44971"){
    temp_csv <- read.csv(paste0(temp_path,"35_B1_4497-1_new.csv"))
  }
  if(sample == "CID4535"){
    temp_csv <- read.csv(paste0(temp_path,"35_D1_4535_new.csv"))
  }

  colnames(temp_csv) <- c("Barcode", "Classification")
  temp_csv$patientid <- paste0(unique(se.list[[sample]]@meta.data$patientid))
  # print(temp_sampleID)
  temp_csv_combined <- rbind(temp_csv_combined, temp_csv)
}

levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "(Need to check) folding artefact"] <- "Artefact"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "(Need to check) invasive cancer + lymphocytes"] <- "Invasive cancer + lymphocytes"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "(Need to check) lymphocyte aggregation"] <- "Lymphocytes"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "stroma"] <- "Stroma"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == " stroma"] <- "Stroma"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "invasive cancer + stroma + lymphocytes"] <- "Invasive cancer + stroma + lymphocytes"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "stroma + adipose tissue"] <- "Stroma + adipose tissue"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "Lymphocyte aggregation"] <- "Lymphocytes"
levels(temp_csv_combined[,2])[levels(temp_csv_combined[,2]) == "Normal + stroma + Lymphocytes"] <- "Normal + stroma + lymphocytes"

print(sort(as.vector(unique(temp_csv_combined[,2]))))

for(sample in unique(infoTable$patientid)){
  print(sample)
  # print(head(se.list[[sample]]@meta.data))
  
  temp_csv_combined_subset <- temp_csv_combined[temp_csv_combined$patientid == sample,,drop=F]
  print(head(temp_csv_combined_subset))
}

# append path annotations for CTP samples
for(sample in unique(infoTable$patientid)){
  print(sample)
  temp_csv_combined_subset <- temp_csv_combined[temp_csv_combined$patientid == sample,,drop=F]
  rownames(temp_csv_combined_subset) <- temp_csv_combined_subset$Barcode
  # temp_newid <- sample-6
  rownames(temp_csv_combined_subset) <- gsub(paste0("-1"),
                                             paste0("-1_1"),
                                             rownames(temp_csv_combined_subset))
  colnames(temp_csv_combined_subset) <- c("Barcode", "Classification", "sample")
  
  if(nrow(temp_csv_combined_subset) == nrow(se.list[[sample]]@meta.data)){
    print("rows match")
    temp_csv_combined_subset <- temp_csv_combined_subset[rownames(se.list[[sample]]@meta.data),,drop=F]
  } else {
    print("rows dont match > filtering")
    temp_csv_combined_subset <- temp_csv_combined_subset[rownames(temp_csv_combined_subset) %in% rownames(se.list[[sample]]@meta.data),,drop=F]
  }
  
  if(nrow(temp_csv_combined_subset) == nrow(se.list[[sample]]@meta.data)){
    print("rows now match")
    temp_csv_combined_subset <- temp_csv_combined_subset[rownames(se.list[[sample]]@meta.data),,drop=F]
  } else {
    print("missing barcodes > appending as NA")
    
    temp_missing_barcodes <- rownames(se.list[[sample]]@meta.data)[! rownames(se.list[[sample]]@meta.data) %in% rownames(temp_csv_combined_subset)]
    temp_csv_combined_subset_append <- data.frame(row.names = temp_missing_barcodes,
                                                  Barcode = temp_missing_barcodes,
                                                  Classification = "NA")
    temp_csv_combined_subset_append$sample <- sample
    temp_csv_combined_subset <- rbind(temp_csv_combined_subset,
                                      temp_csv_combined_subset_append) 
  }
  
  if(nrow(temp_csv_combined_subset) == nrow(se.list[[sample]]@meta.data)){
    print("rows now match")
    temp_csv_combined_subset <- temp_csv_combined_subset[rownames(se.list[[sample]]@meta.data),,drop=F]
  } else {
    print("rows still dont match")
  }
  
  if(all.equal(rownames(temp_csv_combined_subset),
               rownames(se.list[[sample]]@meta.data))){
    print("vectors align")
    se.list[[sample]]@meta.data$Classification <- temp_csv_combined_subset[,2]
    # print(table(se.list[[sample]]@meta.data$Classification))
  } else {
    print("vectors dont align")
  }
  
  
  
}


# plot new annotations
dir.create("output/")
dir.create("output/Path_annotations")

for(sample in unique(infoTable$patientid)){
  temp_sampleID <- paste0(unique(se.list[[sample]]@meta.data$patientid),
                          "_rep_",
                          unique(se.list[[sample]]@meta.data$rep))
  
  temp_ggplot <- FeatureOverlay(se.list[[sample]],
                                features = "Classification", 
                                # pt.size = 0.75,
                                # pt.alpha = 0.75,
                                sample.label = F,
                                type = "raw")
  # n <- paste0("temp_ggplot_",sample)
  # assign(n, temp_ggplot)
  
  # temp_png_function <-
  #   function(x) {
  #     png(
  #       file = (x),
  #       width = 8,
  #       height = 6,
  #       res = 600,
  #       units = 'in'
  #     )
  #   }
  # 
  # temp_png_function(paste0("output/08_path_annotations_July_SandyTony/01_FeatureOverlay_Classification_",temp_sampleID,".png"))
  # print(temp_ggplot)
  # dev.off()
  
  temp_pdf_function <-
    function(x) {
      pdf(file = (x),
          width = 8, 
          height = 5)
    }
  temp_pdf_function(paste0("output/Path_annotations/01_FeatureOverlay_Classification_",temp_sampleID,".pdf"))
  print(temp_ggplot)
  dev.off()
  
}





# SAVE RDS ----------------------------------------------------------------

saveRDS(se.list, "RDATA_Visium_brca_objects.Rdata")

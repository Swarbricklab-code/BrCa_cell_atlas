# SURVIVAL ANALYSIS - MACROPHAGES
# SUNNY Z WU
# CODE MODIFIED AND ADOPTED FORM CHENFEI WANG
#
# qrsh -pe smp 8 -l mem_requested=10G -P TumourProgression
# source activate r_seuratdev
# R
#
#
# 
#
# 01: PARSE ARGUMENTS ---------------------------------------------------------

temp_args <-
  commandArgs(trailingOnly = T)

# Number of genes to use per signature - 100
temp_number  <- 
  as.numeric(temp_args[1])

temp_if_run_within_celltype_signatures <- T

# 02: SETUP  ---------------------------------------------------------

library(survival)
library(survminer)
library(RColorBrewer)
library(org.Hs.eg.db)
library(dplyr)

# colour scaled
ncol1 = brewer.pal(12,"Paired")
ncol2 = brewer.pal(8,"Dark2")


# directory on cluster
setwd("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/")
dir.create("METABRIC_macrophage_survival")
setwd("METABRIC_macrophage_survival")


# 03: LOAD METABRIC COHORT DATA ---------------------------------------------------------------

temp_dir <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/cohort_data/"

## S1, read the expression table
expr_METABRIC = as.matrix(read.delim(paste0(temp_dir,"METABRIC.norm_subtract")))     # METABRIC
rownames(expr_METABRIC) = as.vector(unlist(mget(rownames(expr_METABRIC), envir=org.Hs.egSYMBOL, ifnotfound=NA)))


# 04: LOAD GENE SIGNATURES ----------------------------------------------------

temp_dir <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/gene_signatures/"
## S2, classify the patient based on averaged expresion of marker genes from each cell-type
# use within cell type derived cell signatures
# within cell type gene signatures (e.g. iCAFs vs all CAFs)
gene_46 <- NULL
celltype <- "Myeloid"
gene = read.csv(paste0(temp_dir,"within_celltypes_subset/",celltype,"/01_MAST_FINDALLMARKERS_RNA_annotated_celltype.csv"),
                    row.names=1)
gene = gene[gene$p_val_adj < 0.01,]
gene_46 <- rbind(gene_46,gene)
gene_46 <- gene_46[gene_46$cluster %in% c("Myeloid_Macrophage_LAM 2", "Myeloid_Macrophage_LAM 1"),,drop=F]

    
    
# 05: APPEND METABRIC SUBTYPE INFO -----------------------------------------------------

temp_dir <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/cohort_data/"

## S3, add the subtype information to the gene expression matrix
clin_METABRIC = read.delim(paste0(temp_dir,"METABRIC.clinical"))

determine_subtype = function(clin)
{
  clin$subtype = NA
  clin[which(clin$LumA==1),"subtype"] = "LumA"
  clin[which(clin$LumB==1),"subtype"] = "LumB"
  clin[which(clin$Her2==1),"subtype"] = "Her2"
  clin[which(clin$Basal==1),"subtype"] = "Basal"
  return(clin)
}
clin_METABRIC = determine_subtype(clin_METABRIC)

# ADD DISCOVERY VS VALIDATION COHORT AND APPEND MET 

# taken from original study from Curtis et al. https://www.nature.com/articles/nature10983#Sec10
# Table 2: Clinical annotation for the discovery cohort of 997 breast cancer patients
# Table 3: Clinical annotation for the validation cohort of 995 breast cancer patients

temp_dir <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/cohort_data/"
temp_discovery <- read.delim(paste0(temp_dir,"METABRIC_table_S2_revised.txt"), row.names = "METABRIC_ID")
temp_discovery$set <- "discovery"
temp_validation <- read.delim(paste0(temp_dir,"METABRIC_table_S3_revised.txt"),
                              row.names = "METABRIC_ID")
temp_validation$set <- "validation"
temp_combined <- rbind(temp_discovery, 
                       temp_validation)
rownames(temp_combined) <- gsub("-","_",rownames(temp_combined))
temp_combined_filtered <- temp_combined[rownames(temp_combined) %in% rownames(clin_METABRIC),,drop=F]
temp_combined_filtered <- temp_combined_filtered[rownames(clin_METABRIC),,drop=F]
# append to metadata
if(all.equal(rownames(temp_combined_filtered),
             rownames(clin_METABRIC))){
  print('matching')
  clin_METABRIC$set <- temp_combined_filtered$set
}
clin_METABRIC_dis <- clin_METABRIC[clin_METABRIC$set == "discovery",,drop=F]
clin_METABRIC_val <- clin_METABRIC[clin_METABRIC$set == "validation",,drop=F]
expr_METABRIC_dis <- expr_METABRIC[,colnames(expr_METABRIC) %in% rownames(clin_METABRIC_dis),drop=F]
expr_METABRIC_val <- expr_METABRIC[,colnames(expr_METABRIC) %in% rownames(clin_METABRIC_val),drop=F]


# 06: GENERATE AVERAGE EXPRESSION MATRICES --------------------------------------

dir.create("genes_used_for_prognosis")

averaged_gene_signature <- function(inmatrix, gene, top = temp_number, cohort)
{
  outmatrix = NULL
  name = NULL
  
  temp_cluster_ids <- as.vector(unique(gene$cluster))
  # remove IGs from metabric
  if(ncol(inmatrix) == nrow(clin_METABRIC_dis) | ncol(inmatrix) == nrow(clin_METABRIC_val)){
    temp_ig_remove <- grep("^IGH",gene$gene,value=T)
    temp_ig_remove <- c(temp_ig_remove, grep("^IGK",gene$gene,value=T))
    temp_ig_remove <- c(temp_ig_remove, grep("^IGL",gene$gene,value=T))
    gene <- gene[! gene$gene %in% temp_ig_remove,,drop=F]
  }
  
  # filter for shared genes between platforms
  gene <- gene[gene$gene %in% rownames(inmatrix),,drop=F]
  temp_df_combined <- data.frame(row.names = c(1:top),
                                 gene = c(1:top))
  for(i in temp_cluster_ids)
  {
    temp_length_intersect <- length(intersect(rownames(inmatrix), 
                                              as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))))
    # print(i)
    # print(temp_length_intersect)
    
    outmatrix = cbind(outmatrix, 
                      apply(inmatrix[intersect(rownames(inmatrix), 
                                               as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))),],
                            2,
                            mean))
    name = c(name, gsub("-","",gsub("\\+","",gsub("\\ ","",as.character(i)))))
    
    temp_df <- data.frame(cluster = intersect(rownames(inmatrix), 
                                              as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))))
    colnames(temp_df) <- i
    
    write.csv(temp_df, 
              paste0("genes_used_for_prognosis/temp_df_combined_",cohort,"_",i,".csv"))
    
  }
  
  colnames(outmatrix) = name
  return(outmatrix)
}

expr_METABRIC_dis_mean_46 <- averaged_gene_signature(expr_METABRIC_dis, gene_46, cohort= "METABRIC_dis")
expr_METABRIC_val_mean_46 <- averaged_gene_signature(expr_METABRIC_val, gene_46, cohort= "METABRIC_val")

expr_METABRIC_dis_mean_46 = data.frame(clin_METABRIC_dis[rownames(expr_METABRIC_dis_mean_46),], expr_METABRIC_dis_mean_46)
expr_METABRIC_val_mean_46 = data.frame(clin_METABRIC_val[rownames(expr_METABRIC_val_mean_46),], expr_METABRIC_val_mean_46)



# 07: RUN ITERATIVE SURVIVAL ANALYSIS -----------------------------------------

dir.create("Clinical_Mean_48celltype_ALL")

## S4, survival analysis
survival_analysis_celltype <- function(inmatrix, name, celltype, subtype, signature) 
{
  
  # M1, stratify by top 30% and bottom 30%
  inmatrix = inmatrix[order(inmatrix[,celltype], decreasing=T),]
  inmatrix_top = inmatrix[1:(nrow(inmatrix)%/%3),];inmatrix_top$CELLTYPE = 1
  inmatrix_bottom = inmatrix[(nrow(inmatrix)%/%3*2):nrow(inmatrix),];inmatrix_bottom$CELLTYPE = 0
  inmatrix = rbind(inmatrix_top,inmatrix_bottom)
  colnames(inmatrix) <- toupper(colnames(inmatrix))
  surv <- survfit(Surv(OS, EVENT) ~ CELLTYPE, data = inmatrix)
  
  # make p-values plotted (which is from a log-rank test) = to p-value derived from coxprop model
  coxph.surv=Surv(inmatrix[,'OS'],inmatrix[,'EVENT'])
  temp_summary <- as.data.frame(summary(coxph(coxph.surv~inmatrix[,'CELLTYPE']+inmatrix[,'AGE']))$coefficients[1,c(2,5)])
  
  
  ggsurvplot(surv, data = inmatrix,  title = paste(name,celltype), font.title = c(10, "bold", "darkblue"), palette = ncol1, pval = round(temp_summary[2,],digits=8), legend="top", legend.labs = c("bottom 30%","top 30%"), risk.table = F)
  ggsave(paste0("Clinical_Mean_",signature,"celltype_",subtype,"/Clinical_Mean_",celltype,"_",name,".pdf"), width=5,height=5)
  
  return(summary(coxph(coxph.surv~inmatrix[,'CELLTYPE']+inmatrix[,'AGE']))$coefficients[1,c(2,5)])
}

result = NULL
for(celltype in c("Myeloid_Macrophage_LAM1","Myeloid_Macrophage_LAM2"))
{
  result <- rbind(result, survival_analysis_celltype(expr_METABRIC_dis_mean_46, name= "METABRIC_dis", celltype = celltype, "ALL", signature= 48))
  
  result <- rbind(result, survival_analysis_celltype(expr_METABRIC_val_mean_46, name="METABRIC_val", celltype,  subtype ="ALL", signature = 48))
}
result <- as.data.frame(result)
colnames(result) <- c("HR", "Pvalue")
result$celltype <- c("Myeloid_Macrophage_LAM1","Myeloid_Macrophage_LAM2",
                     "Myeloid_Macrophage_LAM1","Myeloid_Macrophage_LAM2")
result$cohort <- c("METABRIC_dis","METABRIC_val",
                     "METABRIC_dis","METABRIC_val")

write.csv(result,
          "coxph.surv_data_summary.csv")

# print number of patients per group
for(cohort in c("expr_METABRIC_dis_mean_46", "expr_METABRIC_val_mean_46")){
  print(cohort)
  inmatrix <- get(paste0(cohort))
  
  # M1, stratify by top 30% and bottom 30%
  inmatrix = inmatrix[order(inmatrix[,celltype], decreasing=T),]
  inmatrix_top = inmatrix[1:(nrow(inmatrix)%/%3),];inmatrix_top$CELLTYPE = 1
  print(nrow(inmatrix_top))
  
  inmatrix_bottom = inmatrix[(nrow(inmatrix)%/%3*2):nrow(inmatrix),];inmatrix_bottom$CELLTYPE = 0
  print(nrow(inmatrix_bottom))
  
  inmatrix = rbind(inmatrix_top,inmatrix_bottom)
  colnames(inmatrix) <- toupper(colnames(inmatrix))
  
}

# [1] "expr_METABRIC_dis_mean_46"
# [1] 310
# [1] 311
# [1] "expr_METABRIC_val_mean_46"
# [1] 177
# [1] 180




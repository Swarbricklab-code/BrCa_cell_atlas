###Processing of bulk RNAseq datasets for Australian set of samples and clustering with TCGA samples
##Quantile normalizing the raw data matrix
#1.This function does QUARTILE Normalization on STAR-SALMON quantified raw gene-data matrix. Change .75 argument to any ything to normalize.
quartileNorm<- function(x,y=NA){
  uqs<- apply(x,2,function(x){ quantile(x[x>0 & !is.na(x)],0.75)})
  if(is.na(y)){
    y<- median(uqs)
  }
  x.norm <- t(apply(x,1,function(x,y){x*y},y/uqs))
  dimnames(x.norm)<- dimnames(x)
  return(x.norm)
}

x=as.matrix((read.csv('Aus-mini-Atlassalmon.matrixforUQ',header =F)))
x.norm <- quartileNorm(x)
write.csv(x=x.norm, file = "quartileNorm_RibozeroAusmatrix.csv")

##Upload the quantile normalized file on Cluster 3.0 for log transformation and median centering
#2.For log transforming - Apply log transform data in the 'Adjust data' tab and save the data matrix
#3.For median centering - Filter data by clicking 80% present genes and at least 1 obs with abs(val) >= 2.0 in the 'Filter data' tab. Accept filter and then apply median centering on the filtered set of genes in the 'Adjust data' tab. Save the data matrix.

#4.Column Standardizing the data matrix for combining with the TCGA data matrix
standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}

x<-read.table("quartileNorm_RibozeroAusMatrix_logtransformedwithOriginalIDsforclustering_medcent.txt", sep="\t", header = TRUE, row.names = 1)
xstd<-standardize(x)

## Combined with the TCGA data matrix which has been upper-quartile normalized, log transformed, filtered, gene median centred and column standardized as above
## Use the Perou Lab Intrinisic gene list 
Combineddata<-merge(TCGA, x, by="row.names", sort = FALSE)

# Cluster this combined matrix in Cluster 3.0 using Correlation centered metric for both genes and arrays in the 'Heirarchical' tab
# View cluster files in Java Treeview



###Running PAM50 on bulk RNAseq datasets using the Zhao ER corrrection


require(Biobase)
require(plyr)

source("./arrayTools_MBSedits_collapseID.R")
load("./ssBCsubtyping.Rdata")

##The original Zhao correction is only for ER; Therefore we perform the following correction for Her2 to run PAM50 on Her2+ samples
#column HER2 to UNC232 phenoData
setwd("./Zhao_ER_HER2_correction/")
x <- read.delim("merge1.txt", header = TRUE)
a <- cbind(UNC232@phenoData@data, x)
a
UNC232@phenoData@data <- a
UNC232@phenoData@data 

#HER2 label to varMetadata
b <- read.delim("metadata.txt", header = TRUE, row.names = 1)
b
UNC232@phenoData@varMetadata <- b
UNC232@phenoData@varMetadata

#only for HER2 positive disease: rerun to estimated subgroup specific quantile by:
setwd("./data/")
gene.sigma = getSigma.generic(exprs(UNC232), pData(UNC232))
Pam50 = signature.pam50(gene.sigma = gene.sigma)

##Performing PAM50 subtyping on the Australian bulk RNA-seq samples
##Use the upper-quartile normalized and log transformed data matrix generated from Step 2. above
##Separate samples into ER+, Her2+ and TNBC based on their clinical IHC information. Separate dual ER+ and Her2+ separately from only ER+ and only Her2+
##Only use the PAM50 50 genes in the input gene-matrix

#ERpos_HER2neg
x <- read.delim("PAM50_logtransformedwithOriginalIDsforclustering_ERpos.txt", header = TRUE, row.names = 1)
subtypes = pam50.symbol2symbol.v2(x, s = "ERpos_HER2neg")
write.table(subtypes, "PAM50calls_Swarbrick_ERpos_HER2neg.txt", sep = "\t", col.names = NA)
#Erneg_Her2neg
x <- read.delim("PAM50_logtransformedwithOriginalIDsforclustering_TNBC.txt", header = TRUE, row.names = 1)
subtypes = pam50.symbol2symbol.v2(x, s = "ERneg_HER2neg")
write.table(subtypes, "PAM50calls_Swarbrick_TNBC.txt", sep = "\t", col.names = NA)
#HErpos_ERpos
x <- read.delim("PAM50_logtransformedwithOriginalIDsforclustering_Her2Erpos.txt", header = TRUE, row.names = 1)
subtypes = pam50.symbol2symbol.v2(x, s = "HER2pos_ERpos")
write.table(subtypes, "PAM50calls_Swarbrick_HER2pos_ERPos.txt", sep = "\t", col.names = NA)
#HER2pos_ERneg
x <- read.delim("PAM50_logtransformedwithOriginalIDsforclustering_Her2pos.txt", header = TRUE, row.names = 1)
subtypes = pam50.symbol2symbol.v2(x, s = "HER2pos_ERneg")
write.table(subtypes, "PAM50calls_Swarbrick_HER2pos_ERneg.txt", sep = "\t", col.names = NA)
#Merge all the result files

#The aim of the script is transfer of CITE data from known to unknown samples
#with standard Seurat protocol: https://satijalab.org/seurat/archive/v3.0/integration.html
#and within specific cell types only.

library(Seurat)
library(dplyr)
library(Matrix)
library(sctransform)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(plyr)


#directory structure
homeDir="../"
RobjectsDir=paste0(homeDir,"project_results/Robjects/")
#figDir<-paste0(homeDir,"project_results/figures/")
#system(paste("mkdir -p",figDir))


#transfer of CITE information will be performed for specific cell types
annotations<-c("B_cells","Myeloid","STROMAL","T_cells")
cellTypesBarcodes<-list()

#get a list of all barcodes used
for(cellType in annotations){
	cat(cellType)
	cat("\n")
	integratedObject<-paste0(RobjectsDir,"RDATA_02_PROCESSED_FILTERED_",cellType,".Rdata")
	citeCells<-readRDS(integratedObject)
	cellTypesBarcodes[[cellType]]<-colnames(citeCells)
}

#fetch the list of sample names from the barcodes
BCs<-unlist(cellTypesBarcodes)
sampleNames<-unique(gsub("_.*","",BCs))

#select sample names for samples with CITE information
sampleNames_CITE<-sampleNames[grepl("3838|3946|4040|4378|4515",sampleNames)]

#create individual whitelists to be used for processing of raw CITE data
for(sampleName in sampleNames_CITE){
	cat(sampleName)
	cat("\n")
	#get all barcodes
	BCs_CITE<-BCs[grepl(sampleName,BCs)]
	BCs_CITE<-gsub(".*_","",BCs_CITE)
	#export them
	shortName<-gsub("CID","",sampleName)
	outWhitelist<-paste0(homeDir,"raw_files/",shortName,"/CITE_fastqs/",shortName,"_whitelist_miniatlas.txt")
	write(bcs,outWhitelist)
}


######## Take only those CITE ADTs from individual samples that can distinguish biological clusters
ABsList<-list()
functionalCITE<-list()
rawCiteCounts<-list()

for(sampleName in gsub("CID","",sampleNames_CITE)){
	cat(sampleName)
	cat("\n")
	#for each sample extract Seurat-normalised CITE
	#data from previously performed independent analyses
	load(paste0(RobjectsDir,sampleName,"_v3_CITE.Rdata"))

	#decide which antibody will be used by loading the Robjects and selecting 
	adt.markers <- FindAllMarkers(cells, assay = "ADT", only.pos = TRUE,verbose=F)
	adt.markers <- adt.markers[!grepl("unmapped",adt.markers$gene),]
 	functional_ADTs <- unique(adt.markers$gene)
 	cat("Total marker antibodies:\t")
	cat(length(unique(functional_ADTs)))
	cat("\n")
	ABsList[[sampleName]]<-functional_ADTs

	#subselect only the functional data
	citeRaw<-cells@assays$ADT@counts
	colnames(citeRaw)<-paste0("CID",sampleName,"_",colnames(citeRaw))
	functionalCITE[[sampleName]]<-citeRaw[row.names(citeRaw) %in% functional_ADTs,]
}

#combine data across all samples
functionalCITEMelt<-lapply(functionalCITE,function(x){melt(as.matrix(x))})
functionalCITENarrow<-do.call(rbind,functionalCITEMelt)
functionalCITEWide<-acast(functionalCITENarrow, Var1 ~ Var2)
functionalCITEWide[is.na(functionalCITEWide)]=0

######## Perform data transfer
annotations<-c("T_cells","B_cells","Myeloid","STROMAL")
imputedSeurat4celltypes<-list()


imputeCITE<-function(cells,cellType){
	cite<-functionalCITEWide[,colnames(functionalCITEWide) %in% colnames(cells)]
	cite<-as.data.frame(cite)
	cite$id<-row.names(cite)
	allBCs<-colnames(cells@assays$RNA@data)

	#fix missing zeroes so that all cells and all ADTs have information
	citeAll<-as.data.frame(matrix(0,dim(cite)[1],length(allBCs)))
	row.names(citeAll)<-row.names(cite)
	colnames(citeAll)<-allBCs
	citeAll$id<-row.names(citeAll)
	citeL<-list(cite=cite,citeAll=citeAll)
	matrix.df <- ldply(citeL, melt)
	#some values are NA, so turn them into 0
	matrix.df$value[is.na(matrix.df$value)]<-0
	sum.matrix <- acast(matrix.df, id ~ variable, sum)
	citeAll<-sum.matrix
	citeAll<-as.data.frame(citeAll)
	citeAll<-citeAll[colnames(citeAll) %in% colnames(cells),]

	#perform CLR analysis
	cells[["ADT"]] <- CreateAssayObject(counts = citeAll)
	cells <- NormalizeData(cells, assay = "ADT", normalization.method = "CLR")
	cells <- ScaleData(object = cells, assay = "ADT")

	#perform data transfer
	allCells<-colnames(cells@assays$RNA@data)
	nonQuery<-allCells[!(allCells %in% colnames(cite))]
	reference <- subset(cells,cells=colnames(cite) )
	query <- subset(cells,cells=nonQuery )

	query.anchors <- FindTransferAnchors(reference = reference, query = query, project.query=T)
	predictions1 <- TransferData(anchorset = query.anchors, refdata = reference@assays$ADT@counts)
	predictionsDF<-as.data.frame(predictions1@data)
	cite<-cite[order(row.names(cite)),]
	predictionsDF<-predictionsDF[row.names(predictionsDF) %in% row.names(cite),]
	predictionsDF<-predictionsDF[order(row.names(predictionsDF)),]

	cite<-cbind(cite,predictionsDF)
	cite<-cite[,colnames(cite) %in% colnames(cells@assays$RNA@data)]

	#perform CLR on the transferred data and return the object
	cells[["ADT"]] <- CreateAssayObject(counts = cite)
	cells <- NormalizeData(cells, assay = "ADT", normalization.method = "CLR")
	cells <- ScaleData(object = cells, assay = "ADT")
}


#perform data transfer for each cell type
for(cellType in annotations){
	cat(cellType)
	cat("\n")
	integratedObject<-paste0(RobjectsDir,"RDATA_02_PROCESSED_FILTERED_",cellType,".Rdata")
	citeCells<-readRDS(integratedObject)
	for(norm in c(FALSE)){
		cat(norm)
		cat("\n")
		cells<-imputeCITE(citeCells,cellType)
		save(cells,file=paste0(RobjectsDir,"imputed_",cellType,".Rdata"))
	}
}



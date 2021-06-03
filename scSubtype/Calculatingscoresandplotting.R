# Calculating sc50 scores
#Aatish Thennavan Perou Lab

# load packages
library(Seurat)
library(cowplot)

#Reading in the datasets
TNBCtraining<-read.table("TNBCmerged_Training_SwarbrickInferCNV.txt", sep="\t", header = TRUE)
LumAtraining<-read.table("LumAmerged_Training_SwarbrickInferCNV.txt", sep="\t", header = TRUE)
LumBtraining<-read.table("LumBmerged_Training_SwarbrickInferCNV.txt", sep="\t", header = TRUE)
Her2training<-read.table("Her2merged_Training_SwarbrickInferCNV.txt", sep="\t", header = TRUE)
Testingdataset<-read.table("Testingdatasetsmerged_SwarbrickInferCNV.txt", sep="\t", header = TRUE)

#Seurat CCA integration
Group1 <- CreateSeuratObject(counts = LumAtraining, project = "Training_LumA", min.cells = 5)
Group1$stim <- "Group1"
Group1 <- subset(Group1, subset = nFeature_RNA > 500)
Group1 <- NormalizeData(Group1, verbose = FALSE)
Group1 <- FindVariableFeatures(Group1, selection.method = "vst", nfeatures = 3000)
Group2 <- CreateSeuratObject(counts = TNBCtraining, project = "Training_TNBC", min.cells = 5)
Group2$stim <- "Group2"
Group2 <- subset(Group2, subset = nFeature_RNA > 500)
Group2 <- NormalizeData(Group2, verbose = FALSE)
Group2 <- FindVariableFeatures(Group2, selection.method = "vst", nfeatures = 3000)
Group3 <- CreateSeuratObject(counts = Her2training, project = "Training_Her2", min.cells = 5)
Group3$stim <- "Group3"
Group3 <- subset(Group3, subset = nFeature_RNA > 500)
Group3 <- NormalizeData(Group3, verbose = FALSE)
Group3 <- FindVariableFeatures(Group3, selection.method = "vst", nfeatures = 3000)
Group4 <- CreateSeuratObject(counts = LumBtraining, project = "Training_LumB", min.cells = 5)
Group4$stim <- "Group4"
Group4 <- subset(Group4, subset = nFeature_RNA > 500)
Group4 <- NormalizeData(Group4, verbose = FALSE)
Group4 <- FindVariableFeatures(Group4, selection.method = "vst", nfeatures = 3000)
Group5 <- CreateSeuratObject(counts = Testingdataset, project = "Testing", min.cells = 5)
Group5$stim <- "Group5"
Group5 <- subset(Group5, subset = nFeature_RNA > 500)
Group5 <- NormalizeData(Group5, verbose = FALSE)
Group5 <- FindVariableFeatures(Group5, selection.method = "vst", nfeatures = 3000)
Tumor.anchors <- FindIntegrationAnchors(object.list = list(Group1, Group2, Group3, Group4, Group5), dims = 1:50)
Tumor.combined <- IntegrateData(anchorset = Tumor.anchors, dims = 1:50)

#Calculating SC50 scores on the 'RNA' scaled data
DefaultAssay(Tumor.combined) <- "RNA"
Tumor.combined[["percent.mt"]] <- PercentageFeatureSet(object = Tumor.combined, pattern = "^MT.")
plot1 <- FeatureScatter(object = Tumor.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Tumor.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Tumor.combined <- subset(x = Tumor.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
Tumor.combined <- NormalizeData(object = Tumor.combined)
Tumor.combined <- FindVariableFeatures(object = Tumor.combined, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(x = Tumor.combined)
Tumor.combined <- ScaleData(object = Tumor.combined, features = all.genes)
Tumor.combined <- RunPCA(object = Tumor.combined, features = VariableFeatures(object = Tumor.combined))
sigdat <- read.table("SinglecellMolecularSubtypesignaturesonlytumor_SEPTEMBER2019.txt",sep='\t',header=F,row.names=1,fill=T)
tocalc<-as.data.frame(Tumor.combined@assays$RNA@scale.data)
outdat <- matrix(0,nrow=nrow(sigdat),
ncol=ncol(tocalc),
dimnames=list(rownames(sigdat),
colnames(tocalc)))
for(i in 1:nrow(sigdat)){
sigdat[i,!is.na(sigdat[i,])]->module
row <- as.character(unlist(module))
row<-unique(row[row != ""])
genes<-which(rownames(tocalc) %in% row)
temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
outdat[i,]<-as.numeric(temp)
}
final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)
Tumor.combined@assays$RNA@data <- rbind(Tumor.combined@assays$RNA@data,finalm)

#Identifying clusters in the 'Integrated' data
DefaultAssay(object = Tumor.combined) <- "integrated"
Tumor.combined <- ScaleData(Tumor.combined, verbose = FALSE)
Tumor.combined <- RunPCA(Tumor.combined, npcs = 50, verbose = FALSE)
Tumor.combined <- FindNeighbors(object = Tumor.combined, dims = 1:30)
Tumor.combined <- FindClusters(object = Tumor.combined, resolution = 0.2)
Tumor.combined <- RunUMAP(Tumor.combined, reduction = "pca", dims = 1:30)

#Plotting Functions of signatures on the identified clusters
library(RColorBrewer)
colourCount = length(unique(Tumor.combined$orig.ident))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p1 <- DimPlot(object = Tumor.combined, reduction = "umap", group.by = "orig.ident", pt.size = 2.5, cols = getPalette(colourCount))
p2<-FeaturePlot(object = Tumor.combined, features = c("LumA_SC50"), pt.size = 2.5, min.cutoff = 0)
plot_grid(p1, p2)
p2<-FeaturePlot(object = Tumor.combined, features = c("LumB_SC50"), pt.size = 2.5, min.cutoff = 0)
plot_grid(p1, p2)
p2<-FeaturePlot(object = Tumor.combined, features = c("Her2_SC50"), pt.size = 2.5, min.cutoff = 0)
plot_grid(p1, p2)

#First source the code and set all the files needed in the working directory
source('DifferentiationPredictor_functions.R')
source('arrayTools_MBSedits_collapseID.R')

#Read in the files
limModelFile<- "differentiationCentroids_LimDWD_Genes.txt"
lim<-readarray(limModelFile,hr=1)
Singlecell<-readarray("Swarbrick_Tumorcells_Top2000expressedgenes.txt",impute=F, hr=1)

#Calculate Differentiation Score
diffScore<-assignDiffScore.dwd(Singlecell,lim)

#Calculate Proliferation Score
sigdat<-read.table("Proliferation signature_UNC_Genes symbol.txt",sep="\t",header=F,row.names=1,fill=T)
tocalc<-read.table("Swarbrick_Tumorcells_Top2000expressedgenes.txt",sep="\t",header=T, row.names=1)

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
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 3)
finalm<-as.matrix(final) 

finalmt<-t(finalm)

DS_PS<-cbind(finalmt,diffScore)
write.table(DS_PS, "DS_PS_Swarbrick_Alltumorcells.txt", sep="\t")
##Code for Signature scores and highest Calls
#Aatish Thennavan Perou Lab

sigdat <- read.table("SinglecellMolecularSubtypesignaturesincludingSWARBRICKnormals_SEPTEMBER2019.txt",sep='\t',header=F,row.names=1,fill=T)


#Read in the single cell RDS object as 'Mydata'
tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)



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

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

##Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]

##Writing out output files (rownames remain the same for both)
write.table(finalm.sweep.t, "Mydata_Scores.txt", sep="\t")
write.table(Finalnames, "Mydata_CALLS.txt", sep="\t")

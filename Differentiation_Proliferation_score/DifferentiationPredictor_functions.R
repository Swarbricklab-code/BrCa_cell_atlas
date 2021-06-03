
collapseIDs<-function(x,allids=row.names(x),method="mean"){

        allids<-as.vector(allids)
        ids<- levels(as.factor(allids))
        x.col<- NULL

        if(length(ids)==dim(x)[1]){ 
                        dimnames(x)[[1]]<-allids
                        return(x) 
        }
        
        for(i in 1:length(ids)){
                if(sum(allids==ids[i])>1){
                        indices <- allids==ids[i] 
                        if(method=="mean"){
                                vals<-apply(x[indices,],2,mean,na.rm=T)
                        }
                        if(method=="median"){
                                vals<-apply(x[indices,],2,median,na.rm=T)
                        }
                        if(method=="stdev"){   
                                temp<- x[indices,]
                                stdevs<- apply(temp,1,sd,na.rm=T)
                                vals<- temp[match(max(stdevs),stdevs),]
                        }
                        if(method=="iqr"){   
                                temp<- x[indices,]
                                iqrs<- apply(temp,1,function(x){quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T)})
                                vals<- temp[match(max(iqrs),iqrs),]
                        }
                        x.col <- rbind(x.col,vals)
                }else{
                        x.col <- rbind(x.col,x[allids==ids[i],])
                }
        }

        dimnames(x.col)<- list(ids,dimnames(x)[[2]])
        return(x.col)
        
}


readarray<-function(dataFile,designFile=NA,hr=1,impute=T,method="mean"){

	headerRows <- hr

	x<-read.table(dataFile,sep="\t",header=F,fill=T,stringsAsFactors=FALSE)

	if(headerRows==1){
			sampleNames<-as.vector(t(x[1,-1]))
			x<-x[-1,]
			classes<-NULL
			ids<-x[,1]
			xd<-x[,-1]
			xd<-apply(xd,2,as.numeric)
			xd<-collapseIDs(xd,ids,method)	
	}else{
			sampleNames<-as.vector(t(x[1,-1]))
			x<-x[-1,]
			
			classes<-x[1:(headerRows-1),]
			dimnames(classes)[[1]]<-classes[,1]
			classes<-classes[,-1]
			classes[classes==""]<-NA
			classes<-t(classes)
			rownames(classes)<-sampleNames
			classes<-as.data.frame(classes)
						
			xd<-x[(-1:-(headerRows-1)),]
			ids<-as.vector(t(xd[,1]))
			xd<-xd[,-1]
			xd<-apply(xd,2,as.numeric)
			xd<-collapseIDs(xd,ids,method)
	}
	
	features<- dim(xd)[1]
	samples<- dim(xd)[2]
	geneNames<-rownames(xd)
	xd<-apply(xd,2,as.numeric)
	rownames(xd)<-geneNames
	colnames(xd)<-sampleNames

	if(!is.na(designFile)){
		x<-read.table(designFile,sep="\t",header=T,row.names=1,fill=T,,stringsAsFactors=FALSE)
		xd<-xd[,sort.list(colnames(xd))]
		xd<-xd[,colnames(xd) %in% rownames(x)]
		x<-x[rownames(x) %in% colnames(xd),]
		x<-x[sort.list(rownames(x)),]
		classes<-as.data.frame(x)
	}
	
	if(sum(apply(xd,2,is.na))>0 & impute){
		library(impute)
		allAnn<-dimnames(xd)
		data.imputed<-impute.knn(as.matrix(xd))$data
		xd<-data.imputed[1:features,]
		dimnames(xd)<-allAnn
	}
	
	return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}

standardize<-function(x){
	annAll<-dimnames(x)
	x<-scale(x)
	dimnames(x)<-annAll
	return(x)
}


dwd<-function(x,y){

#################
# requires a very specific directory structure in the root directory on C:
# file names and other run time requirements are necessitated by the exe
#
# the executable BatchAdjustSM.exe is available here: https://genome.unc.edu/pubsup/dwd/DWD.zip
# and more information is available here: https://genome.unc.edu/pubsup/dwd/
#
# the executable should be installed in C:/DWD/lib
# create the directory C:/DWD/DWDdata
#################

	y<-factor(y,levels=c("1","2"))
	nas<-sum(is.na(y))
	if(nas>0){ print("ERROR: Class labels must be '1' and '2'")
						 return(NA) }
	
	y<-as.numeric(as.vector(t(y)))

	x<-rbind(y,x)
	
	cwd<-getwd()
	setwd("C:/DWD/DWDdata")
	write.table(x,"DWD_Input.txt",sep="\t",col.names=F,row.names=F)

	twd<-getwd()
	write(twd,"fileDir.txt")

	y<-0
	write(y,"DWD_MeanAdjustType.txt")

	setwd("C:/DWD/lib")
	system("BatchAdjustSM.exe",invisible=T)

	setwd("C:/DWD/DWDdata")
	unlink("fileDir.txt")
	unlink("DWD_MeanAdjustType.txt")
	unlink("DWD_Input.txt")

	setwd("C:/DWD/DWDdata")
	dwdVec<-read.table("DWD_Vec.txt",sep="\t",header=F)

	setwd(cwd)

	return(t(dwdVec))
}

medianCtr<-function(x){
	annAll <- dimnames(x)
	medians <- apply(x,1,median,na.rm=T)
	x <- t(scale(t(x),center=medians,scale=F))
	dimnames(x) <- annAll
	return(x)
}

overlapSets<-function(x,y){

	# subset the two lists to have a commonly ordered gene list
	x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
	y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]

	#and sort such that thing are in the correct order
	x<-x[sort.list(row.names(x)),]
	y<-y[sort.list(row.names(y)),]

	return(list(x=x,y=y))
}


assignDiffScore.dwd<-function(x,y){
	both<-overlapSets(y$xd,x$xd) 	# get the overlap of genes
	both$x<- apply(both$x,2,function(x){sign(x)*sqrt(x^2/sum(x^2))}) 	# set the distance to the origin to 1
	both$y<- apply(both$y,2,function(x){sign(x)*sqrt(x^2/sum(x^2))})	# set the distance to the origin to 1
	msproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,1])	# project the samples on the MS-pL axis
	mlproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,2])	# project the samples on the pL-mL axis
	diffScore<- mlproj - msproj 
	return( diffScore )	# return the point on the differentiation axis
}

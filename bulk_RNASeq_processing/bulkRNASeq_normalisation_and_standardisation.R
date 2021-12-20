#QuantileNormalization Function from Joel Parker

# Step 1: Source the code
# Step 2: Save your data matrix ( no gene names or sample ids- order is preserved) as .csv file
# Step 3:x=as.matrix((read.csv('TCGA_Matrix.csv',header =F)))
# Step 4:x.norm <- quartileNorm(x)
# Step 5:write.csv(x=x.norm, file = "quartileNorm_TCGAMatrix.csv")

#This function does QUARTILE Normalization. Change .75 argument to any ything to normalize.
quartileNorm<- function(x,y=NA){
  uqs<- apply(x,2,function(x){ quantile(x[x>0 & !is.na(x)],0.75)})
  if(is.na(y)){
    y<- median(uqs)
  }
  x.norm <- t(apply(x,1,function(x,y){x*y},y/uqs))
  dimnames(x.norm)<- dimnames(x)
  return(x.norm)
}

## Then, with x as a matrix of counts, call it like this to normalize to the overall upper quartile:
# x.norm <- quartileNorm(x)
# Write out the file


#Column Standardize
#
x<-read.table("TCGA_filtered-normalized-log-mediancentered.txt", sep="\t", header = TRUE, row.names = 1)
standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}
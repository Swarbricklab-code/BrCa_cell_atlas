
pam50.symbol2symbol.v2 = function(x, s, cutoff=NULL, TS=NULL, ER=NULL, HER2=NULL, PR=NULL, GeneSignature=Pam50, gene.sigma = Pam50$subgroupQuantile, QT = F, export = F){
	# Pam50 (Parker et al 2009: Agilent human 1Av2 microarrays or custom-designed Agilent human 22k 2 channel; symbols in centroids) 
	# classify data indexed by gene symbols
	# Calculate ROR scores in Parker et al 2009 
    # ** Applied in BC cohort consist only centain population, e.g. Triple negative BCs (80% overlap basal)
	# x: expression matrix rownames are probes, colnames are arrayIDs
	# gene.sigma: dataframe, subpopulation quantile per gene computed from original Parker et al 2009 UNC traningn set (instead of centering to median, center to the gene-specic quantile)
    # new sample assigned to the subtype with highest centroids correlation (spearman),  "unclassifed" when corr < cutoff
	# TS: numeric vector. Tumor size for samples in x for ROR-C model. Encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
	# ER: ER status. 1: positive, 0:negative
	# HER2: HER2 status. 1: positive, 0:negative
	# PR: PR status. 1: positive, 0:negative
	# Author: Xi Zhao. July 16, 2012
        
    
	GeneSig = GeneSignature$GeneSignature
	orgGenes = rownames(GeneSig)
    
	# cat(">> Extracting matched PAM50 gene expression matrix ... ")
	overlapgenes = intersect(orgGenes, rownames(x))
	x.m = x[overlapgenes, , drop = F]
	GeneSig.m = GeneSig[overlapgenes, ]
    
	if(QT){
        nTotal = length(unique(rownames(GeneSig)))
        nMatch = length(unique(rownames(x.m)))
        cat("   mapping covarage: ", nMatch,"/",nTotal, "≈", round(nMatch/nTotal*100, 0),"%", sep="","\n")
    } 	
	
    gene.sigma.o = gene.sigma[rownames(x.m), s]
    x.sigma = unlist(llply(1:nrow(x.m), function(i) quantile(x.m[i,], probs = gene.sigma.o[i], na.rm = T)))
    x.m = sweep(x.m, 1, x.sigma)
  
    write.csv(x.m, "Geneexpressionmatrix.csv", sep="")
    
	# cat(">> Calibrating subtypes for single sample ... \n")
	out = getSSP(x.m, GeneSig.m, cutoff, method="spearman")
	
	# cat(">> ROR scoring by Relapse risk Prediction Models (Parker et al 2009) ... \n")
	out = data.frame(out, getPam50ROR(CentroidCorrMat=out[, 1:5], TS)[, -c(1:5)], stringsAsFactors = F)
	
	#cat(">> Producing quality-check statistics ...\n")
    if(QT) getSubtypeQC(x.m, out, ER, HER2, PR, subtypeCode=GeneSignature$subtypeCode, signature="pam50")
	#cat("[DONE]", "\n")
	
	if(export) write.csv(out, paste("Subtyping_output_subpopulationSpecific_", s, ".csv", sep=""), row.names = F)
	return(out)	
}


#-------------------------------------

getSSP = function(x.m, GeneSig.m, cutoff = NULL, method){
	# get "Single Sample Prediction" from Intrinsic or Pam50 signature
	# x.m: matched expression matrix, no duplicate samples
	# GeneSig.m: matched GeneSignature
	# revised on Oct 24 2012
	
	if(nrow(GeneSig.m)!=nrow(x.m)) stop("num of gene in GeneSig.m differed in num of genes in x.m")
	
	centroidCorr = vector()
	out = data.frame(matrix(NA, ncol=ncol(GeneSig.m), nrow = ncol(x.m)))
	colnames(out) = paste("Corr_",colnames(GeneSig.m),sep="")

	for(i in 1:ncol(GeneSig.m)){
		if(method=="pearson")   for(j in 1:ncol(x.m)) centroidCorr[j] = cor(x.m[,j], GeneSig.m[,i], method="pearson", use="na.or.complete")
		if(method=="spearman")  for(j in 1:ncol(x.m)) centroidCorr[j] = cor(x.m[,j], GeneSig.m[,i], method="spearman",use="pairwise.complete.obs")
		if(method=="euclidean") for(j in 1:ncol(x.m)) centroidCorr[j] = -c(dist(rbind(x.m[,j],GeneSig.m[,i]), method= "euclidean"))
		out[,i] = centroidCorr
	}
	out$call = out$maxCorr = apply(out, 1, max)

	if(!is.null(cutoff)){
		for(i in 1:nrow(out)){
            if(out$call[i]>=cutoff) out$call[i] = sample(colnames(GeneSig.m)[which(out[i,1:ncol(GeneSig.m)]== out$call[i])])[1] # random select 1 highest in case multiple hits
            if(out$call[i] <cutoff) out$call[i] = "unclassifed"
		}
	}else{
        for(i in 1:nrow(out)) out$call[i] = sample(colnames(GeneSig.m)[which(out[i,1:ncol(GeneSig.m)]== out$call[i])])[1] # random select 1 highest in case multiple hits
	}
	
    out = data.frame(sample = colnames(x.m), out, stringsAsFactors = F)

	return(out)
}



#-------------------------------------

getPam50ROR = function(CentroidCorrMat, TS=NULL){
	# ROR score from ROR model (Parker et al 2009): 1. for subtype only model: "ROR-S","ROR-S Group". 2. for combined model (Subtype + Clinical): "ROR-C","ROR-C Group"
	# CentroidCorrMat: correlation matrix for each centroid correlation
	# TS: numeric vector. Tumor size for samples in x for ROR-C model. Encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
	# thresholds from the training set that required no LumA sample to be in the high-risk group and no basal-like sample to be in the low-risk group.
	
	col.basal = grep("Basal", colnames(CentroidCorrMat), ignore.case = T)
	col.her2  = grep("Her2",  colnames(CentroidCorrMat), ignore.case = T)
	col.lumA  = grep("LumA",  colnames(CentroidCorrMat), ignore.case = T)
	col.lumB  = grep("LumB",  colnames(CentroidCorrMat), ignore.case = T)

	glthreshold = -0.15 
	ghthreshold =  0.1
	genomic = 0.04210193*CentroidCorrMat[, col.basal] + 0.12466938*CentroidCorrMat[, col.her2] + -0.35235561*CentroidCorrMat[, col.lumA] + 0.14213283*CentroidCorrMat[, col.lumB]	
	griskgroups = genomic
	griskgroups[genomic>ghthreshold] = "high"
	griskgroups[genomic>glthreshold & genomic<ghthreshold] = "intermediate"
	griskgroups[genomic<glthreshold] = "low"
	genomic = 100*(genomic+0.35)/0.85		
	griskgroups_code = griskgroups; 
	
    labelval = c(1, 0.5, 0); names(labelval) = c("high", "intermediate", "low")
    labelval = sort(labelval[unique(griskgroups)], decreasing=T)
	
    CentroidCorrMat$ROR_S = genomic; CentroidCorrMat$ROR_S_group = griskgroups; 
	CentroidCorrMat$ROR_S_group_code = as.numeric(as.character(factor(griskgroups, labels = labelval)))

	if(!is.null(TS)){
		clthreshold = -0.1
		chthreshold =  0.2
		combined = 0.0442770*CentroidCorrMat[,col.basal] + 0.1170297*CentroidCorrMat[,col.her2] + -0.2608388*CentroidCorrMat[,col.lumA] + 0.1055908*CentroidCorrMat[,col.lumB] + 0.1813751*TS
		criskgroups<-combined
		criskgroups[combined>chthreshold]<-"high"
		criskgroups[combined>clthreshold & combined<chthreshold]<-"intermediate"
		criskgroups[combined<clthreshold]<-"low"		
		combined<- 100* (combined + 0.35 ) / 0.85
		
		CentroidCorrMat$ROR_C = combined; CentroidCorrMat$ROR_C_group = criskgroups
		CentroidCorrMat$ROR_C_group_code = as.numeric(as.character(factor(criskgroups, labels = labelval)))
	}
	return(CentroidCorrMat)
}

#-------------------------------------

getSigma.generic = function(x, d){ 
    # estimate subgroup specific quantile
    suppressMessages(require(plyr))
    quantile.subpopulation.j = function(s, j){
        sigma.i = function(i){
            a = data.frame(x[i,]); 
            colnames(a) = "expr"
            a$subpopulation = d[,j]
            #  a = na.omit(a)
            a.all.med = quantile(a$expr, probs=0.5, na.rm = T)
            x.s = subset(a, subpopulation == s)$expr
            ecdf(x.s)(a.all.med)
        }
        ldply(1:nrow(x), sigma.i)
    }
    
    out = llply(1:ncol(d), function(j){
        subpopulations = levels(factor(d[,j]))
        res = llply(subpopulations, function(s) quantile.subpopulation.j(s=s, j = j))
        res = do.call(cbind, res)
        colnames(res) = subpopulations; rownames(res) = rownames(x)
        res    
    })
    out = do.call(cbind, out)

    out
}

#-------------------------------------

signature.pam50 =  function(signatureFile="Centroids_pam50.txt", gene.sigma = NULL){
	# PAM50 gene signature
	GeneSig = read.delim(signatureFile, row.names=1, as.is=T)
	subtypeCode =  data.frame(level=c("Basal","Her2", "LumA", "LumB","Normal"), 
							  color=c("red", "purple", "darkblue","cyan", "green"), stringsAsFactors=F)
    # subpopulation-specific subyping for PAM50
    if(is.null(gene.sigma)){
        gene.sigma = read.csv("Pam50_subpopulation_specific_gene_quantile.csv", row.names=1, as.is=T, check.names = F) 
    }
    
	pam50 = list(GeneSignature=GeneSig, 
        platform="AgilentExpr44k_2channel", 
        subtypeCode=subtypeCode, 
        subgroupQuantile = gene.sigma)
	return(pam50)
}


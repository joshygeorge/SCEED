library(edgeR)
library(scran)
library(DESeq)
library(Biobase)
library(SCnorm)

normalize_data <- function(eset,method = "tmm")
{
	cnts <- exprs(eset)
	if(method == "tmm")
	{
		s <-edgeR::calcNormFactors(cnts)*colSums(cnts)
  		n.cells<-ncol(cnts)
 		s <- n.cells*(s/sum(s))
 		nor.cnts <-t(t(cnts)/s)
	}
	else if(method == "deseq")
	{
		s <- DESeq::estimateSizeFactorsForMatrix(cnts)
  		n.cells<-ncol(cnts)
		s <- n.cells*(s/sum(s))
 		nor.cnts <-t(t(cnts)/s)
	}
	else if (method =="rpm")
	{
		s <-colSums(cnts)
  		n.cells<-ncol(cnts)
		s <- n.cells*(s/sum(s))
  		nor.cnts <-t(t(cnts)/s)
	}
	else if(method == "scran")
	{
      		clusters <- quickCluster(cnts,min.size = ncol(cnts)/5)
      		factors <- computeSumFactors(counts, cluster=clusters)
		nor.cnts = t(t(cnts)/factors*mean(factors))
    	}
	else if(method == "scnorm")
	{
		nor.cnts = SCnorm(Data = cnts)
    	}
	nor.eset <- ExpressionSet(nor.cnts,featureData = featureData(eset),phenoData = pData(eset))
	return(nor.eset)
}

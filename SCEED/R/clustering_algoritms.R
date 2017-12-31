library(Seurat)
library(dplyr)
library(Matrix)
library(SIMLR)

sceed_kmeans <- function(eset,num_clusters = 2)
{

	cnt <- exprs(eset)
	cnt.log2 <- log2(cnt +1)
	r.sd <- apply(cnt.log2,1,sd)
	o <- order(r.sd,decreasing=T)
	exp.sel <- t(cnt.log2[o[1:500],])
	clus.res <- kmeans(exp.sel,centers=num_clusters)
	return(clus.res$cluster)
}


sceed_seurat <- function(eset,num_clusters=2)
{
	library(Seurat)
	library(dplyr)
	library(Matrix)
	library(SIMLR)
	pbmc <- CreateSeuratObject(raw.data = exprs(eset), min.cells = 3, min.genes = 200, project = "Seurat_SCEED")
	pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000,display.progress=FALSE)
	pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot=FALSE,display.progress=FALSE)
	pbmc <- ScaleData(object = pbmc)
	pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
	pbmc <- FindClusters(object = pbmc, k.param = num_clusters, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
	clus.res <- pbmc@ident
	return(clus.res)
}


sceed_simlr <- function(eset,num_clusters=2)
{
	library(Seurat)
	library(dplyr)
	library(Matrix)
	library(SIMLR)
	cnts <- exprs(eset)
	clus = SIMLR_Large_Scale(cnts, c = num_clusters, k = 10, kk = 100, if.impute = FALSE, normalize = FALSE)
	clus.res = clus$y$cluster
	return(clus.res)
}



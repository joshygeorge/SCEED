library(splatter)
library(checkmate)
library(Biobase)

generateDataSet <- function(nGenes =10000, nCells =100, group.prop =c(0.5,0.5),num.marker.genes.per.group = c(50,50), foldChange = c(2,2))
{
	if(length(group.prop) != length(num.marker.genes.per.group) )
	{
		stop("number of marker genes per group should match the number of groups")
	}
	if(length(group.prop) != length(foldChange) )
	{
		stop("foldChange of marker genes should match the number of groups")
	}
	if(sum(group.prop) != 1 )
	{
		stop("sum of the proportions of cells should be equal to one")
	}
	num.marker.genes <- sum(num.marker.genes.per.group)
	if(num.marker.genes > nGenes )
	{
		stop("total number of marker genes should be less than the total genes")
	}
	splat = newSplatParams(nGenes = 10000, batchCells = nCells, dropout.present = T,seed=sample(1:2000,1))
	sim.data = splatSimulate(splat, verbose = F)
	num.cells <- integer(length = length(group.prop))
	for(i in 1:length(group.prop))
	{
		num.cells[i] <- round(nCells*group.prop[i])
	}
	num.cells[1] <- nCells - sum(num.cells[2:length(group.prop)])
	cnts = counts(sim.data)
	marker.status = character(length = nGenes)
	group.idx <- list()
	genes.all <- 1:nGenes
	for(j in 1:length(num.marker.genes.per.group))
	{
		group.idx[[j]] <- sample(genes.all, num.marker.genes.per.group[j])
		genes.all <- setdiff(genes.all,group.idx[[j]])
		marker.status[group.idx[[j]]] <- paste("Group",j,sep="")
	}
	
	probe.anns <- data.frame(primerid = rownames(cnts), markerof = marker.status)
	rownames(probe.anns) <- rownames(cnts)
	cell.group <- NULL 
	for(j in 1:length(group.prop))
	{
		cell.group <- c(cell.group,rep(paste("Group",j,sep=""),num.cells[j]))
	}
	cell.anns <- data.frame(wellKey =colnames(cnts), group = cell.group)
	rownames(cell.anns) <- colnames(cnts)
	pheno.data = new("AnnotatedDataFrame",cell.anns)
	probe.data = new("AnnotatedDataFrame",probe.anns)

	# Add marker gene expression 
	groups = unique(cell.group)
	logCnts <- log2(cnts +1)
	for(j in 1:length(groups))
	{
		logFC <- log2(foldChange[j])
		genes.idx <- which(marker.status == groups[j])
		cells.idx <- which(cell.group == groups[j])
		num <- length(genes.idx)*length(cells.idx)
		x <- rnorm(num,mean=logFC,sd=0.1)
		logCnts[genes.idx,cells.idx] <- logCnts[genes.idx,cells.idx] + x	
	}
	cnts <- 2**logCnts
	eset <- ExpressionSet(cnts,featureData = probe.data,phenoData = pheno.data)
	return(eset)
}




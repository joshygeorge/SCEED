library(MAST)
library(data.table)

find.de.genes <- function(eset)
{
	mastOBJ = FromMatrix(exprs(eset),cData = pData(eset), fData = fData(eset) )
	zlmCond <- zlm(~group,mastOBJ)
	summaryCond <- summary(zlmCond, doLRT='groupGroup2')
	summaryDt <- as.data.frame(as.list(summaryCond)$datatable)
	contrast <- as.character(summaryDt$contrast)
	component <- as.character(summaryDt$component)
        pval <- summaryDt[contrast=='groupGroup2' & component=='H',]
        logFC <- summaryDt[contrast=='groupGroup2' & component=='logFC',]
	all(logFC$primerid == pval$primerid)
	{
        	fcHurdle <- cbind(pval,logFC)
	}

	FDR <- p.adjust(fcHurdle[,4],"fdr")
	fcHurdle <- data.frame(fcHurdle[,c(1,4,15)],FDR)
	rownames(fcHurdle) <- as.character(fcHurdle[,1])
	probe.anns <- as(featureData(eset),"data.frame")
	ids <- intersect(rownames(fcHurdle),rownames(probe.anns))
	fcHurdle.ord <- fcHurdle[ids,]
	probe.anns.ord <- probe.anns[ids,]
	probe.stats <- data.frame(probe.anns.ord, fcHurdle.ord)
	return(probe.stats)
}


de.stats <- function(eset)
{
	mastOBJ = FromMatrix(exprs(eset),cData = pData(eset), fData = fData(eset) )
	zlmCond <- zlm(~group,mastOBJ)
	summaryCond <- summary(zlmCond, doLRT='groupGroup2')
	dt <- summaryCond[[1]]
	return(dt)
}



outputAttractors <- function(x, min.basins=2, strength.pos=10, outputGeneNumber=100, write2File=FALSE, fileName="attractors.txt"){
	m <- nrow(x)
	outMat <- NULL
	k <- 1
	for(i in 1:m){
		o <- order(x[i,], decreasing=TRUE)
		if(rownames(x)[i] %in% colnames(x)[o[1:min.basins]]) next
		outMat <- cbind(outMat,colnames(x)[o[1:outputGeneNumber]], round(x[i,o[1:outputGeneNumber]], 4))
		aname <- paste("Attractor", sprintf("%03d",k), sep="")
		colnames(outMat)[c(2*k-1, 2*k)] <- c(paste(aname, ":Genes",sep="" ), paste(aname, ":MI",sep="") )
		k <- k+1
	}
	rownames(outMat) <- 1:outputGeneNumber
	if(write2File){
		write.table(outMat, file=fileName, sep='\t', quote=F)
	}
	return (outMat)
}

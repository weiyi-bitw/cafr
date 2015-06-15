FindGLAttractor <- function(data,seed, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE){
  if(! seed %in% rownames(genome)) stop("Cannot find seed gene in genome file rownames!")
  genes.genome <- rownames(genome)
  idx.seed <- which(rownames(genome)==seed)
  data <- data[intersect(genes.genome, rownames(data)),] # sort the rownames according to genomic location
  n <- ncol(data)
  m <- nrow(data)
  mg <- length(genes.genome)
  
  idxrange <- max(1, idx.seed-windowSize):min(mg, idx.seed + windowSize) 
  generange <- genes.genome[idxrange]
  data <- data[intersect(generange, rownames(data)),]
  
  out <- attractorScanningGL(data, genome, alpha, windowSize, maxIter, epsilon, bin, so, score.position, num.output, negateMI, verbose, saveAllScore=TRUE)
  sumOut <- out$summary
  if(nrow(sumOut) > 1){
    idxtop <- sapply(sumOut[,2], function(g){which(genes.genome==g)})
    select <- which.min(abs(idxtop - idx.seed))
  }else{
    select <- 1
  }
  genes <- sumOut[select,1:num.output]
  mis <- out$scoremat[,rownames(sumOut)[select]]
  
  names(mis) <- genes
  
  return(mis)
  
}

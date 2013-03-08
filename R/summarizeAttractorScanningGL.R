summarizeAttractorScanningGL <- function(out, genes.genome, windowSize=50){
  topIdx <- sapply(out$attractome[1,], function(x){which(genes.genome==x)})
  o <- order(out$score, decreasing=T)
  sumOut <- NULL
  center <- NULL
  while(length(o) > 0){
    o.top <- o[1]
    #print(length(o))
    #print(o.top)
    sumOut <- rbind(sumOut,c(out$attractome[,o[1]], out$score[o[1]]))
    center <- c(center, colnames(out$attractome)[o[1]])
    idxDiff <- abs(topIdx - topIdx[o.top])
    killIdx <- which(idxDiff <= (windowSize))
    if(length(killIdx) == 0) break
    o <- setdiff(o, killIdx)
  }
  rownames(sumOut) <- center
  return (sumOut)
}

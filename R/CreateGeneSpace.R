CreateGeneSpace <- function(ge, oncogenes, map){
  ng <- length(oncogenes)
  gSpace <- matrix(NA, nrow=ng, ncol=ncol(ge))
  dimnames(gSpace) <- list(oncogenes, colnames(ge))
  mappedGenes <- map[rownames(ge), "Gene.Symbol"]
  for(i in 1:ng){
    g <- oncogenes[i]
    idx <- which(mappedGenes == g)
    if(length(idx)==1){gSpace[i,] <- ge[idx,]}
    else{gSpace[i,] <- apply(ge[idx,], 2, mean)}
  }
  return (gSpace)
}

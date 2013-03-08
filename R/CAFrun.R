CAFrun <- function(data, vec, a=5, maxIter = 100, epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  negateMI = TRUE, verbose=TRUE, sorting=TRUE){
  m <- nrow(data)
  n <- ncol(data)
  
  if(rankBased){
    vec <- rank(vec)
    dataIn <- t( apply(data, 1, rank) )
  }
  miOut <- .Call("cafR2C", data, vec, a, maxIter, epsilon, m, n, bin, so, as.integer(negateMI), as.integer(verbose))
  
  if(miOut[1] == -999) return (NULL)
  
  names(miOut) <- rownames(data)
  if(sorting){
    return (sort(miOut, decreasing=T))
  }else{
    return (miOut)
  }
}

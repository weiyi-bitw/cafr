FindCorrAttractor <- function(data, vec, a=10, maxIter=100, epsilon=1E-14, rankBased=FALSE){
  dataIn <- data
  if(rankBased){
    vec <- rank(vec)
    for(i in 1:nrow(data)){
      dataIn[i,] <- rank(data[i,])
    }
  }
  rs <- getAllCorWz(dataIn, vec, rankBased=rankBased)
  prers <- rs
  w <- abs(rs)^a / sum(abs(rs)^a)
  w[rs < 0] <- -w[rs < 0]
  metagene <- w %*% data
  if(rankBased) metagene <- rank(metagene)
  c <- 0
  while(c < maxIter){
    rs <- getAllCorWz(dataIn, metagene)
    delta <- sum((rs - prers)^2)
    cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
    print(rs[order(rs, decreasing=T)[1:20]])
    flush.console()
    if(delta < epsilon){
      break
    }
    prers <- rs
    w <- abs(rs)^a / sum(abs(rs)^a)
    w[rs < 0] <- -w[rs < 0]
    metagene <- w %*% data
    if(rankBased) metagene <- rank(metagene)
    c <- c + 1
  }
  if(c >= maxIter) return (NULL)
  return (sort(rs, decreasing=T))
}

getAllCCDIWz <- function(x, observations, sorted=FALSE){
  n <- ncol(x)
  m <- nrow(x)
  out <- .C("concordance_index_all", predictions=as.double(x), observations=as.double(observations), mIn = as.integer(m), nIn = as.integer(n), score = rep(0, m))
  score <- out$score
  names(score) <- rownames(x)
  if(sorted) score <- sort(score, decreasing=T)
  return (score)
}

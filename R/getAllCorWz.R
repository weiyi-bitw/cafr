getAllCorWz <- function(x, vec, rankBased=FALSE, sorting=FALSE){
  m <- nrow(x)
  n <- ncol(x)
  if(length(vec) != n){stop("legnth of two vectors are different!")}
  
  if(rankBased){
    for(i in 1:m){
      x[i,] <- rank(x[i,])
    } 
    vec <- rank(vec)
  }
  garbage <- rep(-999, m)
  out <- .C("getAllCorWz", data = as.double(x), vec = as.double(vec), m = as.integer(m), n = as.integer(n), rs=as.double(garbage))
  
  rs <- out$rs
  names(rs) <- rownames(x)
  if(sorting) rs <- sort(rs, decreasing=TRUE)
  return (rs)
}

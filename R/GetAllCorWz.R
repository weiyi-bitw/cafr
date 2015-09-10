GetAllCorWz <- function(data, x, rankBased=FALSE, sorting=FALSE){
  m <- nrow(data)
  n <- ncol(data)
  if(length(x) != n){stop("legnth of two vectors are different!")}
  
  if(rankBased){
    for(i in 1:m){
      data[i,] <- rank(data[i,])
    } 
    x <- rank(x)
  }
  garbage <- rep(-999, m)
  out <- .C("GetAllCorWzC", data = as.double(data), vec = as.double(x), m = as.integer(m), n = as.integer(n), rs=as.double(garbage))
  
  rs <- out$rs
  names(rs) <- rownames(data)
  if(sorting) rs <- sort(rs, decreasing=TRUE)
  return (rs)
}

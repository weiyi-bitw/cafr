Corr <- function(x, y, rankBased = FALSE){
  n <- length(x)
  if(length(y) != n){stop("length of two vectors are different!")}
  if(rankBased){
    x <- rank(x)
    y <- rank(y)
  }
  
  out <- .C("corR", x=as.double(x), y=as.double(y), n=as.integer(n), cOut=as.double(0))
  
  r <- out$cOut
  return(r)
}

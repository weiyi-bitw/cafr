getAllMIWz <- function(data, x, bin=6, so=3, rankBased=FALSE, normalize = TRUE, sorting = FALSE, negateMI = FALSE){
  m <- nrow(data)
  n <- ncol(data)
  if(length(x) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  if(rankBased){
    for(i in 1:m){
      data[i,] <- rank(data[i,])
    } 
    x <- rank(x)
  }
  garbage <- rep(-1, m)
  out <- .C("getAllMIWz_R", data = as.double(data), vec = as.double(x), mi = as.double(garbage), m = as.integer(m), n = as.integer(n), bin = as.integer(bin),so = as.integer(so), norm=as.integer(normalize), negateMI = as.integer(negateMI))
  
  mi <- out$mi
  names(mi) <- rownames(data)
  if(sorting) mi <- sort(mi, decreasing=TRUE)
  return (mi)
}

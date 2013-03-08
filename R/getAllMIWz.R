getAllMIWz <- function(x, vec, bin=6, so=3, rankBased=FALSE, normalize = TRUE, sorting = FALSE, negateMI = FALSE){
  m <- nrow(x)
  n <- ncol(x)
  if(length(vec) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  if(rankBased){
    for(i in 1:m){
      x[i,] <- rank(x[i,])
    } 
    vec <- rank(vec)
  }
  garbage <- rep(-1, m)
  out <- .C("getAllMIWz_R", data = as.double(x), vec = as.double(vec), mi = as.double(garbage), m = as.integer(m), n = as.integer(n), bin = as.integer(bin),so = as.integer(so), norm=as.integer(normalize), negateMI = as.integer(negateMI))
  
  mi <- out$mi
  names(mi) <- rownames(x)
  if(sorting) mi <- sort(mi, decreasing=TRUE)
  return (mi)
}

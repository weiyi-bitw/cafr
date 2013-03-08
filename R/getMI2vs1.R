getMI2vs1 <- function(x, y, z, bin=6, so=3, normalize=FALSE){
  n <- length(x)
  if(length(y) != n | length(z) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  out <- .C("mi2vs1", x = as.double(x), y=as.double(y), z=as.double(z), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), mi = 0, norm=as.integer(normalize))
  
  mi <- out$mi
  return (mi)
  
}

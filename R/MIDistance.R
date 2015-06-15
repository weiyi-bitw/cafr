MIDistance <- function(x, y, bin=6, so=3){
  n <- length(x)
  if(length(y) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  out <- .C("mid2R", x = as.double(x), y=as.double(y), n=as.integer(n), 
            bin=as.integer(bin), so=as.integer(so), mi_out=0)
    
  mid <- out$mi_out
  return (mid)
  
}

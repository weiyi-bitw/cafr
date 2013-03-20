getMI <- function(x, y, bin=6, so=3, rankBased=FALSE, normalize=TRUE, negateMI = FALSE){
  n <- length(x)
  if(length(y) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  if(is.factor(y)){
    yph <- y
    idx <- which(!is.na(yph))
    x <- x[idx]
    yph <- yph[idx]
    
    y <- as.numeric(yph)
    l <- levels(yph)
    nl <- length(l)
    if(is.ordered(yph)){soph <- 3}
    else {soph <- 1}
    out <- .C("mi2DiffBins", x = as.double(x), y = as.double(y), n = as.integer(length(yph)), binx = as.integer(bin), biny = as.integer(nl), sox = as.integer(so), soy = as.integer(soph), miOut = 0, norm = as.integer(normalize), negateMI = as.integer(negateMI))
  }else{
    
    if(rankBased){
      x <- rank(x)
      y <- rank(y)
    }
    out <- .C("mi2R", x = as.double(x), y=as.double(y), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), miOut = 0, norm = as.integer(normalize), negateMI = as.integer(negateMI))
    
  }
  
  mi <- out$miOut
  return (mi)
  
}

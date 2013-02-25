getCCDIdx <- function(predictions, observations){
  n <- length(predictions)
  out <- .C("concordance_index", predictions=as.double(predictions), observations=as.double(observations), nIn = as.integer(n), c = as.double(0))
  c <- out$c
  return (c); 
}

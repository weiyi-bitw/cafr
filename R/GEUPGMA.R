GEUPGMA <- function(sim, ge) {
  trc <- .Call("GEUPGMACC", sim, ge)
  cut <- which(abs(trc[,1] + 999) < 1E-10)[1]
  trc <- trc[1:(cut-1), ]
  return (trc)
}

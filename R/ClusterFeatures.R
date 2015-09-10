ClusterFeatures <- function(ge, trc, cutoff=-1) {
  if (cutoff == -1) {
    grad <- diff(trc[,3])
    end.idx <- whic.min(grad)
  } else {
    end.idx <- min(which(trc[,3] < cutoff))-1
  }
  cutoff.out <- trc[end.idx, 3]

  clusters <- as.list(rownames(ge))
  names(clusters) <- 0:(nrow(ge)-1)
  for (i in 1:end.idx) {
    qo <- as.character(trc[i,])
    clusters[[qo[1]]] <- c(clusters[[qo[1]]], 
                               clusters[[qo[2]]])
    clusters[[qo[2]]] <- NULL
  }
  return (list(clusters=clusters, cutoff=cutoff.out))
}

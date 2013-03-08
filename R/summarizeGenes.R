summarizeGenes <- function(ge, map, sumfun=median){
  genes <- sort(unique(map))
  #genes = genes[-which(genes==""),1]
  geids <- rownames(ge)
  mapids <- rownames(map)
  m <- length(genes)
  n <- dim(ge)[2]
  out <- matrix(0, nrow=m, ncol=n)	
  cnt <- 0
  glist <- list()	
  for(g in genes){
    ids <- mapids[which(map==g)]
    idx <- which(geids %in% ids)
    if (length(idx) > 0){
      cnt <- cnt + 1
      glist[cnt] <- g
      # cat(idx,"\n")
      # flush.console()
      if(length(idx)==1){
        out[cnt,] <- ge[idx,]
      }else{
        out[cnt,] <- apply(ge[idx,],2,sumfun)
      }
    }
  }
  out <- out[1:cnt,]
  rownames(out) <- unlist(glist)
  colnames(out) <- colnames(ge)
  return (out)
}

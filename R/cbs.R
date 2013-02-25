## SEEMS INCOMPLETE?

cbs <- function(ge, map, sumfun=median, corTh = 0.9){
  genes <- sort(unique(map))
  geids <- rownames(ge)
  mapids <- rownames(map)
  m <- length(genes)
  n <- dim(ge)[2]
  out <- NULL
  cnt <- 0
  glist <- list()
  for(g in genes){
    cnt <- cnt + 1
    ps <- rownames(map)[which(map==g)]
    nps <- length(ps)
    if(nps==1){
      out <- rbind(out, ge[ps,])
      glist[cnt] <- g
      plist[cnt] <- ps
    }else{
      #cc = 
      
    }
  }
  
}

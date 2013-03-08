findCytoband <- function(sumOut, genome){
  genes.genome <- rownames(genome)
  cytoband <- apply(sumOut, 1, 
                    function(a, genes.genome, genome){
                      genes <- a[1:(length(a)-1)]
                      idx <- which(genes.genome %in% genes)
                      cys <- genome[idx, "Cytoband"]
                      cys <- cys[!is.na(cys)]
                      start <- cys[1]
                      end <- cys[length(cys)]
                      if(start == end){
                        cy <- paste(genome[idx[1],"Chr"], start, sep="")
                      }else{
                        cy <- paste(genome[idx[1],"Chr"], start, "-", end, sep="")
                      }
                      return(cy)
                    }, genes.genome=genes.genome, genome=genome
  )
  out <- cbind(cytoband, sumOut)
  return (out)
}

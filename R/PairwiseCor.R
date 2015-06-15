PairwiseCor <- function(data, idx.start, buffer.exp=16){
  m <- nrow(data)
  n <- ncol(data)
  all.tasks <- m*(m-1)/2
  buffer.size <- bitwShiftL(1, buffer.exp)
  garbage <- matrix(-999, nrow=3, ncol=buffer.size, 
                    dimnames=list(c("gene1", "gene2", "cor")))
  out <- .C("PairwiseCor", 
            data=as.double(data),
            idx_start=as.integer(idx.start),
            all_tasks=as.integer(all.tasks),
            m=as.integer(m),
            n=as.integer(n),
            buffer_exp=as.integer(buffer.exp),
            out=as.double(garbage))
  cor.table <- t(matrix(out$out, nrow=3, ncol=buffer.size))
  cor.table[,1] <- cor.table[,1]+1
  cor.table[,2] <- cor.table[,2]+1
  return (list(cor.table=cor.table, idx.end=out$idx_start))
}

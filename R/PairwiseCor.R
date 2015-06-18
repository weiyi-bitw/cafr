PairwiseCor <- function(data, idx.start=0, buffer.exp=16){
  m <- nrow(data)
  n <- ncol(data)
  all.tasks <- m*(m-1)/2
  buffer.size <- min(bitwShiftL(1, buffer.exp), all.tasks-idx.start)
  out <- .Call("PairwiseCor", 
            data, idx.start, all.tasks, m, n, buffer.exp)
  cor.table <- t(matrix(out, nrow=3, ncol=buffer.size))
  cor.table[,1] <- cor.table[,1]+1
  cor.table[,2] <- cor.table[,2]+1
  return (list(cor.table=cor.table, idx.end=idx.start + buffer.size))
}

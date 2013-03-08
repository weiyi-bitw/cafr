loadGs <- function(file, rownames, colnames, sep="\t", header=T){
  gs <- as.matrix(read.table(file, sep=sep, header=header))
  n <- dim(gs)[1]
  B <- matrix(0, nrow=length(rownames), ncol=length(colnames))
  rownames(B) <- rownames
  colnames(B) <- colnames
  for(i in 1:n){
    tf <- gs[i,1]
    tg <- gs[i,2]
    if((tf %in% colnames) & (tg %in% rownames)){
      B[tg, tf] <- 1
    }
  }
  return (B)
}

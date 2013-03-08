loadTextMat <- function(file, sep="\t"){
  mset <- read.table(file=file, sep=sep, as.is=TRUE)
  colNames <- mset[1, -1]
  rowNames <- mset[-1, 1]
  mset <- mset[-1,-1]
  mset <- as.matrix(mset)
  colnames(mset) <- colNames
  rownames(mset) <- rowNames
  mset
}

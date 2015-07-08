AllPairwiseCor <- function(ge) {
  # calculate all pairwise Pearson correlation of a gene expression matrix
  # result stored in SQLite database table
  # 
  # Args:
  #   ge: gene expression matrix with genes in rows and samples in columns
  #
  # Return:
  #   A vector of all pair wise correlation, with the elements mapped to
  #   0: Cor(0:1), 1: Cor(0:2) ... etc
  #
  .Call("AllPairwiseCorCC", ge) 
}

varFilter <- function(mset,dim=1, numFilterIn=10000){
  if(dim==1){
    mset <- t(mset)
  }
  msetSd <- sd(mset)
  idx <- sort(msetSd, decreasing=TRUE, index.return=TRUE)$ix[1:numFilterIn]
  mset <- mset[,sort(idx)]
  if(dim==1){
    mset <- t(mset)
  }
  mset
}

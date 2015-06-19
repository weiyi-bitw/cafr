HierarchicalClustering <- function(ge, sim.db, cut.off=NULL) {
  # Performing hierarchical clustering using average linkage based on similariry 
  # database
  # 
  # Args:
  #   ge: gene expression matrix with genes in rows and samples in columns
  #   sim.db: name of similarity database with the following table
  #     allPairwiseSim and indexMap, as described in PairwiseCor.R
  #
  # Return:
  #   A list of:
  #   matrix hclust.trace, containing:
  #     1. idx1, idx2: indices to be merged
  #     2. sim: similarity measurement between idx1 and idx2
  #     3. diff: similarity drop compared with previous merge
  #   a list of clustered features based on the cut.off similarity, if no 
  #     cut.off is given, using the greatest drop in similarity
  #
  #

  
  cat( "Copying similarity database ... \n" )
  flush.console()
  db.tmp <- paste(sim.db, "tmp", sep=".")
  file.copy(sim.db, db.tmp)

  conn <- dbConnect(SQLite(), dbname=db.tmp)
  
  m <- nrow(ge)
  index.map <- dbReadTable(conn=conn, "indexMap", row.names=1)
  hclust.trace <- matrix(NA, nrow=m, ncol=4,
                         dimnames=list(1:m, c("idx1", "idx2", "sim", "diff")))
  ge.copy <- ge
  
  cat("Clear low similarity pairs")
  dbGetQuery(conn, "DELETE FROM allPairwiseSim WHERE sim < 0")

  que <- "DELETE FROM allPairwiseSim WHERE idx1=? OR idx2=?"
  cc <- 1
  while ( TRUE ) {
    qo <- dbGetQuery(conn, 
                     "SELECT * FROM allPairwiseSim ORDER BY sim DESC LIMIT 1")
    if ( nrow(qo) == 0 ) {
      break
    } else if ( qo[1,3] <= 0 ) {
      break
    }
    hclust.trace[cc,c("idx1", "idx2", "sim")] <- as.numeric(qo[1,])
    
    idx1 <- qo[1,1]
    idx2 <- qo[1,2]
    to.delete <- data.frame(idx1=c(idx1, idx2), idx2=c(idx2, idx1))

    cat(sprintf("(%d) Merge %d into %d\n", cc, idx2, idx1))
    flush.console()

    ge.copy[idx1,] <- ge.copy[idx1,] + ge.copy[idx2,]
    ge.copy[idx2,] <- 0
    cor.out <- cbind(c(1:idx1, rep(idx1, m-idx1)), 
                     c(rep(idx1, idx1), (idx1+1):m), 
                     GetAllCorWz(ge.copy, ge.copy[idx1,]))
    good.idx <- !is.nan(cor.out[,3]) & cor.out[,3]>0 & (cor.out[,1] != cor.out[,2])
    if (sum(good.idx) < 1) {
      break
    }
    cor.out <- as.data.frame(matrix(cor.out[good.idx,], nrow=sum(good.idx)))

    dbBegin(conn)
    dbGetPreparedQuery(conn, que, bind.data=to.delete)
    dbCommit(conn)
    
    dbWriteTable(conn, "allPairwiseSim", cor.out, append=TRUE, 
                 overwrite=FALSE, row.names=FALSE)
    
    cc <- cc + 1

  }
  dbDisconnect(conn)
  hclust.trace <- hclust.trace[1:cc,]
  hclust.trace[2:cc, "diff"] <- diff(hclust.trace[,"sim"])
  file.remove(db.tmp)
  if ( !is.null(cut.off) ) {
    cutoff.index <- min(which(as.numeric(hclust.trace[,"sim"]) < cut.off))
  } else {
    cutoff.index <- which.min(as.numeric(hclust.trace[,"diff"]))
  }

  clusters <- as.list(index.map[,"name"])
  names(clusters) <- rownames(index.map)
  
  for ( i in 1:(cutoff.index-1) ) {
    qo <- as.character(hclust.trace[i,])
    clusters[[qo[1]]] <- c(clusters[[qo[1]]], clusters[[qo[2]]])
    clusters[[qo[2]]] <- NULL
  }

  return (list(clusters=clusters, hclust.trace=hclust.trace))
}

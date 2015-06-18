AllPairwiseCor <- function(ge, db="cafr", buffer.exp=16) {
  # calculate all pairwise Pearson correlation of a gene expression matrix
  # result stored in SQLite database table
  # 
  # Args:
  #   ge: gene expression matrix with genes in rows and samples in columns
  #   db: name of the SQLite database
  #
  # Return:
  #   An SQLite database named after argument db, containing following tables
  #   allPairwiseSim: table of all pairwise correlation in form of 
  #     (idx1, idx2, sim), where sim is the similarity measure
  #   indexMap: mapping between index and row names
  #

  conn <- dbConnect(SQLite(), dbname=db)
  
  dbSendQuery(conn=conn, 
              "CREATE TABLE allPairwiseSim 
               (idx1 INTEGER, idx2 INTEGER, sim DOUBLE)" )

  dbSendQuery(conn=conn, 
              "CREATE TABLE indexMap (idx INTEGER, name TEXT)")
  index.map <- as.data.frame(row.names=1:nrow(ge), rownames(ge))
  dbWriteTable(conn=conn, name="indexMap", value=index.map, 
               row.names=TRUE, overwrite=FALSE, append=TRUE)

  idx.start <- 0
  m <- nrow(ge)
  tasks <- m*(m-1)/2
  while (idx.start < tasks) {
    cat(sprintf("Start index: %d / %d\n", idx.start, tasks))
    flush.console()
    out <- PairwiseCor(ge, idx.start, buffer.exp)
    idx.start <- out$idx.end
    dbWriteTable(conn=conn, name="allPairwiseSim", 
                 value=as.data.frame(out$cor.table),
                 row.names=FALSE, overwrite=FALSE, append=TRUE)
    gc()
  }
  dbSendQuery(conn=conn, 
              "CREATE INDEX idx_idx1 ON allPairwiseSim (idx1, sim)")
  dbSendQuery(conn=conn, 
              "CREATE INDEX idx_idx2 ON allPairwiseSim (idx2, sim)")
  dbSendQuery(conn=conn, "CREATE INDEX idx_sim ON allPairwiseSim (sim)")
  dbDisconnect(conn)
  return (list(dbname=db, tablenames=c("allPairwiseSim", "indexMap")))
}

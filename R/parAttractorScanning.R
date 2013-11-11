parAttractorScanning <- function(data, taskList=list(1:nrow(data)), wid=1,  a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE){
  cat("Worker", wid, "initialized.\n");flush.console()
  m <- nrow(data)
  n <- ncol(data)
  genes <- rownames(data)
  task <- taskList[[wid]]
  c <- 1
  as <- NULL
  dd <- NULL
  while(length(task) > 0){
    i <- task[1]
    task <- task[-1]
    cat("Worker ", wid, " : ", genes[i], " ( ", c, " / ", length(task), ")\n", sep="");flush.console()
    out <- CAFrun(data, data[i,], a=a, maxIter=maxIter, epsilon=epsilon, bin=bin, so=so, rankBased=rankBased, negateMI=negateMI, verbose=FALSE, sorting=FALSE)
    if(is.null(out)){
      #if(verbose) {cat("not converged.\n");flush.console()}
      next
    }
    
    killIdx <- which(out >= out[i] & out > 0.5)
    task <- setdiff(task, killIdx)
    
    d <- out[order(out)[m]] - out[order(out)[m-1]]
    if(d <= 0.6){
      if(!is.null(as)){
        un <- apply(as, 1, function(x){
          max(abs(x - out)) > 1E-4
        })
        if(prod(un) == 0){
          rownames(as)[which(un==0)] <- genes[i]
          #if(verbose) {cat("done!\n");flush.console()}
          next
        }
      }
      as <- rbind(as, out)
      rownames(as)[c] <- genes[i]
      c <- c + 1
      #if(verbose) {cat("done!\n");flush.console()}
    }else{
      dd <- c(dd, i)
      #if(verbose) {cat("dominant.\n");flush.console()}
    }
  }
  if(!is.null(as)){
    colnames(as) <- rownames(data)
  }
  oo <- list(attractorMatrix=as, dominant=dd)
  cat("Worker", wid, "finished.\n");flush.console()
  return (oo)
}

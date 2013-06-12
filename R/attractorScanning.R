attractorScanning <- function(data, a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE, outputDominant=FALSE){
  m <- nrow(data)
  n <- ncol(data)
  genes <- rownames(data)
  task <- 1:m
  c <- 1
  as <- NULL
  dd <- rep(FALSE, m)
  while(length(task) > 0){
    i <- task[1]
    cat(genes[i], " ( ", c, " / ", length(task), ") ... ", sep="");flush.console()
    out <- CAFrun(data, data[i,], a=a, maxIter=maxIter, epsilon=epsilon, bin=bin, so=so, rankBased=rankBased,negateMI=negateMI, verbose=FALSE, sorting=FALSE)
    if(is.null(out)){
      cat("not converged.\n");flush.console()
      task <- task[-1]
      next
    }
    
    killIdx <- which(out >= out[i])
    task <- setdiff(task, killIdx)
    
    d <- out[order(out)[m]] - out[order(out)[m-1]]
    if(d <= 0.5){
      if(!is.null(as)){
        un <- apply(as, 1, function(x){
          max(abs(x - out)) > 1E-4
        })
        if(prod(un) == 0){
          killIdx <- which(out >= out[i])
          task <- setdiff(task, killIdx)
          rownames(as)[which(un==0)] <- genes[i]
          cat("done!\n");flush.console()
          next
        }
      }
      as <- rbind(as, out)
      rownames(as)[c] <- genes[i]
      c <- c + 1
      cat("done!\n");flush.console()
    }else{
      dd[i] <- TRUE
      cat("dominant.\n");flush.console()
    }
  }
  if(!is.null(as)){
    colnames(as) <- rownames(data)
  }
  if(outputDominant){
    oo <- list(attractorMatrix=as, dominant=dd)
    return (oo)
  }else{
    return (as)
  }
}

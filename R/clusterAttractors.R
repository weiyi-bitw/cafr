clusterAttractors <- function(filePath="./", fileNames,  numGenes=100, strength.pos=10, min.basins=2,  datasetTags=NULL){
  nf <- length(fileNames)
  if(is.null(datasetTags)){
    datasetTags <- paste("Dataset", sprintf("%03d",1:nf))
  }
  if(length(datasetTags) != nf){
    stop("Length of datasetTags and fileNames must be equal!!!")
  }

  # Load in attractors 
  attractorPool <- list()
  env <- new.env()
  for(fn in 1:length(fileNames)){
    can <- fileNames[fn]
    tag <- datasetTags[fn]
    f <- file.path(filePath, can)
    cat("Processing", f, "...\n");flush.console()
    nm <- load(f, env)[1]
    x <- env[[nm]]
    na <- nrow(x)
    for(i in 1:na){
		o <- order(x[i,], decreasing=TRUE)
		if(min.basins > 0){
		# if the attractor has less than 2 attractees, skip
			if(rownames(x)[i] %in% colnames(x)[o[1:min.basins]]) next
		}
		aid <- paste(tag, sprintf("%03d", i), sep="")
		attractorPool[[aid]] <- Attractor$new(
				id = aid,
				numgenes = numGenes,
				a = x[i,],
				genenames=colnames(x), 
				src=tag,
				qt=strength.pos)
    }
  }

  cat(length(attractorPool), "attractors loaded.\n");flush.console()
  attractorPool <- attractorPool[order(sapply(attractorPool, function(a){a$strength}), decreasing=T)]

  # Calculate all pairwise similarities bewteen attractors
  cat("Caculate all pairwise similarities between attractors and attractor sets...\n");flush.console()
  allPairIdx <- combn(names(attractorPool), 2)
  simList <- apply(allPairIdx, 2, function(pr){
    if(attractorPool[[pr[1]]]$src == attractorPool[[pr[2]]]$src) return (NULL)
    sim <- attractorPool[[pr[1]]]$getOverlapNum( attractorPool[[pr[2]]])
    if(sim < 1) return (NULL)
    return (c(pr[1], pr[2],sim))
  })

  simList <- simList[sapply(simList,function(x){!is.null(x)})]
  o <- order(unlist(sapply(simList, function(x){ as.numeric(as.vector(x[3])) })), decreasing=TRUE)
  simList <- simList[o]

  cnt.clust <- 0

  #clustering attractors
  cat("Clustering attractors...\n");flush.console()
  while(length(simList) > 0){
    p <- simList[[1]]
    cat(p, "\n");flush.console()
    simList <- simList[-1]
    as <- AttractorSet$new(paste("AttractorSet", sprintf("%03d", cnt.clust), sep=""), attractorPool[[p[1]]], nf)
    successMerge <- as$add(attractorPool[[p[2]]])
    if(successMerge){
      attractorPool[[ p[1] ]] <- NULL
      attractorPool[[ p[2] ]] <- NULL
      killIdx <- sapply(simList, function(x){x[1]}) %in% p[1:2] | sapply(simList, function(x){x[2]}) %in% p[1:2]
      if(length(killIdx)>0) simList <- simList[!killIdx]
      addList <- lapply(attractorPool, function(a){
			sim <- as$getOverlapNum( a )
			if(sim < 1) return (NULL)
			return (c(as$id, a$id,sim))
			})
      addList <- addList[sapply(addList,function(x){!is.null(x)})]
      simList <- c(simList, addList)
      if(length(simList) == 0) break
      o <- order(unlist(sapply(simList, function(x){ as.numeric(as.vector(x[3])) })), decreasing=TRUE)
      simList <- simList[o]
		
      attractorPool[[ as$id ]] <- as
      cnt.clust <- cnt.clust + 1
    }
  }

  w <- unlist(lapply(attractorPool, 
		function(x){
			if(class(x)=="Attractor") return (0)
			else if(class(x) == "AttractorSet"){
				if(length(x$attractors) >= 0.75 * length(fileNames) ){
					return(floor(0.75*length(fileNames)))
				}else if(length(x$attractors) >= 0.5 * length(fileNames)){
					return(floor(0.5*length(fileNames))) 
				}else{
					return (length(x$attractors))
				}
			}
		})) + 
	unlist(lapply(attractorPool, 
		function(x){
			if(class(x)=="Attractor") return (x$strength)
			else if(class(x) == "AttractorSet") return (x$medStrength)
		}))
  attractorPool <- attractorPool[order(w, decreasing=TRUE)]
  return (attractorPool)
}


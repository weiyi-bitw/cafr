setwd("/home/weiyi/workspace/data/pancan/")
pancan <- c("BLCA", "HNSC", "LUSC")
#, "LUSC", "LUAD", "OV")
nf <- length(pancan)

# Load in attractors 
attractorPool = list()
for(can in pancan){
	f <- paste("/home/weiyi/workspace/data/pancan/", can, "/attractorMatrix.rda", sep="")
	cat("Processing", f, "...\n");flush.console()
	load(f)
	na <- nrow(x)
	for(i in 1:na){
		o <- order(x[i,], decreasing=TRUE)
		if(rownames(x)[i] %in% colnames(x)[o[1:2]]) next
		aid <- paste(can, sprintf("%03d", i), sep="")
		attractorPool[[aid]] <- Attractor$new(
				id = aid,
				a = x[i,],
				genenames=colnames(x), 
				src=can)
	}
}

cat(length(attractorPool), "attractors loaded.\n");flush.console()

# Calculate all pairwise similarities bewteen attractors
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
		simList <- simList[!killIdx]
		addList <- sapply(attractorPool, function(a){
			sim <- as$getOverlapNum( a )
			if(sim < 1) return (NULL)
			return (c(as$id, a$id,sim))
		})
		addList <- addList[sapply(addList,function(x){!is.null(x)})]
		simList <- c(simList, addList)
		o <- order(unlist(sapply(simList, function(x){ as.numeric(as.vector(x[3])) })), decreasing=TRUE)
		simList <- simList[o]
		
		attractorPool[[ as$id ]] <- as
		cnt.clust <- cnt.clust + 1
	}
}

w <- unlist(lapply(attractorPool, 
		function(x){
			if(class(x)=="Attractor") return (0)
			else if(class(x) == "AttractorSet") return (length(x$attractors))
		})) + 
	unlist(lapply(attractorPool, 
		function(x){
			if(class(x)=="Attractor") return (x$strength)
			else if(class(x) == "AttractorSet") return (x$minStrength)
		}))
attractorPool <- attractorPool[order(w, decreasing=TRUE)]


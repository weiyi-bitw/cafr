attractorScanningGL = function(data, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE, saveAllScore=FALSE){
	genes.genome = rownames(genome)
	data = data[intersect(genes.genome, rownames(data)),] # sort the rownames according to genomic location
	m = nrow(data)
	n = ncol(data)

	genes.data = rownames(data)
	mg = nrow(genome)
	alist = list()
	ascore = rep(0, m)
	aalpha = rep(NA, m)
	for(i in 1:m){
		if(verbose & i %% 100 == 0){
			cat(i , "/", m, "\n"); flush.console()
		}
		idx.center = which(genes.genome == genes.data[i])
		winspan = floor(windowSize/2)
		wd = genes.genome[max(1, idx.center - winspan) : min(mg, idx.center+winspan)]
		wd = intersect(wd, genes.data)
		#print(wd)
		if(length(wd) < num.output){
			alist[[i]] = NA
			next
		}
		dataIn = data[wd,]
		mi = getAllMIWz(dataIn, dataIn[genes.data[i],], bin=bin, so=so, negateMI=negateMI)
		premi = mi
		for(a in alpha){
			w = abs(mi)^a / sum(abs(mi)^a)
			w[mi < 0] = 0
			metagene = w %*% dataIn
			c = 0
			while(c < maxIter){
				mi = getAllMIWz(dataIn, metagene, bin=bin, so=so, negateMI=negateMI)
				delta = sum((mi - premi)^2)
				if(delta < epsilon){
					break
				}
				premi = mi
				mi[mi < 0] = 0
				metagene = w %*% dataIn
				c = c + 1
			}
			if (c >= maxIter){ 
				alist[[i]] = NA
			}else{
				aout = sort(mi, decreasing=T)[1:num.output]
				score = aout[score.position]
				if(score >= ascore[i]){
					alist[[i]] = aout
					ascore[i] = score
					aalpha[i] = a
				}
			}
		}
	}
	killidx = which(is.na(alist))
	if (length(killidx) > 0){
		alist[killidx] = NULL
		ascore = ascore[-killidx]
		aalpha = aalpha[-killidx]
		aglist = sapply(alist, names)
		colnames(aglist) = genes.data[-killidx]
		if (saveAllScore){
			 alist = simplify2array(alist)
			colnames(alist) = genes.data[-killidx]
		}else{
			alist=NULL
		}
	}else{
		aglist = sapply(alist, names)
		colnames(aglist) = genes.data
		if (saveAllScore){
			 alist = simplify2array(alist)
			colnames(alist) = genes.data
		}else{
			alist=NULL
		}

	}	
	cat("Summarizing output...");flush.console()
	rownames(alist) = NULL
	out = list(attractome=aglist, score=ascore, bestAlphas = aalpha, scoremat = alist)
	sumOut = summarizeAttractorScanningGL(out, genes.genome)
	if(saveAllScore){
		out$summary=sumOut
	}else{
		out = sumOut
	}
	cat("done.\n");flush.console()
	return (out)
}

summarizeAttractorScanningGL = function(out, genes.genome, windowSize=50){
	topIdx = sapply(out$attractome[1,], function(x){which(genes.genome==x)})
	o = order(out$score, decreasing=T)
	sumOut = NULL
	center = NULL
	while(length(o) > 0){
		o.top = o[1]
		#print(length(o))
		#print(o.top)
		sumOut = rbind(sumOut,c(out$attractome[,o[1]], out$score[o[1]]))
		center = c(center, colnames(out$attractome)[o[1]])
		idxDiff = abs(topIdx - topIdx[o.top])
		killIdx = which(idxDiff <= (windowSize))
		if(length(killIdx) == 0) break
		o = setdiff(o, killIdx)
	}
	rownames(sumOut) = center
	return (sumOut)
}

findCytoband = function(sumOut, genome){
	genes.genome = rownames(genome)
	cytoband = apply(sumOut, 1, 
	function(a, genes.genome, genome){
		genes = a[1:(length(a)-1)]
		idx = which(genes.genome %in% genes)
		cys = genome[idx, "Cytoband"]
		cys = cys[!is.na(cys)]
		start = cys[1]
		end = cys[length(cys)]
		if(start == end){
			cy = paste(genome[idx[1],"Chr"], start, sep="")
		}else{
			cy = paste(genome[idx[1],"Chr"], start, "-", end, sep="")
		}
		return(cy)
	}, genes.genome=genes.genome, genome=genome
	)
	out = cbind(cytoband, sumOut)
	return (out)
}
findAttractor = function(data, vec, a=5, maxIter = 100, epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  negateMI = TRUE, verbose=TRUE){
	dataIn = data
	if(rankBased){
		vec = rank(vec)
		for(i in 1:nrow(data)){
			dataIn[i,] = rank(data[i,])
		} 
	}
	mi = getAllMIWz(dataIn, vec, bin=bin, so=so, negateMI = negateMI)
	premi = mi
	w = abs(mi)^a / sum(abs(mi)^a)
	w[mi < 0] = 0
	metagene = w %*% data
	if(rankBased) metagene = rank(metagene)
	c = 0
	while(c < maxIter){
		mi = getAllMIWz(dataIn, metagene, bin=bin, so=so, negateMI = negateMI)
		delta = sum((mi - premi)^2)
		if(verbose){
			cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
			print(mi[order(mi, decreasing=T)[1:20]])
			flush.console()
		}
		if(delta < epsilon){
			break
		}
		premi = mi
		mi[mi < 0] = 0
		w = abs(mi)^a / sum(abs(mi)^a)
		w[mi < 0] = 0
		metagene = w %*% data
		if(rankBased) metagene = rank(metagene)
		c = c + 1
	}
	if(c >= maxIter) return (NULL)
	return (sort(mi, decreasing=T))
}

findGLAttractor = function(data,seed, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE){
	if(! seed %in% rownames(genome)) stop("Cannot find seed gene in genome file rownames!")
	genes.genome = rownames(genome)
	idx.seed = which(rownames(genome)==seed)
	data = data[intersect(genes.genome, rownames(data)),] # sort the rownames according to genomic location
	n = ncol(data)
	m = nrow(data)
	mg = length(genes.genome)

	idxrange = max(1, idx.seed-windowSize):min(mg, idx.seed + windowSize) 
	generange = genes.genome[idxrange]
	data = data[intersect(generange, rownames(data)),]
	
	out = attractorScanningGL(data, genome, alpha, windowSize, maxIter, epsilon, bin, so, score.position, num.output, negateMI, verbose, saveAllScore=TRUE)
	sumOut = out$summary
	if(nrow(sumOut) > 1){
		idxtop = sapply(sumOut[,1], function(g){which(genes.genome==g)})
		select = which.min(abs(idxtop - idx.seed))
	}else{
		select = 1
	}
	genes = sumOut[select,1:num.output]
	mis = out$scoremat[,rownames(sumOut)[select]]

	names(mis) = genes

	return(mis)

}


CAFrun = function(data, vec, a=5, maxIter = 100, epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  negateMI = TRUE, verbose=TRUE, sorting=TRUE){
	m = nrow(data)
	n = ncol(data)

	if(rankBased){
		vec = rank(vec)
		dataIn = t( apply(data, 1, rank) )
	}
	miOut = .Call("cafR2C", data, vec, a, maxIter, epsilon, m, n, bin, so, as.integer(negateMI), as.integer(verbose))
	
	if(miOut[1] == -999) return (NULL)

	names(miOut) = rownames(data)
	if(sorting){
		return (sort(miOut, decreasing=T))
	}else{
		return (miOut)
	}
}

findCorrAttractor = function(data, vec, a=10,rankBased = FALSE, maxIter=100, epsilon=1E-14){
	dataIn = data
	if(rankBased){
		vec = rank(vec)
		for(i in 1:nrow(data)){
			dataIn[i,] = rank(data[i,])
		}
	}
	rs = getAllCorWz(dataIn, vec, rankBased=rankBased)
	prers = rs
	w = abs(rs)^a / sum(abs(rs)^a)
	w[rs < 0] = -w[rs < 0]
	metagene = w %*% data
	if(rankBased) metagene = rank(metagene)
	c = 0
	while(c < maxIter){
		rs = getAllCorWz(dataIn, metagene)
		delta = sum((rs - prers)^2)
		cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
		print(rs[order(rs, decreasing=T)[1:20]])
		flush.console()
		if(delta < epsilon){
			break
		}
		prers = rs
		w = abs(rs)^a / sum(abs(rs)^a)
		w[rs < 0] = -w[rs < 0]
		metagene = w %*% data
		if(rankBased) metagene = rank(metagene)
		c = c + 1
	}
	if(c >= maxIter) return (NULL)
	return (sort(rs, decreasing=T))
}

attractorScanning = function(data, a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE){
	m = nrow(data)
	n = ncol(data)
	genes = rownames(data)
	task = 1:m
	c = 1
	as = NULL
	while(length(task) > 0){
		i = task[1]
		cat(genes[i], " ( ", c, " / ", length(task), ") ... ", sep="");flush.console()
		out = CAFrun(data, data[i,], verbose=FALSE, sorting=FALSE)
		if(is.null(out)){
			cat("not converged.\n");flush.console()
			task = task[-1]
			next
		}
		
		killIdx = which(out >= out[i])
		task = setdiff(task, killIdx)

		d = out[order(out)[m]] - out[order(out)[m-1]]
		if(d <= 0.5){
			if(!is.null(as)){
				un = apply(as, 1, function(x){
					max(abs(x - out)) > 1E-4
				})
				if(prod(un) == 0){
					killIdx = which(out >= out[i])
					task = setdiff(task, killIdx)
					rownames(as)[which(un==0)] = genes[i]
					cat("done!\n");flush.console()
					next
				}
			}
			as = rbind(as, out)
			rownames(as)[c] = genes[i]
			c = c + 1
			cat("done!\n");flush.console()
		}else{
			cat("dominant.\n");flush.console()
		}
	}
	if(!is.null(as)){
		colnames(as) = rownames(data)
	}
	return (as)
}

parAttractorScanning = function(data, taskList=list(1:nrow(data)), wid=1,  a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE){
	cat("Worker", wid, "initialized.\n");flush.console()
	m = nrow(data)
	n = ncol(data)
	genes = rownames(data)
	task = taskList[[wid]]
	c = 1
	as = NULL
	while(length(task) > 0){
		i = task[1]
		cat("Worker ", wid, " : ", genes[i], " ( ", c, " / ", length(task), ")\n", sep="");flush.console()
		out = CAFrun(data, data[i,], verbose=FALSE, sorting=FALSE)
		if(is.null(out)){
			#if(verbose) {cat("not converged.\n");flush.console()}
			task = task[-1]
			next
		}
		
		killIdx = which(out >= out[i])
		task = setdiff(task, killIdx)

		d = out[order(out)[m]] - out[order(out)[m-1]]
		if(d <= 0.5){
			if(!is.null(as)){
				un = apply(as, 1, function(x){
					max(abs(x - out)) > 1E-4
				})
				if(prod(un) == 0){
					killIdx = which(out >= out[i])
					task = setdiff(task, killIdx)
					rownames(as)[which(un==0)] = genes[i]
					#if(verbose) {cat("done!\n");flush.console()}
					next
				}
			}
			as = rbind(as, out)
			rownames(as)[c] = genes[i]
			c = c + 1
			#if(verbose) {cat("done!\n");flush.console()}
		}else{
			#if(verbose) {cat("dominant.\n");flush.console()}
		}
	}
	if(!is.null(as)){
		colnames(as) = rownames(data)
	}
	cat("Worker", wid, "finished.\n");flush.console()
	return (as)
}


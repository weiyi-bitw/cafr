findAttractor = function(data, vec, a=5, maxIter = 100, epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  negateMI = TRUE){
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
		cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
		print(mi[order(mi, decreasing=T)[1:20]])
		flush.console()
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
	if(c >= maxIter) return (NA)
	return (sort(mi, decreasing=T))
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
	if(c >= maxIter) return (NA)
	return (sort(rs, decreasing=T))
}


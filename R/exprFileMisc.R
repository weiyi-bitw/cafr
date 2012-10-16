# use built-in R function, slow, use load.exp instead
load.exp <- function(file, sep='\t'){
	line = readLines(file)
	tokens = strsplit(line[1], "\t")[[1]]
	n = length(tokens)-1
	m = length(line)-1
	mset = matrix(NA, m, n)	
	cname = tokens[2:(n+1)]
	rname = rep(NA, m)
	b = txtProgressBar(style=3)
	for(i in 1:m){
		tokens = strsplit(line[i+1], "\t")[[1]]
		tokens[tokens=="NA"] = NA		
		rname[i] = tokens[1]
		mset[i,] = as.numeric(tokens[2:(n+1)])
		if(i %% 100 == 0){
			setTxtProgressBar(b, i/m)
		}
	}
	cat("\n")
	colnames(mset) <- cname
	rownames(mset) <- rname
	mset
}

# load text mat, such as clinical file
loadTextMat <- function(file, sep="\t"){
	mset <- read.table(file=file, sep=sep, as.is=TRUE)
	col.names <- mset[1, -1]
	row.names <- mset[-1, 1]
	mset <- mset[-1,-1]
	mset <- as.matrix(mset)
	colnames(mset) <- col.names
	rownames(mset) <- row.names
	mset
}

loadGs <- function(file, rownames, colnames, sep="\t", header=T){
	gs <- as.matrix(read.table(file, sep=sep, header=header))
	n <- dim(gs)[1]
	B <- matrix(0, nrow=length(rownames), ncol=length(colnames))
	rownames(B) <- rownames
	colnames(B) <- colnames
	for(i in 1:n){
		tf <- gs[i,1]
		tg <- gs[i,2]
		if((tf %in% colnames) & (tg %in% rownames)){
			B[tg, tf] <- 1
		}
	}
	return (B)
}

varFilter <- function(mset,dim=1, numFilterIn=10000){
	if(dim==1){
		mset <- t(mset)
	}
	mset.sd <- sd(mset)
	idx <- sort(mset.sd,decreasing=TRUE, index.return=TRUE)$ix[1:numFilterIn]
	mset <- mset[,sort(idx)]
	if(dim==1){
		mset <- t(mset)
	}
	mset
}

cbs <- function(ge, map, sumfun=median, corTh = 0.9){
  genes <- sort(unique(map))
  geids <- rownames(ge)
  mapids <- rownames(map)
  m <- length(genes)
  n <- dim(ge)[2]
  out = NULL
  cnt <- 0
  glist <- list()
  for(g in genes){
    cnt = cnt + 1
    ps = rownames(map)[which(map==g)]
    nps = length(ps)
    if(nps==1){
      out = rbind(out, ge[ps,])
      glist[cnt] = g
      plist[cnt] = ps
    }else{
      #cc = 
      
    }
  }

}


summarizeGenes <- function(ge, map, sumfun=median){
	genes <- sort(unique(map))
	#genes = genes[-which(genes==""),1]
	geids <- rownames(ge)
	mapids <- rownames(map)
	m <- length(genes)
	n <- dim(ge)[2]
	out <- matrix(0, nrow=m, ncol=n)	
	cnt <- 0
	glist <- list()	
	for(g in genes){
		ids <- mapids[which(map==g)]
		idx <- which(geids %in% ids)
		if (length(idx) > 0){
			cnt <- cnt + 1
			glist[cnt] <- g
			# cat(idx,"\n")
			# flush.console()
			if(length(idx)==1){
				out[cnt,] <- ge[idx,]
			}else{
				out[cnt,] <- apply(ge[idx,],2,sumfun)
			}
		}
	}
	out <- out[1:cnt,]
	rownames(out) <- unlist(glist)
	colnames(out) <- colnames(ge)
	return (out)
}

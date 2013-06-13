# R script for merging the results generated from Amazon EC2 cluster
# Wei-Yi Cheng
# 2013.03.07
#
# arguments for script:
#	args[1] : directory of the results (contains exclusively the .rda files)
# Example
# > Rscript mergeResults.R resultDir/
#
#

args = commandArgs(TRUE)

filePath = "./"
if(length(args) > 0){
	filePath = args[1]
}

fileNames = list.files(pattern="*.attractors.*.rda", path = filePath)
nf = length(fileNames)

x = NULL
dom = NULL
for(i in 1:nf){
	cat("Processing", fileNames[i], "...\n");flush.console()
	load(file.path(filePath,fileNames[i]))
	if(is.null(as)) next
	if(is.null(x)){
		x = as
		next
	}
	dom = c(dom, dd)
	m = nrow(x)
	ma = nrow(as)
	attachIdx = NULL
	mergeIdx = NA
#	for(i in 1:ma){
#		a = as[i,]
#		for(j in 1:m){
#			xx = x[j,]
#			dd = max(a-xx)
#			if(dd < 1E-4){
#				mergeIdx = j
#				break
#			}
#		}
#		if(!is.na(mergeIdx)) break
#	}

	mergeIdx = apply(as, 1, function(a){
		dd = apply(x, 1, function(xx){max(a - xx)})
		if(min(dd) < 1E-4) return (which.min(dd))
		else return (NA)
	})
	for(j in 1:ma){
		if(is.na(mergeIdx[j])){ # no same attractor, attach to the attractor matrix
			x = rbind(x, as[j,])
			rownames(x)[nrow(x)] = rownames(as)[j]
		}else{ # there is same attractor, check if the last seed need to be updated
			lastSeed = rownames(x)[mergeIdx[j]]
			asLastSeed = rownames(as)[j]
			if(x[mergeIdx[j],asLastSeed] > as[j, asLastSeed]){
				rownames(x)[mergeIdx[j]] = asLastSeed
			}
		}
	}
}
gc()


strength = apply(x, 1, function(xx){sort(xx, decreasing=T)[10]})
x = x[order(strength, decreasing=T),]
m = nrow(x)
ma = ncol(x)
dd = rep(FALSE, ma)
names(dd) = colnames(x)
dd[dom] = TRUE
  if(!is.null(x)){
    killIdx2 <- rep(FALSE, m )
    for(i in 1:m){
  	if(max(x[i,]) == x[i,rownames(x)[i]]){
		dd[rownames(x)[i]] <- TRUE
		killIdx2[i] <- TRUE
	}
    }
    if(sum(killIdx2) == nrow(x)-1){
      bkup <- rownames(x)[!killIdx2]
      x <- t(x[!killIdx2,])
      rownames(x) <- bkup
    }else{
      x <- x[!killIdx2,]
    }
  }

cat(m, "attractors in total.\n\n")

cat("Generate summarized files...\n");flush.console()

out1 = NULL
k = 1
for(i in 1:m){
	o = order(x[i,], decreasing=TRUE)
	if(rownames(x)[i] %in% colnames(x)[o[1:2]]) next
	out1 = cbind(out1,colnames(x)[o[1:100]], round(x[i,o[1:100]], 4))
	aname = paste("Attractor", sprintf("%03d",k), sep="")
	colnames(out1)[c(2*k-1, 2*k)] = c(paste(aname, ":Genes",sep="" ), paste(aname, ":MI",sep="") )
	k = k+1
}
rownames(out1) = 1:100

save(x, file="attractorMatrix.rda")
write.table(out1, file="attractors.txt", sep='\t', quote=F)



probeSummarization = function(ge, map, method="corr", threshold=0.5, gene.colname="Gene.Symbol", verbose=TRUE){
	m = nrow(ge)
	n = ncol(ge)

	#check for missing value
	mis = sum(is.na(ge))
	if(mis > 0){
		cat("There are missing values in the dataset!!\n")
		cat("Using mean imputation?\n\t 'y' for mean imputation, 'n' for quit (y/N): ");
		ans <- readline()
		if(ans == "y" | ans == "Y"){
			ge = t(apply(ge, 1, function(x){ idx = which(is.na(x)); x[idx] = mean(x, na.rm=T); x}  ))
			#print(dim(ge))
		}else{
			cat("Quit.\n");
			return (NULL)
		}
	}

	mappedGenes = rep(NA, m)
	mappedGenes[rownames(ge) %in% rownames(map)] = map[intersect(rownames(ge), rownames(map)), gene.colname]
	ugenes = sort(unique(mappedGenes[!is.na(mappedGenes)]))
	mg = length(ugenes)
	groupMap = 1:mg
	names(groupMap) = ugenes
	grp = groupMap[mappedGenes] - 1
	grp[is.na(grp)] = -1

	if(method=="mi"){
		useCorr = FALSE
	}else{
		useCorr = TRUE
	}

	out = .Call("probe_summarizationR2C", ge, m, n, grp, mg, useCorr, threshold, as.integer(verbose))
	out = matrix(out, nrow=mg, ncol=n, dimnames=list(ugenes, colnames(ge)))

	cat("Filtering no-probe genes...\n");flush.console()
	sds = apply(out, 1, sd)
	killIdx = which(sds == 0)
	out = out[-killIdx,]

	cat("Done.\n")
	return (out)
}



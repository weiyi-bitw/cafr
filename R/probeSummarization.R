
probeSummarization = function(ge, map, method="corr", threshold=0.6, map.name="Gene.Symbol"){
	m = nrow(ge)
	n = ncol(ge)
	mappedGenes = rep(NA, m)
	mappedGenes[rownames(ge) %in% rownames(map)] = map[intersect(rownames(ge), rownames(map)), map.name]
	ugenes = sort(unique(mappedGenes[!is.na(mappedGenes)]))
	mg = length(ugenes)
	groupMap = 1:mg
	names(groupMap) = ugenes
	grp = groupMap[mappedGenes]

	if(method=="mi"){
		useCorr = FALSE
	}else{
		useCorr = TRUE
	}


	out = .Call("probe_summarizationR2C", ge, m, n, grp, mg, useCorr, threshold)
	return (out)
}



#' AttractorSet
#'
#' An object containing a set of Attractor
#'
#' @author Wei-Yi Cheng
#' @export


AttractorSet <- setRefClass(
	Class="AttractorSet",
	field=list(id="character", 
		attractors = "list",
		capacity = "numeric",
		minStrength = "numeric"),
	methods=list(
		initialize = function(id, a, k){
'Initialize an AttractorSet object\n
Arguments:\n\tid : ID of the attractor set\n
\ta : The initial attractor in the attractor set\n
\tk : the capacity of the attractor set, usually the number of datasets from which the attractors are generated.'
			.self$id <- id
			capacity <<- k
			attractors <<- list()
			.self$add(a)
			if(class(a) == "Attractor"){
				minStrength <<- a$strength
			}else if(class(a) == "AttractorSet"){
				minStrength <<- a$minStrength
			}
			return(.self)
		},
		add = function(a){
'Add another attractor into the attractor set. If the argument is an Attractor, and the attractor set already contains an attractor from the same source, the method returns FALSE and the attractor will not be added. If the argument is AttractorSet and more than one-third of its attractors are from the same source, the method returns FALSE and the two attractor sets will not be merged. Otherwise the method returns TRUE and the Attractor will be added, or the AttractorSet will be merged by choosing from the Attractors from the same source the one with higher strength.\n
Arguments:\n
\ta : An Attractor or AttractorSet object.'
			if(class(a) == "Attractor"){
				if(!is.null(.self$attractors[[a$src]])) return (FALSE)
				attractors[[a$src]]<<- a
				minStrength <<- min(.self$minStrength, a$strength)
				return (TRUE)
			}else if(class(a)=="AttractorSet"){
				ovlpSrc <- sum(unlist(lapply(a$attractors, function(aa){!is.null(.self$attractors[[aa$src]])}) ))
				if(ovlpSrc > (.self$capacity / 3) ) return (FALSE)
				for(aa in a$attractors){
					if(is.null(.self$attractors[[aa$src]])){
						attractors[[aa$src]] <<- aa
					}else{
						myaa <- .self$attractors[[aa$src]]
						if(myaa$strength < aa$strength){
							attractors[[aa$src]] <<- aa
						}
					}
				}
				minStrength <<- min( unlist(lapply(.self$attractors, function(aa){aa$strength})) )
				return (TRUE)
			}
		},
		getOverlapNum = function(a){
'Returns the number of overlapping genes with another Attractor or AttractorSet object.\n 
Arguments:\n
\ta : An object of Attractor or Attractor Set.'
			if(class(a)=="Attractor"){
				allgenes <- names(a$genes)
				allgenes <- c(allgenes, unlist(lapply(.self$attractors, 
						function(aa){
							if(aa$src == a$src){ return (NULL) }
							else return (names(aa$genes))
						}))
						)
				t <- table(allgenes)[names(a$genes)]
				return (sum(t[t>1]))
			}else if(class(a) == "AttractorSet"){
				k <- lapply(.self$attractors, function(aa){aa$getOverlapNum(a)})
				return (sum(unlist(k)))
			}
		},
		getGeneTable = function(...){
'Returns a vector of all genes in the attractor set ranked according to their occurrences in the attractor set.\n'
			allgenes <- unlist(lapply(.self$attractors, function(aa){names(aa$genes)}))
			return (sort(table(allgenes), decreasing=TRUE))
		},
		getConsensus = function(sz=50){ # NOT REAL CONSENSUS, ONLY THE CONSENSUS FROM TOP GENES
'Returns a vector of genes and their MIs of size sz according to their average MI across the attractors in the attractor set.\n
Arguments:\n
\tsz : Number of top genes in the output.\n
NOTE : To produce more accurate consensus, the choice of sz should be much less than the number of genes in the Attractor. '
			allgenes <- unlist(lapply(.self$attractors, function(aa){names(aa$genes)}))
			g <- unique(allgenes)
			mat <- matrix(0, nrow=length(g), ncol=.self$capacity)
			rownames(mat) <- g
			colnames(mat) <- unlist(lapply(.self$attractors, function(aa){aa$src}))
			for(aa in .self$attractors){
				gg <- names(aa$genes)
				mat[gg, aa$src] <- aa$genes
			}
			mis <- apply(mat, 1, mean)
			names(mis) <- g
			return (sort(mis, decreasing=TRUE)[1:sz])
		},
		getGeneMatrix = function(sz=50){
'Return a table of genes from each Attractor in the AttractorSet.\n
Arguments:\n
\t sz : Number of top genes from each Attractor'
			allgenes <- matrix(unlist(lapply(.self$attractors, function(aa){names(aa$genes)})),nrow=length(.self$attractors), byrow=TRUE)[,1:sz]
			rownames(allgenes) <- unlist(lapply(.self$attractors, function(aa){aa$id}))
			allgenes <- cbind(allgenes, unlist(lapply(.self$attractors, function(aa){aa$strength})))
			return (allgenes)
		}
		
	)

)

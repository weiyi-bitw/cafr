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
			allgenes <- unlist(lapply(.self$attractors, function(aa){names(aa$genes)}))
			return (sort(table(allgenes), decreasing=TRUE))
		},
		getConsensus = function(sz=50){ # NOT REAL CONSENSUS, ONLY THE CONSENSUS FROM TOP GENES
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
			allgenes <- matrix(unlist(lapply(.self$attractors, function(aa){names(aa$genes)})),nrow=length(.self$attractors), byrow=TRUE)[,1:sz]
			rownames(allgenes) <- unlist(lapply(.self$attractors, function(aa){aa$id}))
			allgenes <- cbind(allgenes, unlist(lapply(.self$attractors, function(aa){aa$strength})))
			return (allgenes)
		}
		
	)

)

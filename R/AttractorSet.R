#' AttractorSet
#'
#' An object containing a set of Attractor
#'
#' @author Wei-Yi Cheng
#' @export


AttractorSet <- setRefClass(
	Class="AttractorSet",
	field=list(id="character", 
		capacity = "numeric",
		attractors = "list"),
	methods=list(
		initialize = function(id, a, k){
			.self$id <- id
			capacity <<- k
			attractors <<- list()
			.self$add(a)
			return(.self)
		},
		add = function(a){
			if(class(a) == "Attractor"){
				attractors[[length(attractors)+1]]<<- a
				names(attractors)[[length(attractors)]] <<- a$src
			}else if(class(a)=="AttractorSet"){
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
				t <- table(allgenes)
				return (sum(t[t>1]))
			}else if(class(a) == "AttractorSet"){
				k <- lapply(.self$attractors, function(aa){aa$getOverlapNum(a)})
				return (sum(unlist(k)))
			}
		},
		getGeneTable = function(...){
			allgenes <- unlist(lapply(.self$attractors), function(aa){names(aa$genes)})
			return (sort(table(allgenes), decreasing=TRUE))
		},
		getConsensus = function(sz){ # NOT REAL CONSENSUS, ONLY THE CONSENSUS FROM TOP GENES
			allgenes <- unlist(lapply(.self$attractors), function(aa){names(aa$genes)})
			g <- unique(allgenes)
			mat <- matrix(0, nrow=length(g), ncol=.self$capacity)
			rownames(mat) <- g
			colnames(mat) <- unlist(lapply(.self$attractors, function(aa){aa$src}))
			for(aa in .self$attractors){
				gg <- names(aa$genes)
				mat[gg, aa$src] <- aa$genes
			}
			mis <- apply(mat, 1, mean)
			return (sort(mat, decreasing=TRUE)[1:sz])
		}
		
	)

)

#' Attractor
#' 
#' An object containing top genes and their MI with the attractor
#'
#' @author Wei-Yi Cheng
#' @export


Attractor <- setRefClass(
	Class="Attractor",
	fields=list(
		id = "character",
		genes = "numeric",
		strength = "numeric",
		src = "character"),
	methods=list(
		initialize = function(id, a, genenames, src, numgenes=100, qt=10){
			.self$id <- id
			o <- order(a, decreasing=TRUE)[1:numgenes]
			genes <<- a[o]
			names(genes) <<- genenames[o]
			strength <<- genes[qt]
			src <<- src
			return (.self)
		},
		getOverlapNum = function(a){
			if(class(a) == "Attractor"){
				t <- table(c(names(.self$genes), names(a$genes)))
				return (sum(t[t>1]))
			}else if(class(a) == "AttractorSet"){
				allgenes <- names(.self$genes)
				allgenes <- c(allgenes, unlist(lapply(a$attractors, 
							function(aa){
								if(aa$src == .self$src){ return (NULL) }
								else return (names(aa$genes))
							}))
					)
				t <- table(allgenes)[names(.self$genes)]
				return (sum(t[t>1]))
			}
		}
	)
)

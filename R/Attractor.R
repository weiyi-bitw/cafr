#' Attractor
#' 
#' An object containing top genes and their MI with the attractor
#'
#' @author Wei-Yi Cheng
#' @export


Attractor <- setRefClass(
	Class="Attractor",
	fields=list(genes = "numeric",
		src = "character"),
	methods=list(
		initialize = function(a, genenames, src, numgenes){
			o <- order(a, decreasing=TRUE)[1:numgenes]
			genes <<- a[o]
			names(genes) <<- genenames[o]
			src <<- src
			return (.self)
		},
		getOverlapNum = function(a){
			if(class(a) == "Attractor"){
				t <- table(c(names(.self$genes), names(a$genes)))
				return (sum(t[t>1]))
			}else if(class(a) == "AttractorSet"){
				
			}
		}
	)
)

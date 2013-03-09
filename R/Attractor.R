#' Attractor
#' 
#' An object containing top genes and their MI with the attractor
#'
#' @author Wei-Yi Cheng
#' @export


Attractor <- setRefClass(Class="Attractor",
			fields=list(genes = "character",
				mis = "numeric", 
				src = "character"),
			methods=list(
				initialize = function(a, genenames, src, numgenes, qt){
					o <- order(a, decreasing=TRUE)[1:numgenes]
					mis <<- a[o]
					genes <<- genenames[o]
					src <<- src
					return (.self)
				}
			)
)

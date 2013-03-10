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
			attractors[[length(attractors)+1]]<<- a
			names(attractors)[[length(attractors)]] = a$src
		}
	)

)

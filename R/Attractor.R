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
			## Self-documentation
'Initialize an Attractor object\n
Arguments:\n\tid : ID of the attractor\n
\ta : A vector of MIs of the genes in the attractor. As output by the findAttractor function\n
\tgenenames : A vector of all gene symboles in the dataset.\n
\tsrc : The source of the attractor, usually the name of the dataset it is generated.\n
\tnumGenes : Number of top genes stored in the object.\n
\tqt : Which rank should be used to represent the strength of the attractor. Default is the 10th highest MI.'
			.self$id <- id
			o <- order(a, decreasing=TRUE)[1:numgenes]
			genes <<- a[o]
			names(genes) <<- genenames[o]
			strength <<- genes[qt]
			src <<- src
			return (.self)
		},
		getOverlapNum = function(a){
			## Self-documentation of the method
'Returns the number of overlapping genes with another Attractor or AttractorSet object.\n 
Arguments:\n
\ta : An object of Attractor or Attractor Set.'
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

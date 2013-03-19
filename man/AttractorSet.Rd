\name{AttractorSet}
\docType{class}
\alias{AttractorSet}

\title{Class \code{"AttractorSet"}}
\description{
A reference class for an attractor set. This class can be used to store attractors from difference
sources with overlapping genes. It contains methods to give summarized results such as consensus
attractors or a ranked gene list according to their occurrence in different datasets.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
}
\author{
Wei-Yi Cheng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{Attractor}}, \code{\link{clusterAttractors}}
}
\examples{
showClass("AttractorSet")

## ======== Find the CIN attractor in the ov dataset =======
# Load the toy dataset extracted from TCGA OV data
data(ov)

# find attractor using CENPA as a seed
x <- ov["CENPA",]
a <- findAttractor(ov, x)

attractor1 <- Attractor$new("OV001", a, genenames=names(a), src="OV", numgenes=100, qt=10)

## ======== Find the CIN attractor in the brca.pbs dataset======
# load Toy version of Wang et al. breast cancer dataset (GSE2034)
data(brca.pbs)

# download the HGU133A 2.0 annotations
source("http://bioconductor.org/biocLite.R")
biocLite("hgu133a2.db")
library(hgu133a2.db)

# Create map object to fit the format
x <- hgu133a2SYMBOL
map <- cbind(unlist(as.list(x[mappedkeys(x)])))
colnames(map) <- "Gene.Symbol"

# summarize into gene-level expression after eliminating uncorrelated probes
brca <- probeSummarization(brca.pbs, map)

x <- brca["CENPA",]
a <- findAttractor(brca, x)

attractor2 <- Attractor$new("BRCA001", a, genenames=names(a), src="BRCA", numgenes=100, qt=10)

# Find the number of overlapping genes between the two attractors
attractor1$getOverlapNum(attractor2)

# Initialize an AttractorSet object using the OV CIN attractor
attractorSet1 <- AttractorSet$new(id="AttractorSet001", a=attractor1, k=2)
# add the BRCA CIN attractor
attractorSet1$add(attractor2)
# rank the genes in the attractor set according to their occurrence
head(attractorSet1$getGeneTable())
# create a consensus attractor using the two CIN attractors
attractorSet1$getConsensus()
# shows the matrix of top genes in each attractor in the attractor set
attractorSet1$getGeneMatrix()
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{id}:}{ID of the attractor set. }
    \item{\code{attractors}:}{A list of Attractor in the attractor set.}
    \item{\code{capacity}:}{Capacity of the AttractorSet, usually the number of datasets from which the attractors are generated.}
    \item{\code{minStrength}:}{The minimum strength among the Attractors in the set.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{getGeneMatrix(sz)}:}{
Return a table of genes from each Attractor in the AttractorSet.

Arguments:

  \code{sz} : Number of top genes from each Attractor
}
    \item{\code{getConsensus(sz)}:}{
Returns a vector of genes and their MIs of size sz according to their average MI across the attractors in the attractor set.

Arguments:

  \code{sz} : Number of top genes in the output.

NOTE : To produce more accurate consensus, the choice of sz should be much less than the number of genes in the Attractor.
}
    \item{\code{getGeneTable(...)}:}{
Returns a vector of all genes in the attractor set ranked according to their occurrences in the attractor set.
}
    \item{\code{getOverlapNum(a)}:}{
Returns the number of overlapping genes with another Attractor or 
AttractorSet object.

Arguments:

  \code{a} : An object of Attractor or Attractor Set.
}
    \item{\code{add(a)}:}{
Add another attractor into the attractor set. If the argument is an 
Attractor, and the attractor set already contains an attractor from 
the same source, the method returns \code{FALSE} and the attractor 
will not be added. If the argument is AttractorSet and more than 
one-third of its attractors are from the same source, the method returns 
\code{FALSE} and the two attractor sets will not be merged. Otherwise 
the method returns \code{TRUE} and the Attractor will be added, or the 
AttractorSet will be merged by choosing from the Attractors from the same
 source the one with higher strength.

Arguments:

  \code{a} : An Attractor or AttractorSet object.
    }
    \item{\code{initialize(id, a, k)}:}{
Initialize an AttractorSet object

Arguments:

  \code{id} : ID of the attractor set

  \code{a} : The initial attractor in the attractor set

  \code{k} : the capacity of the attractor set, usually the number of datasets 
  from which the attractors are generated.
    }
  }
}

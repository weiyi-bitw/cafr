\name{Attractor}
\docType{class}
\alias{Attractor}

\title{Class \code{"Attractor"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
A reference class for an attractor object. This object is used to summarize attractors found
in difference datasets. And also a basic building block for class \code{AttractorSet}.
It contains top genes in an attractor, and a method \code{getOverlapNum} to return the number 
of overlapping genes with other attractors. 
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

\seealso{
\code{\link{AttractorSet}}, \code{\link{findAttractor}}
}
\examples{
showClass("Attractor")

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
    \item{\code{id}:}{ID of the Attractor.}
    \item{\code{genes}:}{A vector of top genes and their MIs in the attractor.}
    \item{\code{strength}:}{The strength of the attractor, defined by the \code{qt}-th highest MI in the attractor as assigned when initialized.}
    \item{\code{src}:}{The source of the attractor, usually the name of the dataset from which the attractor is generated.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{getOverlapNum(a)}:}{
Returns the number of overlapping genes with another Attractor or AttractorSet object.

Arguments:

  \code{a} : An object of Attractor or Attractor Set.
    }
    \item{\code{initialize(id, a, genenames, src, numgenes, qt)}:}{
Initialize an Attractor object

Arguments:

  \code{id} : ID of the attractor

  \code{a} : A vector of MIs of the genes in the attractor. As output by the findAttractor function

  \code{genenames} : A vector of all gene symboles in the dataset.

  \code{src} : The source of the attractor, usually the name of the dataset it is generated.

  \code{numgenes} : Number of top genes stored in the object.

  \code{qt} : Which rank should be used to represent the strength of the attractor. Default is the 10th highest MI

    }
  }
}

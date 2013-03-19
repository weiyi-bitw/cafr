\name{cafr-package}
\alias{cafr-package}
\alias{cafr}
\docType{package}
\title{
Attractor Metagenes Finding Algorithm
}
\description{
Mining gene expression profiles has proven valuable for identifying signatures serving as surrogates of cancer phenotypes. However, the similarities of such signatures across different cancer types have not been strong enough to conclude that they represent a universal biological mechanism shared among multiple cancer types. The Attractor Metagenes Finding Algorithm is a computational method for generating signatures using an iterative process that converges to one of several precise attractors defining signatures representing biomolecular events, such as cell transdifferentiation or the presence of an amplicon.
}
\details{
\tabular{ll}{
Package: \tab cafr\cr
Version: \tab 0.3\cr
Date: \tab 2013-03-18\cr
Depends: \tab R (>= 2.13)\cr
License: \tab GPL (>= 2)\cr
Built: \tab R 2.15.2; x86_64-unknown-linux-gnu; 2013-03-18 21:16:10 UTC; unix\cr
}

Index:
\preformatted{
Attractor               Class '"Attractor"'
AttractorSet            Class '"AttractorSet"'
attractorScanning       Find all attractors in the dataset
attractorScanningGL     Find all genomically-localized attractors in
                        the dataset
clusterAttractors       Cluster attractors from difference datasets
createMetageneSpace     Create metagene-level expression matrix using a
                        given list
findAttractor           Finding attractor using the seed gene
findGLAttractor         Finding genomically-localized attractor using
                        the seed gene
getCorr                 Functions to calculate correlation coefficient
getMI                   Functions to calculate mutual information
loadExpr                Load txt files into matrix or dataframe
parAttractorScanning    Parallelized attractor scanning function
probeSummarization      Gene-level expression summarization
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Cheng, Wei-Yi <wc2302@columbia.edu>

Maintainer: Cheng, Wei-Yi <wc2302@columbia.edu>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}

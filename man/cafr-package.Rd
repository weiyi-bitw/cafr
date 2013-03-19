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
}
\author{
Cheng, Wei-Yi <wc2302@columbia.edu>
}
\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
}
\keyword{ package }

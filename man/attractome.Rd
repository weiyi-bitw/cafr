\name{attractome}
\alias{attractome}
\docType{data}
\title{
List of attractors appears in several different cancer types.
}
\description{
The object is a list of attractors obtained from different cancer types 
(breast, colon and ovarian). The list was derived by generating the attractor
metagenes in 6 independant datasets and then compiling them into consensus 
ranking after calculating the average MI of each gene across datasets.
}
\usage{data(attractome)}
\format{
A list of 12 attractor metagenes. Each attractor metagene is represented by
a matrix with first column being the genes ranked by their average MI across
6 datasets.
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
}
\keyword{datasets}

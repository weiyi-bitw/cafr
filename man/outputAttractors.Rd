\name{outputAttractors}
\alias{outputAttractors}
\title{Output a table for attractors in the dataset}
\description{Create a human-readable table list the top genes for every attractor in the dataset and write to file.}


\usage{
outputAttractors(x, min.basins=2, strength.pos=10, outputGeneNumber=100, write2File=FALSE, fileName="attractors.txt")
}

\arguments{
  \item{x}{A numerical matrix created by the \code{attractorScanning} function.}
  \item{min.basins}{Minimal number of basins (attractees) for an attractor to be considered meaningful.}
  \item{strength.pos}{The rank position of the gene list to be considered as the strength of the attractor..}
  \item{outputGeneNumber}{Number of genes in each attractor to be output.}
  \item{write2File}{Whether write the output matrix to a file.}
  \item{fileName}{The file name if write to file.}
}

\details{
  The function output every attractor found by the \code{attractorScanning} function, excluding any 
  attractor that have less than \code{min.basin} basin genes (or as we call, 'attractees'). Each 
  attractor is represented by the top \code{outputGeneNumber} genes, ranked by the gene's mutual 
  information with the converged metagene.
}

\value{
  A character matrix in which each attractor is represented by two columns: first the top genes, followed
  by each gene's mutual information with the converged metagene.
}

\examples{
\dontrun{
# Load the toy dataset extracted from TCGA OV data
data(ov)

# find all attractors in the dataset
as.ov <- attractorScanning(ov)

mat <- outputAttractors(as.ov, outputGeneNumber=20)
}
}

\seealso{
  \code{\link{attractorScanning}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI, cluster}

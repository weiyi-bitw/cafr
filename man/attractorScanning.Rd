\name{attractorScanning}
\alias{attractorScanning}
\title{Find all attractors in the dataset}
\description{Exhaustively search for all attractors in the dataset.}


\usage{
attractorScanning(data, a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE)
}

\arguments{
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{a}{Exponent of the mutual information, used to create weight vector for metagenes. }
  \item{maxIter}{Max number of iterations.}
  \item{epsilon}{Threshold of convergence.}
  \item{bin}{Number of bins used when estimate mutual information (default=6).}
  \item{so}{Spline order used when estimate mutual information (default=3).}
  \item{rankBased}{When \code{TRUE}, convert the expression values into ranks.}
  \item{negateMI}{When \code{TRUE}, negate the mutual information if the two vectors have negative 
                  momentum.}
}

\details{
  \code{attractorScanning} searches for all the attractors in the dataset by applying \code{CAFrun} 
  using every gene as a seed in the dataset. It will ignore the attractors which were dominated by
  the seed gene. As the function finds a new attractor, it removes the genes from the search list
  that are more significant (with higher MI) than the seed, due to that they will lead to the same
  attractor. Therefore the size of the search list decreases exponentially during the process.
  
  During the search process, the program shows its progress by the lines:

  CENPA ( 10 / 300 ) ... 

  which shows the current seed (CENPA), the number of attractor it tries to finds (10th), the number
  of seeds left to be searched (300). 
}

\value{
  Returns a matrix of size k by m, where m is the number of genes (rows) in the dataset, and k the
  number of converged attractors. Each row of the matrix is the MI of all the genes with the converged
  attractor. 
}

\note{ Missing values are not allowed as the input to the function in the current version.}
\examples{
# Load the toy dataset extracted from TCGA OV data
data(ov)

# find attractor using CENPA as a seed
as <- attractorScanning(ov)

# display the top 20 genes in first attractor
sort(as[1,], decreasing=T)[1:20]

}

\seealso{
  \code{\link{parAttractorScanning}} for excuting \code{attractorScanning} on SGE framework. 
  \code{\link{findAttractor}}, \code{\link{attractorScanningGL}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI}

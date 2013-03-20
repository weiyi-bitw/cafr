\name{parAttractorScanning}
\alias{parAttractorScanning}
\title{Parallelized attractor scanning function}
\description{Run attractor scanning on a list of seed genes.}


\usage{
parAttractorScanning(data, taskList=list(1:nrow(data)), wid=1,  a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE)
}

\arguments{
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{taskList}{List of vectors contains indices of seed genes to be run by each worker.}
  \item{wid}{Worker ID to be assinged the tasks.}
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
  \code{parAttractorScanning} performs \code{CAFrun} on a subset of seeds in the dataset. The
  function takes the tasks as assigned to its worker ID. When applied on a parallel framework 
  such as Sun Grid Engine, it becomes a parallelized version of \code{attractorScanning}.
  The function also reduces the size of seed list in the same way as \code{attractorScanning}
  does, so it is recommended that each worker should be given at least 100 seed genes for efficiency.

  The scripts for running parallelized attractor scanning are given in \code{"cafr/parallel"} see \code{"README.txt"} for more information.

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

# Using the createFolds in the caret package to allocate tasks, assuming
# there are 100 workers
require(caret)
taskList <- createFolds(1:nrow(ov), k=100)

# Take the tasks given to worker 1
as <- parAttractorScanning(ov, taskList=taskList, wid=1)

# display the top 20 genes in first attractor
sort(as[1,], decreasing=TRUE)[1:20]

}

\seealso{
  \code{\link{findAttractor}},\code{\link{parAttractorScanning}}, \code{\link{attractorScanningGL}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI, parallel}

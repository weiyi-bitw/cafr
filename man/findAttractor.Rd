\name{findAttractor}
\alias{findAttractor}
\alias{findCorrAttractor}
\title{Finding attractor using the seed gene}
\description{Find a converged attractor based on the seed gene provided.}


\usage{
findAttractor(data, vec, a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI = TRUE, verbose=TRUE)
findCorrAttractor(data, vec, a=10, maxIter=100, epsilon=1E-14, rankBased=FALSE)
}

\arguments{
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{vec}{A vector of expression values of the seed gene.}
  \item{a}{Exponent of the mutual information, used to create weight vector for metagenes. }
  \item{maxIter}{Max number of iterations.}
  \item{epsilon}{Threshold of convergence.}
  \item{bin}{Number of bins used when estimate mutual information (default=6).}
  \item{so}{Spline order used when estimate mutual information (default=3).}
  \item{rankBased}{When \code{TRUE}, convert the expression values into ranks.}
  \item{negateMI}{When \code{TRUE}, negate the mutual information if the two vectors have negative 
                  momentum.}
  \item{verbose}{When \code{TRUE}, it will show the top 20 genes of the metagene in each iteration.}
}

\details{
  Given a seed gene, \code{findAttractor} finds the attractor by iteratively creating metagenes using 
  the weighted average of all genes in \code{data}. The weight of each gene is determined by 
  \eqn{(MI)^a} where \eqn{MI} is the mutual information between the gene and the metagene. The metagene 
  in the initial iteration is the seed itself. The iterative process continues until convergence is 
  reached. In each iteration, the program calculates the sum of squares \eqn{SS} of the differences 
  between the vector of MIs in the current iteration and that of MIs in the previous iteration:

  \eqn{SS = (I_1 - I'_1)^2 + (I_2 - I'_2)^2 + \ldots}

  If \eqn{SS <} \code{epsilon}, the process is deemed converged. If the process cannot reach 
  convergence within \code{maxIter}, \code{NULL} is returned. If \code{negateMI} is \code{TRUE}, MI 
  between the genes with negative momentum with metagene will be negated, and the weight of this 
  gene is set to 0.

  \code{findCorrAttractor} does the same thing as \code{findAttractor}, but using correlation 
  coefficients instead of MI as association measure.
}

\value{
  Returns a vector of MIs between all the genes in the dataset and the converged metagene. If the 
  process does not converge within the maximum number of iteration assigned, it will return \code{NULL}.
}

\note{ Missing values are not allowed as the input to the function in the current version.}
\examples{
# Load the toy dataset extracted from TCGA OV data
data(ov)

# find attractor using CENPA as a seed
x <- ov["CENPA",]
a <- findAttractor(ov, x)

# display the top 10 genes in the attractor
sort(a, decreasing=TRUE)[1:10]

}

\seealso{
  \code{\link{getCorr}}, \code{\link{getMI}},\code{\link{probeSummarization}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI}

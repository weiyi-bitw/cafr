\name{findGLAttractor}
\alias{findGLAttractor}
\title{Finding genomically-localized attractor using the seed gene}
\description{Find a converged genomically-localized attractor based on the seed gene provided.}


\usage{
findGLAttractor(data,seed, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE)
}

\arguments{
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{seed}{A vector of expression values of the seed gene.}
  \item{alpha}{A vector of exponents to try for convergence.}
  \item{windowSize}{Number of genes to be filtered in, centered at the seed gene.}
  \item{maxIter}{Max number of iterations.}
  \item{epsilon}{Threshold of convergence.}
  \item{bin}{Number of bins used when estimate mutual information (default=6).}
  \item{so}{Spline order used when estimate mutual information (default=3).}
  \item{score.position}{Which rank in the attractor will be used as the strength. (default=5)} 
  \item{num.output}{Number of genes to include in the attractor.}
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

  \code{CAFrun} is the same as \code{findAttractor}, but uses C for the whole process (
  \code{findAttractor} only calls C when calculate MI). 
}

\value{
  Returns a vector of MIs between all the genes in the dataset and the converged metagene. If the 
  process does not converge within the maximum number of iteration assigned, it will return \code{NULL}.
}

\note{ Missing values are not allowed as the input to the function in the current version. \code{CAFrun}
 does not have significant advantage than \code{findAttractor} in terms of speed, thus could be obsolete
 in future release.}
\examples{
# Load the toy dataset from TCGA OV data
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

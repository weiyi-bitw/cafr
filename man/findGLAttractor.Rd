\name{findGLAttractor}
\alias{findGLAttractor}
\title{Finding genomically-localized attractor using the seed gene}
\description{Find a converged genomically-localized attractor based on the seed gene provided.}


\usage{
findGLAttractor(data,seed, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE)
}

\arguments{
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{seed}{The gene symbol of the seed gene.}
  \item{genome}{A matrix with genome information. The rownames of the matrix must
               contains gene symbols in their genomical order. See \code{data(genome)}
               for an example.}
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
  The genomically-localized attractor usually points to an amplicon or delecon that contains several
  genes being amplified or deleted simultaneously. Such amplification and deletion within the region
  reflects on the co-expression of genes in the neighborhood. \code{findGLAttractor} performs a 
  mini-scaled \code{attractorScanningGL} to find the most prominent genomically-localized attractor 
  in the region of the seed gene. It uses the same algorithm as \code{findAttractor}, but restricts 
  the gene space to the neighborhood of the seed gene defined by the \code{windowSize}. The 
  significance (or strength) of an attractor is defined by the MI of the \code{score.position}-th 
  gene. Default is the fifth-highest MI in the attractor. The function tries all the possible 
  exponents defined in \code{alpha} and finds the best one that gives highest 
  \code{score.position}-th MI. It then compares all the converged attractors using the seeds in the 
  neighborhood and finds the one with the highest strength.
}

\value{
  Returns a vector of MIs of the top genes with the converged metagene. If the process does not 
  converge within the maximum number of iteration assigned, it will return \code{NULL}. Each column name 
  of the matrix is the seed used to generate the attractor. 
}
\examples{
# Load the toy dataset extracted from TCGA OV data
data(ov)

# Load the genome information
data(genome)

# find a genomically localized attractor in the neighborhood of PUF60 (8q24.3)
seed <- "PUF60"
a <- findGLAttractor(ov, seed, genome)

# display the top 10 genes in the attractor
print(a)

}

\seealso{
  \code{\link{findAttractor}}, \code{\link{attractorScanning}},\code{\link{attractorScanningGL}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI, amplicon, delecon}

\name{attractorScanningGL}
\alias{attractorScanningGL}
\title{Find all genomically-localized attractors in the dataset}
\description{Exhaustively search for all converged genomically-localized attractors in the dataset.}


\usage{
attractorScanningGL(data, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE, saveAllScore=FALSE)
}

\arguments{
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
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
  \item{verbose}{When \code{TRUE}, it shows the progress.}
  \item{saveAllScore}{When \code{TRUE}, the function output every converged attractor using each seed. It could take a lot of memory.}
}

\details{
  The genomically-localized attractor usually points to an amplicon or delecon that contains several
  genes being amplified or deleted simultaneously. Such amplification and deletion within the region
  reflects on the co-expression of genes in the neighborhood. \code{attractorScanningGL} exhaustively 
  uses every gene as a seed in the dataset to find all converged genomically locallized attractors
  in the region of the seed gene. It uses the same algorithm as \code{findAttractor}, but restricts 
  the gene space to the neighborhood of the seed gene defined by the \code{windowSize}. The 
  significance (or strength) of an attractor is defined by the MI of the \code{score.position}-th 
  gene. Default is the fifth-highest MI in the attractor. The function tries all the possible 
  exponents defined in \code{alpha} and finds the best one that gives highest 
  \code{score.position}-th MI. It then compares all the converged attractors using the seeds in the 
  neighborhood and finds the one with the highest strength.
}

\value{
  When \code{saveAllScore} is \code{TRUE}, it returns a list with the following fields:
 
  \item{attractome}{A matrix of all converged attractors. Each column is the genes of the attractor ranked by their MI with the metagene. Each column name is the seed used to generate the attractors.} 
  \item{bestAlphas}{The best exponents that lead to the highest strength given the seed.}
  \item{score}{The strength of each attractor.}
  \item{scoremat}{A matrix of the MIs corresponding to each gene in the \code{attractome} field.}
  \item{summary}{The summarized output after removing the overlapping attractors. The output only contains the attractor with highest strength among all overlapping ones. Each row is the genes of an attractor ranked according to their MIs with the metagene. Each row name is the seed that leads to the attractor. The first column gives the cytoband information. The last column gives the strength of the attractor.}

  If \code{saveAllScore} is \code{FALSE}, the function only output the \code{summary} field above in matrix form.
}
\examples{
# Load the toy dataset extracted from TCGA OV data
data(ov)

# Load the genome information
data(genome)

# find every genomically localized attractor in the toy dataset
out <- attractorScanningGL(ov, genome)

# display the summarized output
print(out)

}

\seealso{
  \code{\link{findAttractor}}, \code{\link{findGLAttractor}}, \code{\link{attractorScanning}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI, amplicon, delecon}

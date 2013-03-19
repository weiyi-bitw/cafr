\name{getMI}
\alias{getMI}
\alias{getAllMIWz}
\title{Functions to calculate mutual information}
\description{Functions to calculate mutual information using B-spline method.}


\usage{
getMI(x, y, bin=6, so=3, rankBased=FALSE, normalize=TRUE, negateMI=FALSE)
getAllMIWz(data, x, bin=6, so=3, rankBased=FALSE, normalize=TRUE, sorting=FALSE, negateMI=FALSE)
}

\arguments{
  \item{x, y}{A vector of expression values}
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{bin}{Number of bins used when estimate mutual information (default=6).}
  \item{so}{Spline order used when estimate mutual information (default=3).}
  \item{rankBased}{When \code{TRUE}, convert the expression values into ranks.}
  \item{normalize}{When \code{TRUE}, normalize the mutual information by dividing the value by 
    the maximum of the self-mutual information of the two variables \emph{i.e.} \code{max(getMI(x, x), 
    getMI(y, y))}.}
  \item{sorting}{Whether sort the output vector of MI's.}
  \item{negateMI}{When \code{TRUE}, negate the mutual information if the two vectors have negative 
    momentum.}
}

\details{
  \code{getMI} 
  evaluate the mutual information between two vectors using B-spline functions. 
  The function calls a C function implemented by Daub \emph{et al.}. The default setting uses
  six bins and third-order B-splines. \code{getAllMIWz} evaluates the mutual information between 
  every row in \code{data} and the given vector \code{x}. The output vector of mutual information
  will be sorted if \code{sorting} is \code{TRUE}.

  The \code{normalize} option normalizes the mutual information by dividing the value by the maximum
  of the self-mutual information of the two variables. It makes all the MI's between 0 and 1 and thus 
  beneficial to calculation of Attractor Metagenes.
}

\value{
  \code{getMI} returns the mutual information between two vectors. \code{getAllMIWz} returns a vector of 
  mutual information between each row of \code{data} and the vector \code{x}.
}

\note{ Missing values are not allowed as the input to the function in the current version. }
\examples{
# Load the toy dataset extracted from TCGA OV data
data(ov)
#==============
# First example: find the mutual information between CENPA and MELK
#==============
x <- ov["CENPA",]
y <- ov["MELK",]
mi <- getMI(x, y)

#==============
# Second example: evaluate the mutual information between every gene in the dataset and CENPA
#==============
mis <- getAllMIWz(ov, x, sorting=TRUE)

}

\seealso{
  \code{\link{findAttractor}} for finding attractor using mutual information.
  \code{\link{getCorr}}, \code{\link{getAllCorWz}} 
}

\references{
  Carsten O Daub, Ralf Steuer, Joachim Selbig and Sebastian Kloska, 
  Estimating mutual information using B-spline functions – an improved similarity measure for 
  analysing gene expression data. BMC Bioinformatics 2004, 5:118.

  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{MI}
\concept{mutual information, spline, information theory}

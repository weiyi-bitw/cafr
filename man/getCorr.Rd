\name{getCorr}
\alias{getCorr}
\alias{getAllCorWz}
\title{Functions to calculate correlation coefficient}
\description{Functions to calculate Pearson / Spearman correlation coefficient.}


\usage{
getCorr(x, y, rankBased = FALSE)

getAllCorWz(data, x, rankBased=FALSE, sorting=FALSE)
}

\arguments{
  \item{x, y}{A vector of expression values}
  \item{data}{An expression matrix with genes in the rows, samples in the columns.}
  \item{rankBased}{When \code{TRUE}, evaluate the Spearman's rank correlation coefficient.}
  \item{sorting}{Whether sort the output vector of correlation coefficients.}
}

\details{
  \code{getCorr} 
  evaluate the correlation coefficient between two vectors. It gives the same result as the generic
  function \code{cor}.\code{getAllMIWz} evaluates the correlation coefficients between every row in 
  \code{data} and the given vector \code{x}. The output vector of correlation coefficients will be 
  sorted if \code{sorting} is \code{TRUE}.

  Spearman's rank correlation coefficients are calculated if \code{rankBased} is set to \code{TRUE}.
}

\value{
  \code{getCorr} returns the correlation coefficient between two vectors. 
  \code{getAllMIWz} returns a vector of correlation coefficients between each row of \code{data} 
  and the vector \code{x}.
}

\note{ Missing values are not allowed as the input to the function in the current version. }
\examples{
# Load the toy dataset from TCGA OV data
data(ov)
#==============
# First example: find the Pearson correlation coefficients between CENPA and MELK
#==============
x <- ov["CENPA",]
y <- ov["MELK",]
mi <- getCorr(x, y)

#==============
# Second example: evaluate the Spearman correlation coefficients between 
# every gene in the dataset and CENPA
#==============
mis <- getAllCorWz(ov, x, rankBased=TRUE, sorting=TRUE)

}

\seealso{
  \code{\link{cor}} for the generic functions to calculate correlation coefficients.
  \code{\link{findCorrAttractor}} for finding the correlation coefficient based attractors.
  \code{\link{getMI}},  \code{\link{getAllMIWz}} 
}

\author{Wei-Yi Cheng}

\keyword{Corr}
\concept{correlation, Pearson, Spearman}

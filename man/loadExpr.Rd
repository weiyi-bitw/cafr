\name{loadExpr}
\alias{loadExpr}
\alias{loadClin}
\title{Load txt files into matrix or dataframe}
\description{Functions to load txt files of gene expression data or clinical data into matrix or dataframe.}


\usage{
  loadExpr(file, sep='\t')
  loadClin(file, sep='\t')  
}

\arguments{
  \item{file}{File path to the text file, in which first row is the header, first column contains row names.}
  \item{sep}{the field separator character. As defined in the \code{read.table} function.}
}

\details{
  \code{loadExpr} and \code{loadClin} use \code{readLines} function to process text files, thus bypass
  some cases where \code{read.table} hangs due to warnings caused by the file. It also provides a 
  progress bar for progress monitoring.
}

\value{
  \code{loadExpr} returns a matrix with expression values. \code{loadClin} returns a matrix with
  text entries. The popular notations for missing values are treated as 'NA' in the matrix.
}

\examples{
# This example requires a Synapse account. It shows how to load expression data and 
# clinical data from Synapse.
\dontrun{
require(synapseClient)
synapseLogin()

# Load TCGA READ IlluminaHiSeq RNASeq data
syn <- loadEntity("syn1446276")
fileName <- file.path(syn$cacheDir, syn$files[[1]])
ge <- loadExpr(fileName)

# Load TCGA READ Clinical_patient file
syn <- loadEntity("syn1446151")
fileName <- file.path(syn$cacheDir, syn$files[[1]])
clnc <- loadClin(fileName)
}
}

\seealso{
  \code{\link{read.table}} 
}

\author{Wei-Yi Cheng}

\keyword{read}
\concept{parse files, read tables}

\name{probeSummarization}
\alias{probeSummarization}
\alias{summarizeGenes}
\title{Gene-level expression summarization}
\description{Summarize probe-set-level expression matrix into gene-level expression matrix.}


\usage{
probeSummarization(ge, map, method="corr", threshold=0.5, gene.colname="Gene.Symbol", verbose=TRUE)
summarizeGenes(ge, map, sumfun=median)
}

\arguments{
  \item{ge}{Gene expression matrix with each row as a gene and each column as a sample.}
  \item{map}{Path to gene symbol annotation file with rownames as probe set IDs. Must contains a column 
             with gene symbols.}
  \item{method}{Methods to evaluate the association between probe sets of the same gene.\code{corr} for
               Pearson correlation coefficients. \code{mi} for mutual information.}
  \item{threshold}{Threshold below which the probe set will not be used for summarization. Default = 0.5}
  \item{gene.colname}{The column name in the \code{map} file that contains the gene symbols.}
  \item{verbose}{If \code{TRUE}, show the summarization process.}
  \item{sumfun}{Summarization function used for simple summarization.}
}

\details{
  When running attractor finding program, it is important to summarize probe-set-level expression 
  into gene-level expression. It can prevent the genes with multiple probe sets in the microarray 
  dominate the rank in the attractors. \code{probeSummarization} achieves this by taking mean
  values of probe sets of the same gene while discarding the probe set with significantly different
  expression pattern. These 'bad probe sets' were identified by calculating the association between 
  every the probe set with the sum of the probe sets. If the association is less than the threshold,
  the probe set was discarded. The remaining probe sets were summarized by their mean values.

  \code{summarizeGenes} simply summarize the probe sets by applying the user-defined function 
  \code{sumfun}. The two functions are expected to be integrated in the future release.

}

\value{
  An gene-level expression matrix with genes at the rows and samples at the columns.
}

\note{ Two functions are expected to be integrated in the future release, in which \code{summarizeGenes}
     will be obsolete. }
\examples{
# load Toy version of Wang et al. breast cancer dataset (GSE2034)
data(brca.pbs)

# download the HGU133A 2.0 annotations
source("http://bioconductor.org/biocLite.R")
biocLite("hgu133a2.db")
library(hgu133a2.db)

# Create map object to fit the format
x <- hgu133a2SYMBOL
map <- cbind(unlist(as.list(x[mappedkeys(x)])))
colnames(map) <- "Gene.Symbol"

# summarize into gene-level expression after eliminating uncorrelated probes
brca <- probeSummarization(brca.pbs, map)

# summarize into gene-level expression using median (default)
brca <- summarizeGenes(brca)
}

\seealso{
  \code{\link{findAttractor}},
  \code{\link{createGeneSpace}},
  \code{\link{createMetageneSpace}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{summarization}
\concept{gene-level expression, gene-level summarization}

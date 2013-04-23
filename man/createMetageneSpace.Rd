\name{createMetageneSpace}
\alias{createMetageneSpace}
\title{Create metagene-level expression matrix using a given list}
\description{Summarize expression matrix into metagene-level expression matrix.}


\usage{
createMetageneSpace(ge, attractome, map=NULL, rownamesMap=FALSE, chosenProbes=NULL, gene.colname="Gene.Symbol")
}

\arguments{
  \item{ge}{Gene expression matrix with each row as a gene and each column as a sample.}
  \item{attractome}{A list of metagenes and their gene mambers. Each metagene in the list is a matrix
                   with gene symbols in the first column. See the format of \code{data(attractome)}. }
  \item{map}{Path to gene symbol annotation file with rownames as probe set IDs. Must contains a column 
             with gene symbols. If \code{NULL}, the program will assume the dataset is already a 
             gene-level expression matrix.}
  \item{rownamesMap}{If \code{TRUE} and no provided \code{map}, use row names of the dataset directly to create metagenes.}
  \item{chosenProbes}{A list of probe set IDs selected by previous run. If this argument is assigned, 
                     the program will not find the best correlated probe sets for summarization. This 
                     is used for probe set consistency among different datasets.}
  \item{gene.colname}{The column name in the \code{map} file that contains the gene symbols.}
}

\details{
	Given the gene members of the attractors, this function transforms the probe-level or gene-level
        expression matrix into the metagene-level expression matrix. One can use the metagene expression 
        matrix to build predictive models, observe relations between biomolecular events represented by 
        Attractor Metagenes, create figures, etc.. 

        If \code{ge} is a probe-level expression matrix, the function first summarize the probe-level 
        expression into gene-level expression by taking mean values after discarding the un-correlatted
        probe sets (similar to the \code{probeSummarization} function). Then it create metagenes using
        the mean of the gene-level expression of the genes. If the gene members in a metagene is more 
        than ten, it will only take the top 10 genes in the list. If \code{ge} is gene-level expression
        matrix, the function will directly take the mean of the gene members and create the metagene.

        For consistency of the probe set used, one can use the \code{pbs} field of the output from 
        previous run to create metagene space using the same selected probe sets in the previous run.
}

\value{
  When \code{map} is provided, the function gives a list containing the following field:

  \item{metaSpace}{the metagene-level expression matrix.}
  \item{pbs}{the probe set IDs used to summarize each gene in the metagenes.}

  When \code{map} is \code{NULL}, the function gives a matrix of the metagene-level expression matrix.
}

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

# load the attractor list
data(attractome)

# summarize into metagene-level expression
o <- createMetageneSpace(brca.pbs, attractome, map)
meta <- o$metaSpace

# create metagene expression matrix using the selected probes
meta2 <- createMetageneSpace(brca.pbs, attractome, chosenProbes=o$pbs)

}

\seealso{
  \code{\link{findAttractor}},
  \code{\link{probeSummarization}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{summarization}
\concept{gene-level expression, gene-level summarization, metagene, attractor metagenes, metagene expression}

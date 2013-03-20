\name{clusterAttractors}
\alias{clusterAttractors}
\title{Cluster attractors from difference datasets}
\description{Cluster attractors from difference datasets based on their overlapping in the top genes.}


\usage{
clusterAttractors(filePath="./", fileNames, numGenes=100, datasetTags=NULL)
}

\arguments{
  \item{filePath}{A path to the directory where the RData files were saved. Default is the working directory.}
  \item{fileNames}{Names of the RData files to be loaded.}
  \item{numGenes}{Number of genes load into the Attractor object.}
  \item{datasetTags}{A vector of strings of dataset tags. Used to annotate the Attractor objects.}
}

\details{
  \code{clusterAttractors} loads in the attractor matrices from different datasets produced by 
  \code{attractorScanning} and clusters attractors from different datasets based on the number of overlapping
  genes. 
}

\value{
  A list of AttractorSets and Attractors. Each AttractorSet contains the most similar Attractors from different
  datasets. Attractors that cannot be clustered are left as they are. The contents in the list are ranked first
  by how many Attractors are clustered in the set, then by the minimum strength among the Attractors in the set.
}

\examples{
\dontrun{
# Load the toy dataset extracted from TCGA OV data
data(ov)

# find all attractors in the dataset
as.ov <- attractorScanning(ov)
# Save the result as an RData file
save(as.ov, file="./as.ov.rda")

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

as.brca <- attractorScanning(brca)
save(as.brca, file="./as.brca.rda")

# Cluster attractors from these two datasets, because it's a toy dataset
# with much less genes, only load in top 20 genes
alist <- clusterAttractors(fileNames=c("as.ov.rda", "as.brca.rda"), numGenes=20, datasetTags=c("OV", "BRCA"))

# Show the consensus attractor of the strongest attractor cluster
alist[[1]]$getConsensus(sz=10)

# remove the saved files
unlink("./as.ov.rda")
unlink("./as.brca.rda")

}
}

\seealso{
  \code{\link{AttractorSet}}, \code{\link{Attractor}}
}

\references{
  Wei-Yi Cheng, Tai-Hsien Ou Yang and Dimitris Anastassiou, Biomolecular events in cancer revealed 
  by attractor metagenes, PLoS Computational Biology, Vol. 9, Issue 2, February 2013. 
  }

\author{Wei-Yi Cheng}

\keyword{Attractor Metagenes}
\concept{attractor, metagene, MI, cluster}

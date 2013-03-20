\name{brca.pbs}
\alias{brca.pbs}
\docType{data}
\title{
A toy dataset extracted from Wang et al. breast cancer dataset.
}
\description{
This dataset was sampled from the Wang et al. breast cancer dataset (GSE2034). It contains 338 probe sets
and 200 samples randomly sampled from the original dataset. For demonstration purpose of the 
Attractor Metagene methods, the probe sets include the top 10 genes from the CIN, MES and LYM 
attractor. The rest genes were randomly sampled.
}
\usage{data(brca.pbs)}
\format{
  A matrix of expression values with 338 rows (probe sets) and 200 samples (columns).
}
\source{
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034
}
\references{
    Wang Y, Klijn JG, Zhang Y, Sieuwerts AM et al. Gene-expression profiles to predict distant metastasis of lymph-node-negative primary breast cancer. Lancet 2005 Feb 19-25;365(9460):671-9. PMID: 15721472
}
\keyword{datasets}

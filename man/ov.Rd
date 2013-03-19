\name{ov}
\alias{ov}
\docType{data}
\title{
A toy dataset extracted from TCGA OV dataset.
}
\description{
This dataset was sampled from the official TCGA OV dataset. It contains 455 genes
and 200 samples. For demonstration purpose of the Attractor Metagene methods, the
genes include the top 10 genes from the CIN, MES and LYM attractor, and also the 
genes from the universal chr8q24.3 amplicon. The rest genes were randomly sampled.
}
\usage{data(ov)}
\format{
  A matrix of expression values with 455 rows (genes) and 200 samples (columns).
}
\source{
http://cancergenome.nih.gov
}
\references{
 The Cancer Genome Atlas Research Network (2011) Integrated genomic analyses of ovarian carcinoma. Nature. 474(7353):609-615.
}
\keyword{datasets}

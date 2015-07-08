load("/pred_cello1/analysis/RNAseqAnalysis/workspace/rnaseq_voom/rnaseq.voom.LUAD.rdata")
library(cafr)
x <- rnaseq[1:10,]
out <- AllPairwiseCor(x)
savehistory("cafr.test.R")

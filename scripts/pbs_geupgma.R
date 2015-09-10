#
# R script to be submitted for UPGMA execution in Roche HPC
#
# Wei-Yi Cheng
# Scientist, Data Science III
# 2015/05/08
#
# Args:
#	  args[1] : input gene expression data matrix
#	  args[2] : job codename (optional)
#
# Usage: 
# $ Rscript pbs_geupgma.R data.rda JOBID
#

library(cafr)

setwd(Sys.getenv("PBS_O_WORKDIR"))
args <- commandArgs(TRUE)

if (length(args) > 1) {
	codename <- args[2]
	dirname <- codename
}

dir.create(dirname)

env <- new.env()
nm <- load(args[1], env)[1]
ge <- env[[nm]]
corout <- AllPairwiseCor(ge)
trc <- GEUPGMA(corout, ge)
save(trc, file=file.path(dirname, "upgma.trace.rda"))
cl <- ClusterFeatures(ge, trc)
save(cl, file=file.path(dirname, "upgma.cluster.rda"))



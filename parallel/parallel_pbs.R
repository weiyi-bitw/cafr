# R script to be submitted for parallel execution in Roche HPC
#
# Wei-Yi Cheng
# Scientist, Data Science III
# 2015/05/08
#
# Args:
#	  args[1] : input gene expression data matrix
#	  args[2] : total workers
#	  args[3] : job codename (optional)
#
# Usage: PartialFindAttractor on worker 1 in a 10-worker cluster:
# $ Rscript parallel_pbs.R data.rda 10 JOBID
#

require(cafr)


setwd(Sys.getenv("PBS_O_WORKDIR"))
args <- commandArgs(TRUE)

#================================
# Load data
#===============================

env <- new.env()
nm <- load(args[1], env)[1]
ge <- env[[nm]]
# extract gene names if it's in TCGA format
if (grepl("\\|", rownames(ge)[1])) {
  rownames(ge) <- ExtractGeneNameTCGA(rownames(ge))
}

codename = args[1]
wid = 1
numWorkers=1
dirname <- "output"
if (length(args) > 1) {
	wid <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
	numWorkers <- as.numeric(args[2])
}
if (length(args) > 2) {
	codename <- args[3]
	dirname <- codename
}

cat("=======================================\n")
cat("JOB NAME:\t", codename, "\n", sep="")
cat("FILE NAME:\t", args[1], "\n", sep="")
cat("WID:\t", wid, "\n",sep="")
cat("TOTAL WORKERS:\t",  numWorkers, "\n", sep="")
cat("=======================================\n")
flush.console()

m <- nrow(ge)
set.seed(913503)
taskList <- lapply(1:numWorkers, function(f){
    sample(1:m)[(floor(m/numWorkers*(f-1))+1):floor(m/numWorkers*f)]
})

cat("Finding attractors...\n")
flush.console()

out <- PartialAttractorScanning(ge, maxIter=500 , taskList=taskList, wid=wid)
as <- out$attractorMatrix
dd <- out$dominant
dir.create(dirname)
save(as, dd, 
     file=paste(dirname,"/",codename,".attractors.", 
                sprintf("%04d", wid), ".rda", sep="")
)

cat("Done.\n")
cat("=======================================\n")
flush.console()


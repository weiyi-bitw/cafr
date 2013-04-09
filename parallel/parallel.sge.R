# R script for Amazon EC2 cluster
# Wei-Yi Cheng
# 2013.03.07
#
# arguments for script:
#	args[1] : input gene expression data matrix
#	args[2] : worker id
#	args[3] : total workers
#	args[4] : job codename (optional)
# For example, to run parFindAttractor on worker 1 in a 10-worker cluster:
# > Rscript parallel.sge.R data.rda 1 10 JOBID
#
#

require(devtools)
require(caret)
require(RCurl)
require(cafr)

args <- commandArgs(TRUE)

#================================
# Load data
#===============================

env <- new.env()
nm <- load(args[1], env)[1]
ge <- env[[nm]]

codename = args[1]
wid = 1
numWorkers=1
dirname = "output"
if(length(args)>1){
	wid = as.numeric(args[2])
	numWorkers = as.numeric(args[3])
}
if(length(args)>3){
	codename = args[4]
	dirname=codename
}

cat("=======================================\n")
cat("JOB NAME:\t", codename, "\n", sep="")
cat("FILE NAME:\t", args[1], "\n", sep="")
cat("WID:\t", wid, "\n",sep="")
cat("TOTAL WORKERS:\t",  numWorkers, "\n", sep="")
cat("=======================================\n")
flush.console()

m = nrow(ge)
set.seed(913503)
taskList = createFolds(1:m, k=numWorkers)

cat("Finding attractors...\n");flush.console()

as = parAttractorScanning(ge, maxIter=500 , taskList=taskList, wid=wid)
dir.create(dirname)
save(as, file=paste(dirname,"/",codename,".attractors.", sprintf("%04d", wid), ".rda", sep=""))

cat("Done.\n")
cat("=======================================\n");flush.console()


# R script for Amazon EC2 cluster
#
# arguments for script:
#	args[1] : input gene expression data matrix
#	args[2] : worker id
#	args[3] : total workers
#	args[4] : job codename (optional)
#
#

require(devtools)
require(parallel)
require(caret)
require(RCurl)

# Install cafr package
install_github(repo="cafr", username="weiyi-bitw", ref="dev")
library(cafr)

args <- commandArgs(TRUE)

#================================
# Load data
#===============================

load(args[1])
codename = args[1]
wid = 1
numWorkers=1
dirname = "output"
if(length(args)>1){
	wid = as.numeric(args[2])+1
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

as = parAttractorScanning(ge, taskList=taskList, wid=wid)
dir.create(dirname)
save(as, file=paste(dirname,"/",codename,".attractors.", sprintf("%04d", wid), ".rda", sep=""))

cat("Done.\n")
cat("=======================================\n");flush.console()


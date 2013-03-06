# R script for Amazon EC2 cluster
#
# arguments for script:
#	args[1] : input gene expression data matrix
#	args[2] : worker id
#	args[3] : total workers
#
#

require(devtools)
require(parallel)
require(caret)

# Install cafr package
install_github(repo="cafr", username="weiyi-bitw", ref="dev")
library(cafr)

args <- commandArgs(TRUE)

#================================
# Load data
#===============================

load(args[1])
wid = as.numeric(args[2])
numWorkers = as.numeric(args[3])

cat("FILE NAME:\t", args[1], sep="")
cat("WID:\t", wid, sep="")
cat("TOTAL WORKERS:\t", numWorkers, sep="")
flush.console()

m = nrow(ge)
taskList = createFolds(1:m, k=numWorkers)

cat("Finding attractors...\n");flush.console()

as = parAttractorScanning(ge, taskList=taskList, wid=wid)
save(as, file=paste("attractors.", sprintf("%04d", wid), ".rda", sep=""))

cat("Done.\n");flush.console()
